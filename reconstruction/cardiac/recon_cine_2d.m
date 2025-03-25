clear; clc; close all;

DATA_PATH = '../../data/cardiac/';
TRAJ_SEARCH_PATH = DATA_PATH;
USC_RECON_PATH = '../usc_dynamic_reconstruction/';
ISMRMRDPATH = [DATA_PATH, 'h5/'];
NOISEPATH = [DATA_PATH, 'noise/'];

addpath("../ismrmrd/matlab/")
addpath("../ismrm_sunrise_matlab/")
addpath("../girf/")
addpath([USC_RECON_PATH, 'encoding/']);
addpath([USC_RECON_PATH, 'optim/']);
addpath([USC_RECON_PATH, 'utility/']);
addpath(USC_RECON_PATH);

% options   
PHASE_SENSITIVE_IR = 1;
USE_GPU = 0;

if USE_GPU
    gpuDevice(1)
end

file_paths = dir(fullfile([ISMRMRDPATH, '*.h5']));

for file_idx = 1:length(file_paths)

    fprintf('File %d of %d\n', file_idx, length(file_paths));
     
     clearvars -except PHASE_SENSITIVE_IR NOISEPATH ISMRMRDPATH TRAJ_SEARCH_PATH file_paths file_idx USE_GPU
    
    ismrmrdfile = file_paths(file_idx).name;
    nfile = [NOISEPATH, 'noise_', file_paths(file_idx).name];

    dset = ismrmrd.Dataset([ISMRMRDPATH,ismrmrdfile]);
    noise_dset = ismrmrd.Dataset(nfile);

    header = ismrmrd.xml.deserialize(dset.readxml());

    disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
    disp(['Reconstructing: ' header.measurementInformation.protocolName]); disp(' ');
    disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
    
    fprintf('Retreiving k-space trajectory...');
    traj_path = fullfile(TRAJ_SEARCH_PATH, [header.userParameters.userParameterString(2).value, '.mat']);
    load(traj_path);
    fprintf('Done.\n');
    
    %Let's load the data
    slide_size = param.interleaves;
    raw_data = dset.readAcquisition();  % read all the acquisitions
    n_arms_fully_sample = param.interleaves;
    
    dmtx = noise_adjust(nfile, header, raw_data.head.sample_time_us(1)*1e-6);

    %% set up frame indices

    disp(' ');disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
    disp('### Cardiac Binning ###');
    disp(' ');disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
    
    %% GIRF corrections
    % set-up GIRF object
    R = [raw_data.head.read_dir(:,1), raw_data.head.phase_dir(:,1), raw_data.head.slice_dir(:,1)  ]; % Rotation matrix
    tRR = 0; % legacy object for sub-dwelltime delay correction
    sR.R = R; % pack rotation matrix and field info in to one object for backwards compatability..
    sR.T = header.acquisitionSystemInformation.systemFieldStrength_T;

    % identify magnet for specifix GIRF
    sR.systemVendor = header.acquisitionSystemInformation.systemVendor;
    sR.systemModel  = header.acquisitionSystemInformation.systemModel;
    dt = param.dt;
    gx = [zeros([1, size(kx, 2)]); diff(kx)] * 1e3 / dt / 42.58e6;
    gy = [zeros([1, size(ky, 2)]); diff(ky)] * 1e3 / dt / 42.58e6;
    
    g_nom = cat(3, gx, gy);
    g_nom(:, :, 3) = 0;
    
    k_pred = apply_GIRF(g_nom, dt, sR, tRR );
    kx_girf = k_pred(:,:,1);
    ky_girf = k_pred(:,:,2);
    
    matrix_size = double([ceil(param.fov*10 / param.spatial_resolution), ...
                   ceil(param.fov*10 / param.spatial_resolution)]);

    % FOV oversampling
    matrix_size = floor(matrix_size * 1.5);

    % ensure matrix size even.
    matrix_size(1) = 2*ceil(matrix_size(1)/2);
    matrix_size(2) = 2*ceil(matrix_size(2)/2);

    kmax = max(abs(kx_girf(:) + 1i * ky_girf(:)));
    kx_girf = 0.5 .* (kx_girf / kmax);
    ky_girf = 0.5 .* (ky_girf / kmax);
    kx = kx_girf;
    ky = ky_girf;

    n_cardiac_frame=size(kx,2) / param.interleaves;
    ns = size(kx, 1);
    n_coil = size(raw_data.data{1}, 2);

    kxx = zeros(ns, param.interleaves, n_cardiac_frame);
    kyy = zeros(ns, param.interleaves, n_cardiac_frame);
    datadata = zeros(ns, param.interleaves, n_cardiac_frame, n_coil);

    data = cat(3, raw_data.data{:});
    data = data(param.pre_discard+1:end, :, :, :);

    for tr_idx = 1:size(data, 3)
        encodestep = raw_data.head.idx.kspace_encode_step_1(tr_idx)+1;
        phase = mod(floor((tr_idx-1)/8), n_cardiac_frame)+1;

        tr_idx_mod = mod(tr_idx-1,size(kx,2))+1;

        kxx(:,encodestep , phase) = kx(:, tr_idx_mod);
        kyy(:, encodestep, phase) = ky(:,tr_idx_mod);
        datadata(:, encodestep, phase, :) = squeeze(datadata(:, encodestep, phase, :)) + data(:,:,tr_idx);
    end

    kx = kxx;
    ky = kyy;
    data = datadata;

    data = ismrm_apply_noise_decorrelation_mtx(data, dmtx);
    
    data = data .* sqrt(w');
    F = Fnufft_2D(kx, ky, n_coil, [matrix_size(1), matrix_size(1)], USE_GPU, w', 1, [3,3]);
    output_images = F' * data;

    fprintf('Estimating Sensitivity Map\n');
    % get sensitivity maps
    if USE_GPU
        sens = gpuArray(get_sens_map(gather(output_images), '2D'));
    else
        sens = get_sens_map(output_images, '2D');
    end
    C = C_2D(size(output_images), sens, USE_GPU);

    E = F * C;
    image_cc = E' * data;

     opt.Nmaxiter = 20;
     opt.Nlineiter = 20;
     [image_stcr, cost] = STCR_2D_NCG_TTV_STV(E, image_cc, data, 0.01, 0, opt);

    image_cc = gather(image_cc);
    image_stcr = gather(image_stcr);

    % phase sensitive inversion recovery (only for grid)
    if PHASE_SENSITIVE_IR
             [stcr_ir, ~, ir_window] = phaseSensitiveIR(image_stcr, 20);
             image_stcr = stcr_ir;

             [cc_ir, ~, ir_window] = phaseSensitiveIR(image_cc, 20);
            image_cc = cc_ir;            
    end

    %% export to video
    filename_out = [ismrmrdfile(1:end-3)];
    image_cc = rectify_pk(flipud(image_cc));
    image_stcr = rectify_pk(flipud(image_stcr));
    mkdirIfNotExist('recon');
    save(['recon/', filename_out, '_STCR.mat'], 'image_cc', 'image_stcr', 'param', 'header');
end