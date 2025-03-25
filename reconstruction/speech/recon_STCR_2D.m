clear; clc; close all;

DATA_PATH = '../../data/speech/';
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
    gpuDevice(1);
end

% File paths
file_paths = dir(fullfile([ISMRMRDPATH, '*.h5']));
nfile = length(file_paths);
USE_GPU = 1;

for file_idx = [1:nfile]
     fprintf('File %d of %d\n', file_idx, length(file_paths));
     
     clearvars -except FOV ISMRMRDPATH NOISEPATH TRAJ_SEARCH_PATH file_paths file_idx USE_GPU

    ismrmrdfile = file_paths(file_idx).name;
    nfile = [NOISEPATH, 'noise_', file_paths(file_idx).name];

    dset = ismrmrd.Dataset([ISMRMRDPATH,ismrmrdfile]);
    header = ismrmrd.xml.deserialize(dset.readxml());
    raw_data    = dset.readAcquisition;

    disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
    disp(['Reconstructing: ' header.measurementInformation.protocolName]); disp(' ');
    disp('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'); disp(' ');
    
    dmtx = noise_adjust(nfile, header, raw_data.head.sample_time_us(1)*1e-6);

    fprintf('Retreiving k-space trajectory...');
    traj_path = fullfile(TRAJ_SEARCH_PATH, [header.userParameters.userParameterString(2).value, '.mat']);
    load(traj_path);
    fprintf('Done.\n');

    spatial_res = [param.spatial_resolution, param.spatial_resolution,8];
    matrix_size = param.matrix_size;
    matrix_size = floor(double(matrix_size)/2) *2;

    % accumulate data
    for j = 1:length(raw_data.data)
        if isempty(raw_data.data{j})
            len_raw_data = j-1;
            break
        end
        data(:,:,j) = raw_data.data{j};
        len_raw_data = j;
    end
    
    % apply noise decorrelation matrix.
    data = data(:,:,1:len_raw_data);
    data = ismrm_apply_noise_decorrelation_mtx(permute(data,[1,3,2]), dmtx);
    data = permute(data, [1, 3, 2]);
    
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
    
    kx = k_pred(:,:,1);
    ky = k_pred(:,:,2);
   
    % discard data
    data = data(param.pre_discard+1:end,:, :);
    data = single(data);
    nrep = floor(size(data,3) / size(kx,2));
    
    nrep = min(10, nrep); % PK hack for max 10 reps
    
    % repeat kx, ky, kz as needed.
    kx = repmat(kx, [1, nrep]);
    ky = repmat(ky, [1, nrep]);
    data = data(:,:,1:size(kx,2));


    % Trim TRS.
    TRToTrim = 0;
    data = data(:,:,TRToTrim+1:end);
    kx = kx(:,TRToTrim+1:end);
    ky = ky(:,TRToTrim+1:end);
    

    % scale kx ky for rthawk
    kmax = max(abs(kx(:) + 1i * ky(:)));
    kx =  0.5 .* (kx / kmax);
    ky = 0.5 .* (ky / kmax);
    
    % construct kx, ky
    n_TRs = size(kx, 2);
    n_arms_per_frame = 7;
    n_read = size(kx, 1);
    
    n_frame = floor(n_TRs / n_arms_per_frame);
    n_coil = size(data, 2);

    kx = reshape(kx, n_read, n_arms_per_frame, n_frame);
    ky = reshape(ky, n_read, n_arms_per_frame, n_frame);
    kspace = reshape(data, n_read, n_coil, n_arms_per_frame, n_frame);
  
    kspace = kspace * 1000;
    kspace = permute(kspace, [1, 3, 4, 2]);
    
    % Encoding operator (using sqrt(w)).
    w = w + eps; % fixing dcf (numerically?)
    rootw = sqrt(w');
    kspace = kspace .* rootw;
    kspace = gpuArray(kspace);
    E = Fnufft_2D(kx, ky, n_coil, matrix_size, USE_GPU, rootw, 1.5, [4, 4]);
    x0_ = E' * kspace;
    
    % sensitivity operator.
    sens = get_sens_map(x0_, '2D');
    C = C_2D(size(x0_), sens, USE_GPU);
    
    % --------- adjoint test on the operator C (optional). --------------------
    test_fatrix_adjoint(C);

    %% First Estimate to the solver, gridding + coil combination
    first_estimate = C' * x0_;
    scale = max(abs(first_estimate(:)));

    %% STCR (functionized)
    tTV = 0.01;
    sTV = 0.001;
    
    [xf, cost] = STCR_2D_NCG_TTV_STV(E*C, first_estimate, kspace, tTV, sTV);
    
    xf = phaseSensitiveIR(xf, 10);
    image_stcr = rectify_pk(xf);
    
    mkdirIfNotExist('recon');
    save(['recon/', file_paths(file_idx).name(1:end-3), '.mat'], 'image_stcr', 'cost', 'header');

end