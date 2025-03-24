function [dmtx, n_stats] = noise_adjust(nfile, i_MRD_h,  data_samp_time)

dsetin      = ismrmrd.Dataset(nfile, 'dataset');
n_MRD_h     = ismrmrd.xml.deserialize(dsetin.readxml);

disp(['Required Sens Map: ' i_MRD_h.measurementInformation.measurementDependency.measurementID ', Noise ID: ' n_MRD_h.measurementInformation.measurementID]);

%% coil check

n_channels = length(i_MRD_h.acquisitionSystemInformation.coilLabel);

coil_label = cell(n_channels,4);
coil_flag = zeros(n_channels,2);

for ii = 1:n_channels
    
    coil_label{ii,1} = i_MRD_h.acquisitionSystemInformation.coilLabel(ii).coilNumber;
    coil_label{ii,2} = i_MRD_h.acquisitionSystemInformation.coilLabel(ii).coilName;
    
    coil_label{ii,3} = n_MRD_h.acquisitionSystemInformation.coilLabel(ii).coilNumber;
    coil_label{ii,4} = n_MRD_h.acquisitionSystemInformation.coilLabel(ii).coilName;
    
    coil_flag(ii,1) = coil_label{ii,1} == coil_label{ii,3};
    coil_flag(ii,2) = strcmp(coil_label{ii,2}, coil_label{ii,4});
end

Data_CoilNum = coil_label(:,1);
Data_CoilName = coil_label(:,2);
Noise_CoilNum = coil_label(:,3);
Noise_CoilName = coil_label(:,4);

disp(' '); disp(table(Data_CoilNum, Data_CoilName, Noise_CoilNum, Noise_CoilName)); disp(' ');

if sum(coil_flag(:,1)) ~= n_channels
    disp('### WARNING ### ! Coil order mismatch! (not critical?) '); disp(' ');
elseif sum(coil_flag(:,2)) ~= n_channels
    disp('### WARNING ###  ! Coil name mismatch!'); disp(' ');
end

%% acquire noise

% assuming Siemens using 2 averages of 128 readouts:
noise_data  = dsetin.readAcquisition(1,256);
% nhlbi_toolbox.plot_experiment(noise_data)

% readout bw filter factor
Siemens_rBW = n_MRD_h.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;

% noise scan info
n_samples   = double(noise_data.head.number_of_samples(1));
n_channels  = double(noise_data.head.active_channels(1));
noise_ind   = length(noise_data.data);

nt2 = zeros(n_samples, noise_ind, n_channels);
for i = 1:noise_ind
    nt2(:,i,:) =  double( noise_data.data{i});
end

% Plot coil noise
if nargout == 2
    nt3 = reshape(nt2,size(nt2,1)*size(nt2,2), size(nt2,3));
    n_stats = std(nt3,[],1);
    
    figure, hold on, title('Coil Noise SD (m +/- std)'); xlabel('Coils'); ylabel('Noise SD');
    plot(1:n_channels, n_stats, 'bo');
    plot([ones(n_channels,1).*mean(n_stats)], 'r-');
    plot([ones(n_channels,1).*(mean(n_stats)-std(n_stats)) ones(n_channels,1).*(mean(n_stats)+std(n_stats))], 'r--');
    
    str = {['Mean ' num2str(mean(n_stats))] ['Std  ' num2str(std(n_stats))]};dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
end

% calculate whitening matrix
n_scaling   = Siemens_rBW * data_samp_time / (noise_data.head.sample_time_us(1)*1e-6);
dmtx        = ismrm_calculate_noise_decorrelation_mtx(nt2, n_scaling ); % figure,imagesc(abs(dmtx));

end