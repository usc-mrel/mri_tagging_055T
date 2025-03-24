function [wfm, t, wfm_header] = get_wfm_from_rth_header(rth_header, t_points, upsampling)
    
    dt = t_points(2) - t_points(1);
    dt = dt/upsampling;
    t = 0:dt:dt*(length(t_points)*upsampling-1);
    
    if isfield(rth_header, 'user_ADCstart')        
        adc_start_time = rth_header.user_ADCstart;
    else
        error ('header file does not contain "user_ADCstart" field');
    end
     
    wfm_header = get_wfm_header_from_rth_header(rth_header);
    
    wfm = zeros(length(wfm_header), length(t));
    wfm_baseline_t = zeros(length(wfm_header), 5);
    wfm_baseline_g = zeros(length(wfm_header), 5);
    for i = 1:length(wfm_header)
        wfm_baseline_t(i,1) = 0;
        wfm_baseline_t(i,end) = t(end);
        wfm_baseline_t(i,2:4) = [wfm_header(i).start, (wfm_header(i).start+wfm_header(i).end)/2, wfm_header(i).end] - adc_start_time;
        wfm_baseline_t(i,2:4) = wfm_baseline_t(i,2:4) / 1000; % time in [us]
        wfm_baseline_g(i,3) = wfm_header(i).maxG;
        wfm(i,:) = interp1( wfm_baseline_t(i,:), wfm_baseline_g(i,:), t, 'linear');
    end
    wfm = wfm/100;
end





function [wfm_header] = get_wfm_header_from_rth_header(rth_header)
% rth_header should contain specific keys such as "user_grad10_duration",
% "user_grad10_start", "user_grad10_end", "user_grad10_maxG"

header_fields = fieldnames(rth_header);
wfm_header = struct;

for i = 1:length(header_fields)
    each_header = header_fields{i};
    if contains(each_header, 'user_grad')
        each_header_split = split(each_header,"_");
        amplitude_index = each_header_split{2};
        amplitude_index = str2num(amplitude_index(5:end));
        wfm_header(amplitude_index).(each_header_split{3}) = rth_header.(each_header);
    end
end

end

