clear; close all; clc;
load('measured2.mat')

%% plot the CNR
% set up export variables
fas = [];
timescales = {};

fitted_slopes = {};

figure;
set(gcf, 'Position', [0,0,1000, 800]) %,'Color', 'k')
hold on
final_label = {};

labelcount = 1;
for i =1:length(file_names)
    fprintf("------ plotting CNR %i of %i --------\n", i, length(file_names));
    timescale = times{i};

    if isempty(timescale)
        continue
    end


    % force timescale to end at 1.5
    timescale = timescale(1:find(timescale > 1500,1));
     
    flip_angle = regexp(labels{i}, '\d+', 'match');
    flip_angle = str2num(flip_angle{2});

    cnr_line = contrast_per_fa{i};

    timescales{labelcount} = timescale;
    final_label{labelcount} = sprintf('FA = %i %c', flip_angle, char(176));

     labelcount = labelcount + 1;
     n_start = 17;
     % let's get the tag persistence.
     persistence{labelcount} = timescale(find(cnr_line(n_start:end) < 6,1));
     plot(timescale(n_start:end), cnr_line(n_start:length(timescale)), 'LineWidth', 4);

     %plot(timescale(n_start:10:end), cnr_line(n_start:10:length(timescale)), 'x', 'MarkerSize', 10, 'LineWidth', 1.5);
     xlim([0, 1500]);
     ylim([0,50])

     xlabel('ms'); %, 'Color', 'w')
     ylabel('contrast (a.u.)'); %, 'Color', 'w')
end

l = legend(final_label);
%l.Color = 'k';l.TextColor = 'w';l.EdgeColor='none';
%legend(final_label, 'Location', 'best outside')
grid on;
grid minor;
gray = [0.7, 0.7, 0.7];
set(gca,'FontName','Times','FontSize',36,...
    'XColor','k','YColor', 'k','GridColor',gray,...
    'LineWidth',1.5);
title({'Tagging contrast over time',''},'FontSize',48)

export_fig('gre_bssfp_measured.pdf')

% Print each element pair on a new line
for ival = 2:length(persistence)
    ival_ = ival -1;
    str1 = num2str(persistence{ival});
    str2 = final_label{ival-1};
    % Print the pair
    fprintf('%s\t\t%s\n', str1, str2);
end



%%
load('simulation.mat')

num_plots = size(seq_files, 1);

plothandles = gobjects(1, num_plots*2);
legends = cell(1, num_plots*2);

plotctr = 0;

figure;
set(gcf, 'Position', [0,0,1000, 800]) %,'Color', 'k')
hold on
for i = 1:num_plots

    str_seq = seq_files(i,:);

    tokens = regexp(str_seq, 'FA(\d+)', 'tokens');

    if contains(str_seq, 'trueFISP')
    else
        continue
    end

    % Extract the flip angle
    if ~isempty(tokens)
        flip_angle = str2double(tokens{1}{1});
        %disp(['Flip angle: ', num2str(flip_angle)]);
    else
        disp('Flip angle not found.');
    end

    if flip_angle == 20 || flip_angle == 80

    else
        continue
    end

    plotctr = plotctr + 1;
    % find the index that matches.
    faNumbers = cellfun(@(c) sscanf(c, 'FA = %d'), final_label);
    index = find(faNumbers == flip_angle);

     timescale_measured = timescales{index};
     measured_contrast = contrast_per_fa{index};
    

    TR = 5.19;
    start_sig = bssfp(flip_angle, TR, 0, 683, 77, 500/TR);
    stdev = real(start_sig / 10); % this is a simulated stdev to match measured data.


    tr_idx_sim = floor(timescale_measured(1)/TR);
    tr_scale = contrast_per_angle(tr_idx_sim, i) / measured_contrast(1);

    hline = plot((1:1000)*TR, contrast_per_angle(:,i)/stdev, 'LineWidth', 4);
    lineColor = hline.Color;
    hmeas = plot(timescale_measured(12:6:end)-(TR*10.5), measured_contrast(12:6:length(timescale_measured)),  'x', 'MarkerSize', 15, 'Color', lineColor, 'LineWidth', 5 );

    xlim([0,1500]);
    ylim([0,50]);
    hold on;

    legends{2*plotctr} = ['Measured FA = ', num2str(flip_angle), char(176)];
    legends{2*plotctr-1} =  ['Simulation FA = ', num2str(flip_angle), char(176)];
    
    plothandles(2*plotctr) = hmeas;
    plothandles(2*plotctr-1) = hline;

    disp( ['FA: ', num2str(flip_angle), 'Persistence: ', num2str(find(contrast_per_angle(:,i)/stdev < 6, 1)*TR)]);
end


legend(plothandles(1:(2*plotctr)), legends{1:(2*plotctr)}, 'Location', 'best');
ylabel('CNR')
xlabel('ms')
grid on;
grid minor;
gray = [0.7, 0.7, 0.7];
set(gca,'FontName','Times','FontSize',36,...
    'XColor','k','YColor', 'k','GridColor',gray,...
    'LineWidth',1.5);
title({'bSSFP Tagging Contrast ',''},'FontSize',48)

export_fig('plot_simulation.pdf')