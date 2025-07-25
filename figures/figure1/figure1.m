load('measured.mat')
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
    else
        disp('Flip angle not found.');
    end

    
    if flip_angle == 20 || flip_angle == 80
    else
        continue
    end
    

    plotctr = plotctr + 1;
    % find the index that matches.
    faNumbers = cellfun(@(x) str2double(regexp(x, 'FA(\d+)', 'tokens', 'once')), final_label);
    %faNumbers = cellfun(@(x) str2double(regexp(x, 'FA = (\d+)', 'tokens', 'once')), final_label);
    index = find(faNumbers == flip_angle);

     n_f = size(contrast_per_fa, 2);
     n_per_t = length(timescales) / n_f;

     timescale_measured = timescales( (((index-1)*n_per_t)+1):(index*n_per_t));
     measured_contrast = contrast_per_fa(:,index);

    TR = 4.95;
    start_sig = bssfp(flip_angle, TR, 0, 655, 43, 500/TR);

    stdev = real(start_sig / 15);

    tr_idx_sim = floor(timescale_measured(1)/TR);
    %tr_idx_sim = floor(timescale_measured{1}(1)/TR);
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

%% plot the CNR
% set up export variables
fas = [];
timescales = {};

fitted_slopes = {};

figure;
set(gcf, 'Position', [0,0,1000, 800]) %,'Color', 'k')
hold on
labels = final_label;
final_label = {};

labelcount = 1;
for passid = 1:2

for i =1:length(file_names)
    pass1 = passid-1;
    pass1
    if (i == 2 || i == 4  || i == 9 || i == 11) && (pass1 == 0) 
    elseif (i == 7 || i == 8) && (pass1 == 1)
    else
        continue
    end
    

    fprintf("------ plotting CNR %i of %i --------\n", i, length(file_names));
    timescale = times{i};

    % force timescale to end at 1.5
    timescale = timescale(1:find(timescale > 1500,1));
     
    flip_angle = regexp(labels{i}, '\d+', 'match');
    flip_angle = str2num(flip_angle{2});

    cnr_line = contrast_per_fa(:,i);

    timescales{labelcount} = timescale;
    final_label{labelcount} = sprintf('FA = %i %c', flip_angle, char(176));

     labelcount = labelcount + 1;
     n_start = 10;
     % let's get the tag persistence.
     persistence{labelcount} = timescale(find(cnr_line(n_start:end) < 6,1));
     h = plot(timescale(n_start:end), cnr_line(n_start:length(timescale)), 'LineWidth', 4);

     % Get the color of the existing line
     line_color = get(h, 'Color');

     xlim([0, 1500]);
     ylim([0,50])

     xlabel('ms'); %, 'Color', 'w')
     ylabel('contrast (a.u.)'); %, 'Color', 'w')
end
end

for labelidx = 1:length(final_label)
    % HACK. we know the contrast so let's append it'
    fatext = final_label{labelidx};
    fa = fatext(end-3:end-2)
    if str2num(fa) == 20 || str2num(fa) == 40 || str2num(fa) == 60 || str2num(fa) == 80
        fatext = ['bSSFP ', (fa), fatext(end)];
    else
        fatext = ['GRE    ', (fa), fatext(end)];
    end
    final_label{labelidx} = fatext;
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


