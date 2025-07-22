function plot_grid(data_folder, pattern, frames, varargin)
    % Plot grid of trueFISP images across flip angles and frames
    %
    % Inputs:
    %   data_folder - path to folder containing .mat files
    %   pattern - string pattern to search
    %   frames  - list of frames to display
    %   
    % Optional Name-Value pairs:
    %   'show_gridlines' - logical, show thin white gridlines (default: false)
    %   'enhance_contrast' - logical, apply contrast enhancement (default: false)
    %
    % Example usage:
    %   plot_trueFISP_grid('/path/to/data', 1, 4, 20);
    %   plot_trueFISP_grid('/path/to/data', 1, 4, 20, 'show_gridlines', true, 'enhance_contrast', true);
    
    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'show_gridlines', false, @islogical);
    addParameter(p, 'enhance_contrast', false, @islogical);
    parse(p, varargin{:});
    
    show_gridlines = p.Results.show_gridlines;
    enhance_contrast = p.Results.enhance_contrast;

    % Find all trueFISP files for the specified repetition
    files = dir(fullfile(data_folder, pattern));
    
    if isempty(files)
        error('No trueFISP files found for repetition %d in %s', repetition, data_folder);
    end
    
    % Extract flip angles from filenames and load data
    images = {};
    actual_flip_angles = [];
    
    for i = 1:length(files)
        filename = files(i).name;
        
        % Extract flip angle using regex pattern *_FAXX_*
        fa_match = regexp(filename, '_FA(\d+)_', 'tokens');
        if ~isempty(fa_match)
            fa = str2double(fa_match{1}{1});
            
            % Load the file
            filepath = fullfile(data_folder, filename);
            data = load(filepath);
            
            if isfield(data, 'image_stcr')
                images{end+1} = rectify_pk(data.image_stcr);
                actual_flip_angles(end+1) = fa;
                fprintf('Loaded: %s (FA%d)\n', filename, fa);
            else
                warning('Variable "ir" not found in %s', filename);
            end
        else
            warning('Could not extract flip angle from filename: %s', filename);
        end
    end

    %actual_flip_angles = [20,40,60,80];
    
    % Sort by flip angle
    [actual_flip_angles, sort_idx] = sort(actual_flip_angles);
    images = images(sort_idx);
    
    if isempty(images)
        error('No trueFISP images loaded. Check data folder and repetition number.');
    end
    
    % Determine frame indices to display
    max_frames = min(cellfun(@(x) size(x, 3), images));
    frame_indices = frames; %start_frame:frame_gap:(start_frame + (num_frames-1)*frame_gap);
    
    % Ensure we don't exceed available frames
    frame_indices = frame_indices(frame_indices <= max_frames);
    actual_num_frames = length(frame_indices);
    
    
    % Create figure with black background
    figure('Position', [100, 100, 200*actual_num_frames, 200*length(actual_flip_angles)], ...
           'Color', 'black');
    
    % Calculate crop size (slightly more than half of typical image size)
    sample_img = images{1}(:, :, 1);
    [h, w] = size(sample_img);
    crop_size = round([h*0.4, w*0.4]); % 40% of original size instead of 50%
    
    % Store subplot positions for labels
    subplot_positions = {};
    
    % Plot with tight spacing and proper margins for title
    set(gcf, 'Color', 'black');
    
    for fa_idx = 1:length(actual_flip_angles)
        for frame_idx = 1:actual_num_frames
            % Use tight subplot with more top margin for title
            h_subplot = subtightplot(length(actual_flip_angles), actual_num_frames, ...
                        (fa_idx-1)*actual_num_frames + frame_idx, [0.01 0.01], [0.12 0.08], [0.08 0.02]);
            
            % Store position for first subplot of each row
            if frame_idx == 1
                subplot_positions{fa_idx} = get(h_subplot, 'Position');
            end
            
            % Extract and process image
            img = images{fa_idx}(:, :, frame_indices(frame_idx));
            img = rectify_pk(img);
            img = crop_half_FOV(img, crop_size);
            img=flipud(real(img));
            
            % Apply contrast enhancement if requested (with safety checks)
            if enhance_contrast
                img = imadjust(img ./ max(vec(img)), [0 0.3]);
            end
            
            imagesc(img);
            colormap gray;
            axis image off;
            
            % Add frame labels at top (white text)
            if fa_idx == 1
                title(sprintf('Frame %d', frame_indices(frame_idx)), ...
                      'FontSize', 12, 'Color', 'white', 'FontWeight', 'bold');
            end
            
            % Add gridlines if requested
            if show_gridlines
                hold on;
                xlims = xlim;
                ylims = ylim;
                % Right border
                if frame_idx < actual_num_frames
                    plot([xlims(2) xlims(2)], ylims, 'w-', 'LineWidth', 0.5);
                end
                % Bottom border
                if fa_idx < length(actual_flip_angles)
                    plot(xlims, [ylims(2) ylims(2)], 'w-', 'LineWidth', 0.5);
                end
                hold off;
            end
        end
    end
    
    % Add flip angle labels using stored positions
    for fa_idx = 1:length(actual_flip_angles)
        subplot_pos = subplot_positions{fa_idx};
        
        % Calculate text position (left of the first subplot in each row)
        text_x = subplot_pos(1) - 0.04;  % Left of subplot
        text_y = subplot_pos(2) + subplot_pos(4)/2;  % Middle of subplot height
        
        annotation('textbox', [text_x, text_y, 0.03, 0.05], ...
                   'String', sprintf('FA%d°', actual_flip_angles(fa_idx)), ...
                   'FontSize', 14, 'Color', 'white', 'FontWeight', 'bold', ...
                   'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                   'EdgeColor', 'none', 'BackgroundColor', 'none');
    end
    
    % Add overall title with proper spacing
    annotation('textbox', [0, 0.95, 1, 0.05], ...
               'String', sprintf('trueFISP Images'), ...
               'FontSize', 16, 'Color', 'white', 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
               'EdgeColor', 'none', 'BackgroundColor', 'none');
    
    % Print summary
    fprintf('\nGrid created:\n');
    fprintf('- Flip angles: %s\n', mat2str(actual_flip_angles));
    fprintf('- Frames displayed: %s\n', mat2str(frame_indices));
    fprintf('- Grid size: %d × %d\n', length(actual_flip_angles), actual_num_frames);
    fprintf('- Crop size: [%d, %d]\n', crop_size(1), crop_size(2));
    fprintf('- Gridlines: %s\n', mat2str(show_gridlines));
    fprintf('- Contrast enhancement: %s\n', mat2str(enhance_contrast));
end

% Helper function for tight subplot spacing
function h = subtightplot(m, n, p, gap, marg_h, marg_w)
    % Tight subplot with minimal spacing
    % gap: [gap_vert, gap_horiz] - spacing between subplots
    % marg_h: [lower, upper] - vertical margins
    % marg_w: [left, right] - horizontal margins
    
    if nargin < 4, gap = [0.01 0.01]; end
    if nargin < 5, marg_h = [0.05 0.05]; end
    if nargin < 6, marg_w = [0.05 0.05]; end
    
    if numel(gap) == 1
        gap = [gap gap];
    end
    if numel(marg_h) == 1
        marg_h = [marg_h marg_h];
    end
    if numel(marg_w) == 1
        marg_w = [marg_w marg_w];
    end
    
    axh = (1-sum(marg_h)-(m-1)*gap(1))/m;
    axw = (1-sum(marg_w)-(n-1)*gap(2))/n;
    
    py = 1-marg_h(2)-axh;
    
    ii = rem(p-1, n) + 1;
    jj = ceil(p/n);
    
    px = marg_w(1) + (ii-1)*(axw + gap(2));
    py = py - (jj-1)*(axh + gap(1));
    
    h = axes('position', [px py axw axh]);
end