clearvars;
clc;

% load CSV-file and skip first 12 lines
T = readtable('LV-Mask12May2025_14h10m00s_export.csv', 'NumHeaderLines', 12);

% load B-mode video
video = VideoReader('20181130T121536_Bmode_coherent_FIR_Apical 3C_ave-10.mp4');
% create time axis where each element corresponds to the time of each frame
time_axis = 0:1/video.FrameRate:video.NumFrames/video.FrameRate;

% number of frames to process
nFrames = 30;

% create cell arrays to store the variables
all_adjusted_pixel = cell(nFrames, 1);   % pixel coordinates
all_adjusted_mm = cell(nFrames, 1);      % mm coordinates
phase = cell(nFrames, 1);                % phase (inflow or outflow)
polyline_mm_all = cell(nFrames, 1);      % polyline mm coordinates
frames_time = cell(nFrames, 1);          % time of selected frames
normals_all = cell(nFrames, 1);          % normals per frame

% loop over the 30 rows of the table
for i = 1:nFrames
    
    % take string from column 5
    raw_string = T.Var5{i};
    % remove the first "[" and the last "]"
    raw_string = raw_string(2:end-1); 
    % split by commas
    number_strings = strsplit(raw_string, ','); 
    % convert everything to numbers
    numbers = str2double(number_strings);  
    % remove the first number (6)
    numbers = numbers(2:end);
    
    % divide into x and y
    x = numbers(1:2:end);
    y = numbers(2:2:end);    
    % combine into coordinates
    coordinate_pixel = [x(:) y(:)];

    % define the shift (the old starting point of the frame)
    shift_x = 92;
    shift_y = 32;
    % move the coordinates by subtracting (92, 32)
    adjusted_pixel = coordinate_pixel - [shift_x, shift_y];

    % store the moved pixel coordinates per frame in the cell array
    all_adjusted_pixel{i} = adjusted_pixel;

    % take string from column 4
    raw_string = T.Var4{i};
    % remove the first "[" and the last "]"
    raw_string = raw_string(2:end-1);
    % convert everything to numbers
    frame_time = str2double(raw_string);
    % find the index of the frame closest to the specified time
    [~, frame_idx] = min(abs(frame_time - time_axis));
    % store the times for each selected frame in the cell array
    frames_time{i} = time_axis(frame_idx);
    % read the specific video frame based on the frame index  
    frame = read(video, frame_idx);

    % define crop region: [x, y, width, height]
    crop_rect = [92, 32, 555-92, 395-32];
    % crop the frame
    cropped_frame = imcrop(frame, crop_rect);

    % resolution after cropping
    resX = crop_rect(3);
    resY = crop_rect(4);
    
    % real dimensions (in mm)
    realX = 120;   % from -60 to 60 mm
    realY = 120;   % from 0 to 120 mm
    
    % mm per pixel
    mm_per_pixel_x = realX / resX;
    mm_per_pixel_y = realY / resY;

    % generate physical axis values in mm
    x_axis = (1:resX).*mm_per_pixel_x;
    x_axis = x_axis - mean(x_axis);
    y_axis = (0:resY-1).*mm_per_pixel_y;

    % convert pixels to mm
    adjusted_mm = zeros(size(adjusted_pixel));
    adjusted_mm(:,1) = adjusted_pixel(:,1) * mm_per_pixel_x - realX/2;  % x from -60 to 60 mm
    adjusted_mm(:,2) = adjusted_pixel(:,2) * mm_per_pixel_y;            % y from 0 to 120 mm

    % store the moved mm coordinates per frame in the cell array
    all_adjusted_mm{i} = adjusted_mm;

    % get phase of cardiac cycle
    raw_string = T.Var6{i};
    % remove the first "[" and the last "]"
    raw_string = raw_string(2:end-1);
    % split by colon
    strings = strsplit(raw_string, ':');
    % phase is after colon
    strings = strings{2};
    % convert to string
    phase{i} = string(strings(2:end-1));
    
    % create closed poly
    polyline_closed = adjusted_mm;
    if ~all(adjusted_mm(1, :) == adjusted_mm(end,:))
        polyline_closed(end, :) = adjusted_mm(1,:);   
    end

    % select polyline based on flow phase type
    if strcmp(phase{i}, "Inflow")
    % use open polyline for inflow
    in_index = 4:23;
    elseif strcmp(phase{i}, "Outflow")
    % use open polyline for outflow
    in_index = 1:21;
    else
    % use closed polyline as default
    in_index = 1:size(polyline_closed, 1);
    end

    % determine which points to set as NaN
    out_index = setdiff(1:size(polyline_closed,1), in_index);

    % create final polyline with NaNs at out_index
    polyline_mm = polyline_closed(:, :);
    if ~isempty(out_index)
        polyline_mm(out_index, :) = nan(length(out_index), 2);
    end

    % store the polyline mm coordinates per frame in the cell array
    polyline_mm_all{i} = polyline_mm;

    % compute normal vectors for polyline_mm and skip NaNs for open polylines
    valid_idx = all(~isnan(polyline_mm), 2);
    % keep only the valid rows 
    valid_polyline = polyline_mm(valid_idx, :);
    
    % compute normals for the full closed polygon
    normals = compute_polygon_normals(polyline_closed(:,1), polyline_closed(:,2));
    normals = normals(:, :);
    if ~isempty(out_index)
        normals(out_index, :) = nan(length(out_index), 2);
    end
    
    % store normals per frame in the cell array
    normals_all{i} = normals;
    
    % combine coordinates, normals, and time into one matrix
    t = repmat(frames_time{i}, size(polyline_mm,1), 1);
    wall{i} = [polyline_mm, normals, t];

    % show the points in mm coordinates in a separate figure
    if i == 1
    figure;
    im = imagesc(x_axis, y_axis, cropped_frame); 
    axis image
    hold on;
    % poly = scatter(adjusted_mm(:,1), adjusted_mm(:,2), 30, 'b', 'filled');
    poly = plot(polyline_mm(:,1), polyline_mm(:,2), 'b-o', 'LineWidth', 2);
    norms = quiver(polyline_mm(:,1), polyline_mm(:,2), normals(:, 1), normals(:,2), 0.3, 'Color','y');
    hold off
    xlabel('X (mm)');
    ylabel('Y (mm)');
    axis equal;
    set(gca, 'YDir', 'reverse');
    title(['Coordinates in mm']);
    else
        im.CData = cropped_frame;
        % set(poly, 'XData', adjusted_mm(:,1), 'YData', adjusted_mm(:,2))
        set(poly, 'XData', polyline_mm(:,1), 'YData', polyline_mm(:,2))
        set(norms, 'XData', polyline_mm(:,1), 'YData', polyline_mm(:,2), 'UData', normals(:, 1), 'VData', normals(:,2))
        % pause
    end
end

%%
% frame A: closed polyline
frameA = 15;
coordsA = all_adjusted_pixel{frameA};
coordsA_closed = [coordsA; coordsA(1,:)];  % polyline from start to end

% plot figure
figure;
plot(coordsA_closed(:,1), coordsA_closed(:,2), 'g-o', 'LineWidth', 2)
axis equal
set(gca, 'YDir', 'reverse')
title('Closed polyline')

% frame B: open polyline at points 2 and 3 (inflow)
frameB = 20;
coordsB = all_adjusted_pixel{frameB};
coordsB_open = coordsB(4:23, :);  % segment from point 4 to 23

% plot figure
figure;
plot(coordsB_open(:,1), coordsB_open(:,2), 'r-o', 'LineWidth', 2)
axis equal
set(gca, 'YDir', 'reverse')
title('Open polyline at points 2 and 3 (inflow)')

% frame C: open polyline at point 22 (outflow)
frameC = 5;
coordsC = all_adjusted_pixel{frameC};
coordsC_outflow = [coordsC(1:21, :); NaN NaN; coordsC(23:end, :)];  % exclude point 22

% plot figure
figure;
plot(coordsC_outflow(:,1), coordsC_outflow(:,2), 'b-o', 'LineWidth', 2);
axis equal;
set(gca, 'YDir', 'reverse');
title('Open polyline at point 22 (outflow)');

%%
% convert cell array of wall data into a matrix
wallmat = cell2mat(wall');
% number of points per frame in the wall data
num_points = size(wall{1}, 1);
% total number of segmented frames
num_segmented_frames = length(wall);

% get frame time values as a vector
t_ax = cell2mat(frames_time);
% normalize B-mode time axis to range from 0 to 1
t_ax_norm = (t_ax-t_ax(1))./(t_ax(end)-t_ax(1));  

% path to directory
directory = 'D:\San\LVSegmentation';

% extract token (still nested: cell array of 1x1 cell arrays)
tokens = regexp(T.Var4, '\[([^\]]+)\]', 'tokens', 'once');

% flatten to simple cell array of strings
flat = cellfun(@(c) c{1}, tokens, 'UniformOutput', false);
frames_time = str2double(flat);

% determine automatic start time and end time from cell array
t_start_auto = min(frames_time);
t_end_auto = max(frames_time);

% load PIV MAT-file
load(fullfile(directory, '20181130T121536_piv.mat'), 'vfi', 'bmodes');

% PIV time axis
time_axis_piv = 0:1/vfi.frame_rate:size(vfi.vectors,3)/vfi.frame_rate;

% calculate time axis scale factor
tax_sf = time_axis(end)/time_axis_piv(end);

% find frames closest to automatic start time and end time
[~, frame_start] = min(abs(time_axis_piv*tax_sf - t_start_auto));
[~, frame_end] = min(abs(time_axis_piv*tax_sf - t_end_auto));

% number of frames in the selected HFR segment
num_frames = frame_end - frame_start;

% time step between frames in HFR video
dt_HFR = 0.0081405;
% create normalized time axis for HFR video segment
t_ax_hfr_norm = [0:num_frames-1]./num_frames;
% calculate start and end times for HFR video segment
t_start_HFR = frame_start * dt_HFR;
t_end_HFR = frame_end * dt_HFR;

% generate time vector over HFR segment duration
t_hfr = linspace(t_start_HFR, t_end_HFR, num_frames)';

% extract all x coordinates from wall and interpolate over hfr time axis
x_all = cell2mat(cellfun(@(x) x(:, 1), wall, 'UniformOutput', false))';
x_hfr = 1e-3.*interp1(t_ax_norm', x_all, t_ax_hfr_norm, "linear");

% extract all y coordinates from wall and interpolate over hfr time axis
y_all = cell2mat(cellfun(@(x) x(:, 2), wall, 'UniformOutput', false))';
y_hfr = 1e-3.*interp1(t_ax_norm', y_all, t_ax_hfr_norm, "linear");

% extract all x components of normal vectors and interpolate over hfr time axis
nx_all = cell2mat(cellfun(@(x) x(:, 3), wall, 'UniformOutput', false))';
nx_hfr = interp1(t_ax_norm', nx_all, t_ax_hfr_norm, "linear");

% extract all y components of normal vectors and interpolate over hfr time axis
ny_all = cell2mat(cellfun(@(x) x(:, 4), wall, 'UniformOutput', false))';
ny_hfr = interp1(t_ax_norm', ny_all, t_ax_hfr_norm, "linear");

% repeat each HFR timestamp for all points in the frame
t_hfr = repmat(t_hfr, 1, num_points);

% final wall matrix
new_wall = cat(3, t_hfr, x_hfr, y_hfr, nx_hfr, ny_hfr);

directory = 'D:\San\LVSegmentation';
save(fullfile(directory, 'new_wall.mat'), 'new_wall');

% show polyline and normal vectors over time
for i = 1:size(x_hfr,1)
    if i ==1
        figure()
        poly = plot(x_hfr(i,:), y_hfr(i,:), 'b-o', 'LineWidth', 2);
        hold on
        norms = quiver(x_hfr(i,:), y_hfr(i,:), nx_hfr(i,:), ny_hfr(i,:), 0.3);
        axis equal
    else
        set(poly, 'XData', x_hfr(i,:), 'YData', y_hfr(i,:))
        set(norms, 'XData', x_hfr(i,:), 'YData', y_hfr(i,:), 'UData', nx_hfr(i,:), 'VData', ny_hfr(i,:))

    end
    pause(0.1)
end