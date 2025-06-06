clearvars;
clc;

% load file and skip first 12 lines
T = readtable('LV-Mask12May2025_14h10m00s_export.csv', 'NumHeaderLines', 12);

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

    % load video
    video = VideoReader('20181130T121536_Bmode_coherent_FIR_Apical 3C_ave-10.mp4');
    % create time axis where each element corresponds to the time of each frame
    time_axis = 0:1/video.FrameRate:video.NumFrames/video.FrameRate;
   
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