%%
clearvars;
clc;

% path to directory
directory = 'D:\San\LVSegmentation';

% load CSV-file and skip first 9 lines
T = readtable(fullfile(directory, 'LV-Mask12May2025_14h10m00s_export.csv'), 'NumHeaderLines', 9, 'Delimiter', ',');

% Step 1: data cleaning
% identify rows where 'spatial_coordinates' is exactly the string '[]'
isEmptyCoord = strcmp(T.spatial_coordinates, '[]');
% remove those rows
T(isEmptyCoord, :) = [];

% Step 2: extract times
temporal_coordinates = cellfun(@(s) sscanf(s(2:end-1), '%f'), T.temporal_coordinates);

% Step 3: extract and clean coordinates from the remaining rows
cleaned_coordinates = cellfun(@(s) sscanf(s(2:end-1), '%f,').', T.spatial_coordinates, 'UniformOutput', false);
cleaned_coordinates = cellfun(@(v) v(2:end), cleaned_coordinates, 'UniformOutput', false);

% convert each interleaved coordinate vector to an Nx2 matrix
coordinate_matrices = cellfun(@(v) reshape(v, 2, []).', cleaned_coordinates, 'UniformOutput', false);

% Step 4: extract labels
phase = cellfun(@(s) extractBetween(s, ':"', '"}'), T.metadata, 'UniformOutput', false);
% convert to flat cell array of strings
phase = vertcat(phase{:});

%% Step 5: load video and time axis
video = VideoReader(fullfile(directory, '20181130T121536_Bmode_coherent_FIR_Apical 3C_ave-10.mp4'));
% create time axis where each element corresponds to the time of each frame
time_axis = 0:1/video.FrameRate:video.NumFrames/video.FrameRate;

% determine automatic start time and end time from cell array
t_start_auto = min(temporal_coordinates);
t_end_auto = max(temporal_coordinates);

% load PIV MAT-file
load(fullfile(directory, '20181130T121536_piv.mat'), 'vfi', 'bmodes');
time_axis_vfi = 0:1/vfi.frame_rate:size(vfi.vectors,3)/vfi.frame_rate;
movie2vfi_sf = time_axis(end)/time_axis_vfi(end);

% find frames closest to automatic start time and end time
[~, frame_start_idx] = min(abs(time_axis_vfi*movie2vfi_sf - t_start_auto));
[~, frame_end_idx] = min(abs(time_axis_vfi*movie2vfi_sf - t_end_auto));

%% Step 6: convert pixel coordinates to mm
% define the shift (the old starting point of the frame)
shift_x = 92;
shift_y = 32;

% extent of image in mm
extent = [120, 120];          % x, y
origin = [-extent(1)/2, 0];   % center top is origin

% define crop region: [x, y, width, height]
crop_rect = [shift_x, shift_y, 555-shift_x, 395-shift_y];
mm_per_pixel = extent./crop_rect(3:4);

% move the coordinates by subtracting (92, 32)
adjusted_pixel = cellfun(@(x) x - [shift_x, shift_y], coordinate_matrices, 'UniformOutput', false);
adjusted_mm = cellfun(@(x) x.*mm_per_pixel + origin, adjusted_pixel, 'UniformOutput', false);

%% Step 7: condition polygon based on cardiac phase
for i = 1:numel(adjusted_mm)
    % create closed poly
    polyline_closed = adjusted_mm{i};
    if ~all(polyline_closed(1, :) == polyline_closed(end,:))
        polyline_closed(end, :) = polyline_closed(1,:);   
    end
    % select polyline based on flow phase type
    switch phase{i}
        case "Inflow"
            % use open polyline for inflow
            in_index = 4:23;
        case "Outflow"
            % use open polyline for outflow
            in_index = 1:21;
        otherwise
            % use closed polyline as default
            in_index = 1:size(polyline_closed, 1);
    end
    
    % determine which points to set as NaN
    out_index = setdiff(1:size(polyline_closed,1), in_index);
    polyline = polyline_closed;
    if ~isempty(out_index)
        polyline(out_index, :) = nan(length(out_index), 2);
    end
    
    % store the polyline mm coordinates per frame in the cell array
    polyline_mm_all{i} = polyline;

    % Step 8: compute normals to polyline
    % compute normal vectors for polyline_mm and skip NaNs for open polylines
    valid_idx = all(~isnan(polyline), 2);
    % keep only the valid rows 
    valid_polyline = polyline(valid_idx, :);
    
    % compute normals for the full closed polygon
    normals = compute_polygon_normals(polyline_closed(:,1), polyline_closed(:,2));
    normals = normals(:, :);
    if ~isempty(out_index)
        normals(out_index, :) = nan(length(out_index), 2);
    end
    
    % store normals per frame in the cell array
    normals_all{i} = normals;

    % Step 9: create a matrix containing the wall polygon with normals per time point
    % combine coordinates, normals, and time into one matrix
    t = repmat(temporal_coordinates(i), size(polyline,1), 1);
    wall{i} = [polyline, normals, t];
end

% Step 10: check that it fits the image
figure;
% find the index of the frame closest to the specified time
[~, frame_idx] = min(abs(temporal_coordinates(1) - time_axis));
frame = read(video, frame_idx);
frame = imcrop(frame, crop_rect);

% generate physical axis values in mm
x_axis = (1:crop_rect(3)).*mm_per_pixel(1);
x_axis = x_axis - mean(x_axis);
y_axis = (0:crop_rect(4)-1).*mm_per_pixel(2);
im = imagesc(x_axis, y_axis, frame); 
axis image;
hold on;
wall_ = wall{1};
% poly = scatter(adjusted_mm(:,1), adjusted_mm(:,2), 30, 'b', 'filled');
poly = plot(wall_(:,1), wall_(:,2), 'b-o', 'LineWidth', 2);
norms = quiver(wall_(:,1), wall_(:,2), wall_(:, 3), wall_(:,4), 0.3, 'Color','y');
hold off;
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
set(gca, 'YDir', 'reverse');
title('Coordinates in mm');
for i = 1:numel(temporal_coordinates)
    [~, frame_idx] = min(abs(temporal_coordinates(i) - time_axis));
    frame = read(video, frame_idx);
    frame = imcrop(frame, crop_rect);
    im.CData = frame;
    wall_ = wall{i};
    % set(poly, 'XData', adjusted_mm(:,1), 'YData', adjusted_mm(:,2))
    set(poly, 'XData', wall_(:,1), 'YData', wall_(:,2));
    set(norms, 'XData', wall_(:,1), 'YData', wall_(:,2), 'UData', wall_(:, 3), 'VData', wall_(:,4));
    pause(0.1);
end

% save figure frames as video
v = VideoWriter(fullfile(directory, 'wall_with_normals.mp4'), 'MPEG-4');
v.FrameRate = 10;
open(v);

for i = 1:numel(temporal_coordinates)
    [~, frame_idx] = min(abs(temporal_coordinates(i) - time_axis));
    frame = read(video, frame_idx);
    frame = imcrop(frame, crop_rect);
    im.CData = frame;
    wall_ = wall{i};
    set(poly, 'XData', wall_(:,1), 'YData', wall_(:,2));
    set(norms, 'XData', wall_(:,1), 'YData', wall_(:,2), ...
        'UData', wall_(:,3), 'VData', wall_(:,4));
    
    drawnow;
    frame_out = getframe(gcf); 
    writeVideo(v, frame_out);
end

close(v);

%% Step 11: trim PIV data and interpolate wall onto each HFR frame
vfi_frame_fudge_factor = 14;
vectors_all = permute(vfi.vectors(:,:,frame_start_idx-vfi_frame_fudge_factor:frame_end_idx-vfi_frame_fudge_factor,[2,1]), [2 1 3 4]);  % in vfi, vectors [z x t component] -component [vz vx]
grid = permute(vfi.grid(:,:,[2,1]), [2 1 3]);   % in vfi.grid [z x component] -component [z x]
dt = 1/vfi.frame_rate;   % vfi.dt;
vfi_dt = vfi.dt;
frame_rate = vfi.frame_rate; 

% trim HFR B-mode frames
bmode = bmodes.imagedata(:,:, frame_start_idx:frame_end_idx);

% dimensions
[nx, ny, nt, ~] = size(vectors_all);
n_points = ny * nx;

% positions
pos = single(reshape(grid, [], 2));

% vectors
vecs = single(zeros(size(vectors_all,3), n_points, 2));
for t = 1:size(vectors_all,3)
    vframe = squeeze(vectors_all(:,:,t,:));  % (nx x ny x 2)
    vecs(t,:,1) = reshape(vframe(:,:,1), [], 1);   % vx
    vecs(t,:,2) = reshape(vframe(:,:,2), [], 1);   % vy
end

% wall
% convert cell array of wall data into a matrix
wallmat = cell2mat(wall');
% number of points per frame in the wall data
num_points = size(wall{1}, 1);
% total number of segmented frames
num_segmented_frames = length(wall);
% normalize B-mode time axis to range from 0 to 1
t_ax_norm = (temporal_coordinates-temporal_coordinates(1))./(temporal_coordinates(end)-temporal_coordinates(1));  
t_ax_vfi_norm = linspace(0,1,size(vectors_all,3));
% calculate start and end times for HFR video segment
t_start_HFR = frame_start_idx * dt;
t_end_HFR = frame_end_idx * dt;

% Step 12: interpolate wall
% extract all x coordinates from wall and interpolate over hfr time axis
x_all = cell2mat(cellfun(@(x) x(:, 1), wall, 'UniformOutput', false))';
x_hfr = 1e-3.*interp1(t_ax_norm, x_all, t_ax_vfi_norm, "linear");

% extract all y coordinates from wall and interpolate over hfr time axis
y_all = cell2mat(cellfun(@(x) x(:, 2), wall, 'UniformOutput', false))';
y_hfr = 1e-3.*interp1(t_ax_norm', y_all, t_ax_vfi_norm, "linear");

% extract all x components of normal vectors and interpolate over hfr time axis
nx_all = cell2mat(cellfun(@(x) x(:, 3), wall, 'UniformOutput', false))';
nx_hfr = interp1(t_ax_norm', nx_all, t_ax_vfi_norm, "linear");

% extract all y components of normal vectors and interpolate over hfr time axis
ny_all = cell2mat(cellfun(@(x) x(:, 4), wall, 'UniformOutput', false))';
ny_hfr = interp1(t_ax_norm', ny_all, t_ax_vfi_norm, "linear");

% repeat each HFR timestamp for all points in the frame
tax = linspace(t_start_HFR, t_end_HFR, numel(t_ax_vfi_norm))';
t_hfr = repmat(tax, 1, num_points);

% interpolate phase
phase_ind = single(strcmp(phase, 'Inflow'));
phase_ind_hfr = interp1(t_ax_norm', phase_ind, t_ax_vfi_norm, "nearest");
phase_hfr = repmat("Inflow", size(phase_ind_hfr));
phase_hfr(phase_ind_hfr==0) = "Outflow";

% final wall matrix
wall_hfr = cat(3, t_hfr, x_hfr, y_hfr, nx_hfr, ny_hfr);
num_wall_points = size(wall_hfr,2);

% Step 13: plot to check they are the same
figure();
im = pcolor(bmodes.x.*1e3, bmodes.z.*1e3, bmodes.imagedata(:,:,1));
shading flat;
colormap gray;
axis image;
hold on;
wall_ = wall{1};
poly_hfr = plot(wall_hfr(1,:,2).*1e3, wall_hfr(1,:,3).*1e3, 'b', 'LineWidth', 2);
norms_hfr = quiver(wall_hfr(1,:,2).*1e3, wall_hfr(1,:,3).*1e3, ...
                   wall_hfr(1,:,4), wall_hfr(1,:,5), ...
                   0.3, 'Color','y');
quiv_hfr = quiver(grid(:,:,1).*1e3, grid(:,:,2).*1e3, vectors_all(:,:,1,1), vectors_all(:,:,1,2), 3);
poly = plot(wall_(:,1), wall_(:,2), 'r--', 'LineWidth', 1);

hold off;
xlabel('X (mm)');
ylabel('Y (mm)');
set(gca, 'YDir', 'reverse');
title(['Coordinates in mm']);
for i = 1:size(wall_hfr,1)
    [~, frame_idx] = min(abs(t_ax_vfi_norm(i) - t_ax_norm));
    im.CData = bmodes.imagedata(:,:,i);
    wall_ = wall{frame_idx};
    % set(poly, 'XData', adjusted_mm(:,1), 'YData', adjusted_mm(:,2))
    set(poly, 'XData', wall_(:,1), 'YData', wall_(:,2));
    set(poly_hfr, 'XData', wall_hfr(i,:,2).*1e3, 'YData', wall_hfr(i,:,3).*1e3);
    set(norms_hfr, 'XData', wall_hfr(i,:,2).*1e3, 'YData', wall_hfr(i,:,3).*1e3, ...
                   'UData', wall_hfr(i,:, 4), 'VData', wall_hfr(i,:,5));
    set(quiv_hfr, 'UData', vectors_all(:,:, i,1), 'VData', vectors_all(:,:, i,2));

    % save figure frames as video
    drawnow;
    frame_out = getframe(gcf);
    writeVideo(v, frame_out);
end
close(v);

%% Final step: save data in MAT-file

% shapes
in_shape = [nt, ny, nx];
out_shape = [nt, 64, 64];

% extents
x_vals = grid(:,:,1);
y_vals = grid(:,:,2);
extents = [max(tax), range(x_vals(:)), range(y_vals(:))];

% bounding box
bounding_box = [min(x_vals(:)), max(x_vals(:)), min(y_vals(:)), max(y_vals(:))];

% build structure
in_struct.tax = tax;
in_struct.pos = pos;
in_struct.vecs = vecs;
in_struct.wall = wall_hfr;
in_struct.in_shape = in_shape;
in_struct.out_shape = out_shape;
in_struct.extents = extents;
in_struct.bounding_box = bounding_box;
in_struct.bmodes = bmodes;
in_struct.wall_shape = size(wall_hfr);
in_struct.phase = phase_hfr;

% save MAT-file
save(fullfile(directory, 'patient18_variables.mat'), 'in_struct');