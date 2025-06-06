clearvars;
clc;

% load file and skip first 12 lines
T = readtable('LV-Mask12May2025_14h10m00s_export.csv', 'NumHeaderLines', 12);

% number of frames to process
nFrames = 30;

% create cell array to store the time of selected frames
frames_time = cell(nFrames, 1);

% loop over the 30 rows of the table
for i = 1:nFrames
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
end

all_times = cell2mat(frames_time);
% determine automatic start time and end time from cell array
t_start_auto = min(all_times);
t_end_auto = max(all_times);

% load video
video = VideoReader('20181130T121536_Bmode_coherent_FIR_Apical 3C_ave-10.mp4');
% create time axis where each element corresponds to the time of each frame
time_axis = 0:1/video.FrameRate:video.NumFrames/video.FrameRate;

% find frames closest to automatic start time and end time
[~, frame_start_idx] = min(abs(time_axis - t_start_auto));
[~, frame_end_idx] = min(abs(time_axis - t_end_auto));

% load mat-file
load('20181130T121536_piv.mat');

% automatically determined start and end frame
start_frame = frame_start_idx;
end_frame = frame_end_idx;

% PIV variabelen
vectors_all = permute(vfi.vectors(:,:,:,[2,1]), [2 1 3 4]); % in vfi, vectors [z x t component] -component [vz vx]
grid = permute(vfi.grid(:,:,[2,1]), [2 1 3]);               % in vfi.grid [z x component] -component [z x]
dt = vfi.dt;
frame_rate = vfi.frame_rate; 
mPoly = vfi.mPoly;

% dimensions
[ny, nx, nt_total, ~] = size(vectors_all);
n_points = ny * nx;

% select frames in the automatic range
vectors = vectors_all(:,:,start_frame:end_frame,:);
nt = end_frame - start_frame + 1;

% time axis
tax = (start_frame:end_frame) * dt;

% positions
pos = single(reshape(grid, [], 2));

% vectors
vecs = single(zeros(nt, n_points, 2));
for t = 1:nt
    vframe = squeeze(vectors(:,:,t,:)); % (ny x nx x 2)
    vecs(t,:,1) = reshape(vframe(:,:,1), [], 1); % vx
    vecs(t,:,2) = reshape(vframe(:,:,2), [], 1); % vy
end

% wall
if size(mPoly, 1) > 1
    norm_wall = compute_polygon_normals(mPoly(:,1), mPoly(:,2));
    wall = single([mPoly norm_wall(1:end-1,:)]);  % (N-1 x 4)
else
    wall = single(zeros(0,4));
end

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
struct.tax = tax;
struct.pos = pos;
struct.vecs = vecs;
struct.wall = wall;
struct.in_shape = in_shape;
struct.out_shape = out_shape;
struct.extents = extents;
struct.bounding_box = bounding_box;

% save mat-file
save('patient18.mat', 'struct');