clearvars;
clc;

% path to directory
directory = 'D:\San\LVSegmentation';

% load CSV-file and skip first 12 lines
T = readtable(fullfile(directory, 'LV-Mask12May2025_14h10m00s_export.csv'), 'NumHeaderLines', 12);

% extract token (still nested: cell array of 1x1 cell arrays)
tokens = regexp(T.Var4, '\[([^\]]+)\]', 'tokens', 'once');

% flatten to simple cell array of strings
flat = cellfun(@(c) c{1}, tokens, 'UniformOutput', false);
frames_time = str2double(flat);

% load video
video = VideoReader(fullfile(directory, '20181130T121536_Bmode_coherent_FIR_Apical 3C_ave-10.mp4'));
% create time axis where each element corresponds to the time of each frame
time_axis = 0:1/video.FrameRate:video.NumFrames/video.FrameRate;

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
[~, start_frame] = min(abs(time_axis_piv*tax_sf - t_start_auto));
[~, end_frame] = min(abs(time_axis_piv*tax_sf - t_end_auto));

% PIV variabelen
vectors_all = permute(vfi.vectors(:,:,:,[2,1]), [2 1 3 4]); % in vfi, vectors [z x t component] -component [vz vx]
grid = permute(vfi.grid(:,:,[2,1]), [2 1 3]);               % in vfi.grid [z x component] -component [z x]
dt = 1/vfi.frame_rate; % vfi.dt;
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
% load interpolated wall
load(fullfile(directory, 'new_wall.mat'), 'new_wall');

% extract time column
wall_time = new_wall(:,1,1);

% determine number of wall points per frame
num_wall_points = size(new_wall,2);

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
struct.wall = reshape(new_wall, [], 5);
struct.in_shape = in_shape;
struct.out_shape = out_shape;
struct.extents = extents;
struct.bounding_box = bounding_box;
struct.bmodes = bmodes;

% save MAT-file
save(fullfile(directory, 'patient18_data.mat'), 'struct');