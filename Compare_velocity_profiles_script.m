clearvars;
clc;

% path to directory
directory = 'D:\San\LVSegmentation';

% load raw and smoothed data
raw = load(fullfile(directory, 'patient18_variables.mat'));
smoothed = load(fullfile(directory, 'patient18_variables_smoothed.mat'));

% load MAT-file
frames = load(fullfile(directory, 'inflow_outflow_frames.mat'));
% extract inflow and outflow frames
inflow_frames = frames.inflow_frames;
outflow_frames = frames.outflow_frames;

% extract and reshape raw data
vecs_raw = raw.in_struct.vecs;
in_shape = squeeze(raw.in_struct.in_shape);
ntime = in_shape(1);
ny = in_shape(2);
nx = in_shape(3);

vecs_raw = permute(vecs_raw, [2,1,3]);
vecs_raw = reshape(vecs_raw, [nx, ny, ntime, 2]);

% extract and reshape smoothed data
vecs_smooth = smoothed.data.V;
vecs_smooth = permute(vecs_smooth(:,:,:,1:2), [2,3,1,4]);

% compute velocity magnitudes
mag_raw_all = sqrt(vecs_raw(:,:,:,1).^2 + vecs_raw(:,:,:,2).^2);
mag_smooth_all = sqrt(vecs_smooth(:,:,:,1).^2 + vecs_smooth(:,:,:,2).^2);

% compute inflow and outflow averages
mag_inflow_raw = mean(mag_raw_all(:,:,inflow_frames), 3);
mag_inflow_smooth = mean(mag_smooth_all(:,:,inflow_frames), 3);
mag_outflow_raw = mean(mag_raw_all(:,:,outflow_frames), 3);
mag_outflow_smooth = mean(mag_smooth_all(:,:,outflow_frames), 3);

% position grids
pos_raw = raw.in_struct.pos;
pos_raw = reshape(pos_raw, size(vecs_raw,[1,2,4]));
pos_smooth = smoothed.data.grid;
pos_smooth = permute(pos_smooth, [2,3,1]);

% draw inflow line
figure;
pcolor(double(pos_smooth(:,:,1)), double(pos_smooth(:,:,2)), double(mag_inflow_smooth'));
shading flat;
colormap jet; 
axis image;
set(gca, 'YDir', 'reverse');
title('Draw Inflow Line');
h_in = drawline;
in_pos = h_in.Position;

% extract inflow profiles
inflow_raw = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(mag_inflow_raw'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)));

inflow_smooth = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(mag_inflow_smooth'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)));

% draw outflow line
figure;
pcolor(double(pos_smooth(:,:,1)), double(pos_smooth(:,:,2)), double(mag_outflow_smooth'));
shading flat; 
colormap jet; 
axis image;
set(gca, 'YDir', 'reverse');
title('Draw Outflow Line');
h_out = drawline;
out_pos = h_out.Position;

% extract outflow profiles
outflow_raw = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(mag_outflow_raw'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)));

outflow_smooth = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(mag_outflow_smooth'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)));

% plot inflow and outflow profiles
figure;

% inflow profiles
subplot(1,2,1);
plot(inflow_raw, 'r-', 'DisplayName', 'Raw');
hold on;
plot(inflow_smooth, 'b-', 'DisplayName', 'Smoothed');
title('Inflow Region Velocity Magnitude');
xlabel('Line Position');
ylabel('Velocity Magnitude');
legend; grid on;

% outflow profiles
subplot(1,2,2);
plot(outflow_raw, 'r-', 'DisplayName', 'Raw');
hold on;
plot(outflow_smooth, 'b-', 'DisplayName', 'Smoothed');
title('Outflow Region Velocity Magnitude');
xlabel('Line Position');
ylabel('Velocity Magnitude');
legend; grid on;