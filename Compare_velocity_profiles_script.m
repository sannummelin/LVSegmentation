clearvars;
clc;

% path to directory
directory = 'D:\San\LVSegmentation';

% load raw and smoothed data
raw = load(fullfile(directory, 'patient18_variables.mat'));
smoothed = load(fullfile(directory, 'patient18_variables_smoothed.mat'));

% find indices of according frames
inflow_frames = find(strcmp(raw.in_struct.phase, "Inflow"));
outflow_frames = find(strcmp(raw.in_struct.phase, "Outflow"));

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

% path to MAT-file with drawn lines for inflow and outflow
try
draw_lines_file = fullfile(directory, 'draw_lines.mat');
catch
    disp('Could not find draw_lines.mat')
end

if exist(draw_lines_file, 'file')
    % load previously drawn lines
    load(draw_lines_file, 'in_pos', 'out_pos');
    fprintf('Loaded existing inflow/outflow lines from %s\n', draw_lines_file);
else
    % show smooth or raw inflow data (true = smooth, false = raw)
    show_smooth = false;

    % inflow data
    if show_smooth
        X = double(pos_smooth(:,:,1));
        Y = double(pos_smooth(:,:,2));
        C = double(mag_inflow_smooth');
    else
        X = double(pos_raw(:,:,1));
        Y = double(pos_raw(:,:,2));
        C = double(mag_inflow_raw);
    end

    % draw inflow line
    figure;
    pcolor(X, Y, C);
    shading flat;
    colormap jet;
    axis image;
    set(gca, 'YDir', 'reverse');
    title('Draw Inflow Line');
    h_in = drawline;
    in_pos = h_in.Position;

    % outflow data
    if show_smooth
        X = double(pos_smooth(:,:,1));
        Y = double(pos_smooth(:,:,2));
        C = double(mag_outflow_smooth');
    else
        X = double(pos_raw(:,:,1));
        Y = double(pos_raw(:,:,2));
        C = double(mag_outflow_raw);
    end

    % draw outflow line
    figure;
    pcolor(X, Y, C);
    shading flat; 
    colormap jet; 
    axis image;
    set(gca, 'YDir', 'reverse');
    title('Draw Outflow Line');
    h_out = drawline;
    out_pos = h_out.Position;

    % save in_pos and out_pos in MAT-file
    save(fullfile(directory, 'draw_lines.mat'), 'in_pos', 'out_pos');
    fprintf('Saved drawn lines to %s\n', draw_lines_file);
end

% number of sample points along the lines
n_pts = 10;

% extract inflow profiles
inflow_raw = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(mag_inflow_raw'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

inflow_smooth = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(mag_inflow_smooth'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

% extract outflow profiles
outflow_raw = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(mag_outflow_raw'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

outflow_smooth = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(mag_outflow_smooth'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

%%
% plot inflow and outflow profiles
figure;

% calculate points along the drawn inflow line
inflow_line_x = linspace(in_pos(1,1), in_pos(2,1), n_pts);
inflow_line_y = linspace(in_pos(1,2), in_pos(2,2), n_pts);
inflow_length = linspace(0, norm(diff(in_pos)), n_pts);

% inflow profiles
subplot(1,2,1);
plot(inflow_length, inflow_raw, 'r-', 'DisplayName', 'Raw');
hold on;
plot(inflow_length, inflow_smooth, 'b-', 'DisplayName', 'Smoothed');
title('Inflow Region Velocity Magnitude');
xlabel('Line Position');
ylabel('Velocity Magnitude');
legend; grid on;

% calculate points along the drawn outflow line
outflow_line_x = linspace(out_pos(1,1), out_pos(2,1), n_pts);
outflow_line_y = linspace(out_pos(1,2), out_pos(2,2), n_pts);
outflow_length = linspace(0, norm(diff(out_pos)), n_pts);

% outflow profiles
subplot(1,2,2);
plot(outflow_raw, 'r-', 'DisplayName', 'Raw');
hold on;
plot(outflow_smooth, 'b-', 'DisplayName', 'Smoothed');
title('Outflow Region Velocity Magnitude');
xlabel('Line Position');
ylabel('Velocity Magnitude');
legend; grid on;

%%
% normal vector
in_norm = diff([-in_pos(:,2), in_pos(:,1)]);     % inflow
out_norm = diff([-out_pos(:,2), out_pos(:,1)]);  % outflow

% time-resolved inflow raw profile and projection on normal
inflow_tr_raw = zeros(n_pts, 1);
v_proj_n_in_raw = zeros(n_pts, 1);

for i = 1:length(inflow_frames)
    mag = mag_raw_all(:,:,inflow_frames(i));
    inflow_tr_raw(:,i) = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(mag'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

    vx = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(vecs_raw(:,:,i,1)'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

    vy = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(vecs_raw(:,:,i,2)'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

    v_proj_n_in_raw(:,i) = [vx, vy] * in_norm' ./ norm(in_norm);
end

% time-resolved inflow smooth profile and projection on normal
inflow_tr_smooth = zeros(n_pts, 1);
v_proj_n_in_smooth = zeros(n_pts, 1);

for i = 1:length(inflow_frames)
    mag = mag_smooth_all(:,:,inflow_frames(i));
    inflow_tr_smooth(:,i) = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(mag'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

    vx = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(vecs_smooth(:,:,i,1)'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

    vy = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(vecs_smooth(:,:,i,2)'), ...
    double(in_pos(:,1)), ...
    double(in_pos(:,2)), ...
    n_pts);

    v_proj_n_in_smooth(:,i) = [vx, vy] * in_norm' ./ norm(in_norm);
end

% time-resolved outflow raw profile and projection on normal
outflow_tr_raw = zeros(n_pts, 1);
v_proj_n_out_raw = zeros(n_pts, 1);

for i = 1:length(outflow_frames)
    mag = mag_raw_all(:,:,outflow_frames(i));
    outflow_tr_raw(:,i) = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(mag'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

    vx = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(vecs_raw(:,:,i,1)'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

    vy = improfile(...
    double(pos_raw(:,:,1))', ...
    double(pos_raw(:,:,2))', ...
    double(vecs_raw(:,:,i,2)'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

    v_proj_n_out_raw(:,i) = [vx, vy] * out_norm' ./ norm(out_norm);
end

% time-resolved outflow smooth profile and projection on normal
outflow_tr_smooth = zeros(n_pts, 1);
v_proj_n_out_smooth = zeros(n_pts, 1);
for i = 1:length(outflow_frames)
    mag = mag_smooth_all(:,:,outflow_frames(i));
    outflow_tr_smooth(:,i) = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(mag'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

    vx = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(vecs_smooth(:,:,i,1)'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

    vy = improfile(...
    double(pos_smooth(:,:,1))', ...
    double(pos_smooth(:,:,2))', ...
    double(vecs_smooth(:,:,i,2)'), ...
    double(out_pos(:,1)), ...
    double(out_pos(:,2)), ...
    n_pts);

    v_proj_n_out_smooth(:,i) = [vx, vy] * out_norm' ./ norm(out_norm);
end

%%
% tiled plot
figure;
T = tiledlayout(2,2);

% time axis
inflow_time = raw.in_struct.tax(inflow_frames);
outflow_time = raw.in_struct.tax(outflow_frames);

% meshgrids for line length (L) and time (T)
[L_inflow, T_inflow] = meshgrid(inflow_length, inflow_time);
[L_outflow, T_outflow] = meshgrid(outflow_length, outflow_time);

% plot inflow raw
nexttile
surf(L_inflow, T_inflow, inflow_tr_raw', 'EdgeColor', 'none');  % transpose to match dimensions
title('Inflow Raw');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
colorbar;
shading interp;

% plot inflow smooth
nexttile
surf(L_inflow, T_inflow, inflow_tr_smooth', 'EdgeColor', 'none');
title('Inflow Smooth');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
colorbar;
shading interp;

% plot outflow raw
nexttile
surf(L_outflow, T_outflow, outflow_tr_raw', 'EdgeColor', 'none');
title('Outflow Raw');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
colorbar;
shading interp;

% plot outflow smooth
nexttile
surf(L_outflow, T_outflow, outflow_tr_smooth', 'EdgeColor', 'none');
title('Outflow Smooth');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
colorbar;
shading interp;