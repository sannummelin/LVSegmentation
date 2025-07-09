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

frames = 1:size(vecs_raw);
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
    load(draw_lines_file, 'in_pos', 'out_pos');
    fprintf('Loaded existing inflow/outflow lines from %s\n', draw_lines_file);
catch
    disp('Could not find draw_lines.mat')
    in_pos = [];
    out_pos = [];
end

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
if isempty(in_pos)
    h_in = drawline();
    in_pos = h_in.Position;
else
    h_in = drawline('Position', in_pos);
    in_pos = h_in.Position;
end

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
if isempty(out_pos)
    h_out = drawline;
    out_pos = h_out.Position;
else
    h_out = drawline('Position', out_pos);
    out_pos = h_out.Position;
end

%% 
% save in_pos and out_pos in MAT-file
save(fullfile(directory, 'draw_lines.mat'), 'in_pos', 'out_pos');
fprintf('Saved drawn lines to %s\n', draw_lines_file);

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
% plot inflow and outflow velocity magnitude profiles
figure;

% calculate points along the drawn inflow line
inflow_line_x = linspace(in_pos(1,1), in_pos(2,1), n_pts);
inflow_line_y = linspace(in_pos(1,2), in_pos(2,2), n_pts);
inflow_length = linspace(0, norm(diff(in_pos)), n_pts);

% inflow velocity magnitude profiles
subplot(1,2,1);
plot(inflow_length, inflow_raw, 'r-', 'DisplayName', 'Raw');
hold on;
plot(inflow_length, inflow_smooth, 'b-', 'DisplayName', 'Smoothed');
title('Inflow Region Velocity Magnitude');
xlabel('Length (mm)');
ylabel('Velocity Magnitude (m/s)');
legend; 
grid on;

% calculate points along the drawn outflow line
outflow_line_x = linspace(out_pos(1,1), out_pos(2,1), n_pts);
outflow_line_y = linspace(out_pos(1,2), out_pos(2,2), n_pts);
outflow_length = linspace(0, norm(diff(out_pos)), n_pts);

% outflow velocity magnitude profiles
subplot(1,2,2);
plot(outflow_raw, 'r-', 'DisplayName', 'Raw');
hold on;
plot(outflow_smooth, 'b-', 'DisplayName', 'Smoothed');
title('Outflow Region Velocity Magnitude');
xlabel('Length (mm)');
ylabel('Velocity Magnitude (m/s)');
legend; 
grid on;

% RMS (root mean squared)
% Bereken RMS fout tussen raw en smoothed
rms_inflow = sqrt(mean((inflow_raw - inflow_smooth).^2, 'omitnan'));
rms_outflow = sqrt(mean((outflow_raw - outflow_smooth).^2, 'omitnan'));

% Print resultaten
fprintf('RMS error inflow: %.4f m/s\n', rms_inflow);
fprintf('RMS error outflow: %.4f m/s\n', rms_outflow);

%%
% normal vector
in_norm = diff([-in_pos(:,2), in_pos(:,1)]);     % inflow
out_norm = diff([-out_pos(:,2), out_pos(:,1)]);  % outflow

% time-resolved inflow raw profile and projection on normal
inflow_tr_raw = zeros(n_pts, 1);
v_proj_n_in_raw = zeros(n_pts, 1);

for i = 1:length(frames)
    mag = mag_raw_all(:,:,frames(i));
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

for i = 1:length(frames)
    mag = mag_smooth_all(:,:,frames(i));
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
    fprintf('max = %.2f, %.2f m/s | magmax = %.2f | projmax = %.2f\n', max(abs(vx)), max(abs(vy)), max(inflow_tr_smooth(:,i)), max(abs(v_proj_n_in_smooth(:,i))))
end

% time-resolved outflow raw profile and projection on normal
outflow_tr_raw = zeros(n_pts, 1);
v_proj_n_out_raw = zeros(n_pts, 1);

for i = 1:length(frames)
    mag = mag_raw_all(:,:,frames(i));
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
for i = 1:length(frames)
    mag = mag_smooth_all(:,:,frames(i));
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
% tiled plot of time-resolved profiles and projection on normal
figure;
T = tiledlayout(2,2);

% time axis
inflow_time = raw.in_struct.tax(inflow_frames);
outflow_time = raw.in_struct.tax(outflow_frames);

% meshgrids for line length (L) and time (T)
[L_inflow, T_inflow] = meshgrid(inflow_length.*1e3, raw.in_struct.tax);
[L_outflow, T_outflow] = meshgrid(outflow_length.*1e3, raw.in_struct.tax);

% color limits
max_inflow = max(max(v_proj_n_in_raw(:), [], "omitmissing"), max(v_proj_n_in_smooth(:), [], "omitmissing"));
max_outflow = max(max(v_proj_n_out_raw(:), [], "omitmissing"), max(v_proj_n_out_smooth(:), [], "omitmissing"));
min_inflow = min(min(v_proj_n_in_raw(:), [], "omitmissing"), min(v_proj_n_in_smooth(:), [], "omitmissing"));
min_outflow = min(min(v_proj_n_out_raw(:), [], "omitmissing"), min(v_proj_n_out_smooth(:), [], "omitmissing"));

% plot inflow raw
nexttile
surf(L_inflow, T_inflow, v_proj_n_in_raw', 'EdgeColor', 'none');  % transpose to match dimensions
title('Inflow Raw');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Normal (m/s)';
shading interp;
clim([min_inflow, max_inflow])
axis tight

% plot inflow smooth
nexttile
surf(L_inflow, T_inflow, v_proj_n_in_smooth', 'EdgeColor', 'none');
title('Inflow Smooth');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Normal (m/s)';
shading interp;
clim([min_inflow, max_inflow])
axis tight

% plot outflow raw
nexttile
surf(L_outflow, T_outflow, v_proj_n_out_raw', 'EdgeColor', 'none');
title('Outflow Raw');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Normal (m/s)';
shading interp;
clim([min_outflow, max_outflow])
axis tight

% plot outflow smooth
nexttile
surf(L_outflow, T_outflow, v_proj_n_out_smooth', 'EdgeColor', 'none');
title('Outflow Smooth');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Normal (m/s)';
shading interp;
clim([min_outflow, max_outflow])
axis tight

%%  
% tiled plot of velocity magnitude
figure;
T = tiledlayout(2,2);

% color limits
max_outflow = max(max(outflow_tr_raw(:), [], "omitmissing"), max(outflow_tr_smooth(:), [], "omitmissing"));
max_inflow = max(max(inflow_tr_raw(:), [], "omitmissing"), max(inflow_tr_smooth(:), [], "omitmissing"));

% plot inflow raw magnitude
nexttile
surf(L_inflow, T_inflow, inflow_tr_raw', 'EdgeColor', 'none');  % transpose to match dimensions
title('Inflow Raw Magnitude');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity Magnitude (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Magnitude (m/s)';
shading interp;
clim([0, max_inflow])
axis tight

% plot inflow smooth magnitude
nexttile
surf(L_inflow, T_inflow, inflow_tr_smooth', 'EdgeColor', 'none');
title('Inflow Smooth Magnitude');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity Magnitude (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Magnitude (m/s)';
shading interp;
clim([0, max_inflow])
axis tight

% plot outflow raw magnitude
nexttile
surf(L_outflow, T_outflow, outflow_tr_raw', 'EdgeColor', 'none');
title('Outflow Raw Magnitude');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity Magnitude (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Magnitude (m/s)';
shading interp;
clim([0, max_outflow])
axis tight

% plot outflow smooth magnitude
nexttile
surf(L_outflow, T_outflow, outflow_tr_smooth', 'EdgeColor', 'none');
title('Outflow Smooth Magnitude');
xlabel('Length (mm)');
ylabel('Time (s)');
zlabel('Velocity Magnitude (m/s)');
view(2);
cb = colorbar();
cb.Label.String = 'Velocity Magnitude (m/s)';
shading interp;
clim([0, max_outflow])
axis tight