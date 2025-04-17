% Amanda Syamsul
% March 7th, 2025
% Analysis on IW propagation proximal to OBS

close all; clc
% clear all;

data_dir = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(data_dir);

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');

all_data_2020 = readtable('OIW_2020_data.csv');

OBS_data_2020 = readtable('IW_OBS_2020.csv');
OBS_latband_data_2020 = readtable('IW_OBS_latband_2020.csv');

QC_OBS_data_2020 = readtable('QC_IW_OBS_2020.csv');
QC_OBS_latband_data_2020 = readtable('QC_IW_OBS_latband_2020.csv');
QC_latband_30pt = readtable('QC_IW_OBS_latband_2020_30pt.csv');
fs = 13;

%% Defining datasets
square_around_OBS = QC_OBS_data_2020;
latband_across_OBS = QC_OBS_latband_data_2020;

obs_coords = [21.00116 117.40267];

% Define geographical bounds (longitude and latitude)
lon1 = 116.2; lon2 = 117.7;  % Longitude range of the study area
lat1 = 20.0; lat2 = 21.5;    % Latitude range of the study area

% Define x-y coordinate system bounds (based on satellite image pixels)
x_min = 0; x_max = 500;      % X-coordinates in image space
y_min = -400; y_max = 0;     % Y-coordinates in image space

all_ys = [-400:0];

% Compute the range of coordinates
image_width = x_max - x_min;    % Total width of the image in x-coordinates
image_height = y_max - y_min;   % Total height of the image in y-coordinates

lon_obs = lon1 + (square_around_OBS.x_coord - x_min) / image_width * (lon2 - lon1);
lon_obs_latband = lon1 + (latband_across_OBS.x_coord - x_min) / image_width * (lon2 - lon1);
square_around_OBS.longitude = lon_obs;
latband_across_OBS.longitude = lon_obs_latband;

bathymetry1 = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "elevation");
bathymetry = bathymetry1';
lon = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "lon");
lat = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "lat");

% Find indices within the latitude band
lat_idx = lat >= 21.0012 & lat <= 21.1012;

% Subset bathymetry to only those latitudes
bathymetry_latband = bathymetry(lat_idx, :);

% Now compute the average along latitude dimension (rows)
avg_depth_along_lon = mean(bathymetry_latband, 1, 'omitnan'); 

% % Plot average depth versus longitude
% figure();
% plot(lon, avg_depth_along_lon, 'LineWidth', 2);
% xline(obs_coords(2), 'g-', 'LineWidth',2)
% xline(117.0792445, 'k-', 'LineWidth',2)
% legend('depth', 'OBS', 'Dongsha')
% xlabel('Longitude', 'FontSize', fs)
% ylabel('Average depth (m)', 'FontSize', fs)
% title('Average Depth Along Longitude', 'FontSize', fs, 'FontAngle', 'italic');
% grid on;

%% PART 1: Plot velocity binned by depth

% 1.1 Plot average velocity binned by depth

fig = figure();
subplot(1,2,1)
dataset = latband_across_OBS;

% If only plotting specific times, define the time range
start_date = datetime(2020, 5, 1); 
end_date = datetime(2020, 8, 31); 

% Find indices where datetime falls within the range
range = (dataset.datetime >= start_date) & (dataset.datetime <= end_date);

data_lon = dataset.longitude(range);
data_velocity = dataset.velocity(range);

% Match longitudes to corresponding average depths
idx = knnsearch(lon(:), data_lon(:)); % knnsearch(X,Y) finds the nearest neighbor in X for each point in Y
depths_at_data_lon = avg_depth_along_lon(idx); % Get average depths

% Bin velocities by depth
bin_edges = -1000:25:0; % depth  bins
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[N, edges, bin_indices] = histcounts(depths_at_data_lon, bin_edges);

% Average velocity per depth bin
avg_velocity = arrayfun(@(b) mean(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));
std_velocity = arrayfun(@(b) std(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));

% Define x and y
x = bin_centers;
y = avg_velocity;

% Get rate of slowing in shallow waters (after 900 m depth) & exclude
% anomalous points
valid_idx = (x>-900) & ~isnan(y) & ~isinf(y);% & (y>1) & (y<3.5);
x_clean = x(valid_idx);
y_clean = y(valid_idx);

% Perform polyfit on cleaned data
p = polyfit(x_clean, y_clean, 1);

xline(-900, 'k--', 'LineWidth',2)
hold on;
errorbar(bin_centers, avg_velocity, std_velocity, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
scatter(bin_centers, avg_velocity, 80, 'filled');
yline(3.23, 'b--', 'LineWidth',2)
yline(2.22, 'b--', 'LineWidth',2)
% text(-850, 3.33, 'mean speed in deep basin (Ramp et al., 2010)', 'color', 'b')
% text(-850, 2.33, 'mean speed over cont. slope (Ramp et al., 2010)', 'color', 'b')
% Perform polynomial fit
x_fit = linspace(min(bin_centers), max(bin_centers), 100);
y_fit = polyval(p, x_fit);  % Evaluate polynomial at x_fit

% Plot the fit
plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');

% title(num2str(years(i)), 'FontSize', fs, 'FontAngle', 'italic');
% text(-1100,1.2, ['slope: ' num2str(p(1), '%.5f')], 'color', 'r')
xlim([-800, -300])
ylim([0,5])
grid on;
title('Avg. velocity binned by depth — Latitudinal band over OBS')
subtitle('Latitude 21.0012 - 21.1012')
xlabel('Depth (m)');
ylabel('Velocity (m/s)');
% subtitle('Latitude 21.0012 - 21.1012, Longitude 117.3777 - 117.4777')

all=axes(fig,'visible','off'); 
all.XLabel.Visible='on';
all.YLabel.Visible='on';
xlabel(all, 'Depth (m)', 'FontSize', fs-2)
ylabel(all, 'Average velocity (m/s)', 'FontSize', fs-2)

% 1.2 Plot all velocities binned by depth
ttt = dataset(range, :);

subplot(1,2,2)
hold on;
for b = 1:length(bin_centers)
    % Get all velocity values in the bin
    velocities_in_bin = ttt.velocity(bin_indices == b);
    
    % Scatter plot with jitter
    jitter = (rand(size(velocities_in_bin)) - 0.5) * 5; % Adjust jitter magnitude as needed
    scatter(bin_centers(b) + jitter, velocities_in_bin, 10, 'filled');
    hold on
end

yline(3.23, 'b--', 'LineWidth',2)
yline(2.22, 'b--', 'LineWidth',2)
% text(-850, 3.33, 'mean speed in deep basin (Ramp et al., 2010)', 'color', 'b')
% text(-850, 2.33, 'mean speed over cont. slope (Ramp et al., 2010)', 'color', 'b')

xlabel('Depth (m)');
ylabel('Velocity (m/s)');
title('All velocities binned by depth — Latitudinal band over OBS')
subtitle('Latitude 21.0012 - 21.1012')
ylim([0 5])
grid on;
hold off;

%% PART 2: Correlating velocity with amplitude
filtered_waves = detection_window;
filtered_waves(bad_indices, :) = []; % Removes the row(s)

% filtered_waves = filtered_waves(filtered_waves.amps < 2e5, :); %use this
% if only want a specific amp threshold

% Initialize amplitude column in ttt
ttt.amplitude = NaN(height(ttt), 1);
unique_templates = unique(filtered_waves.template);
max_time_diff = hours(4);  % Max allowed difference for matching

% figure()
% for t = 1:length(unique_templates)
for t=8
    this_template = unique_templates(t);
    filtered_waves_by_temp = filtered_waves(filtered_waves.template == this_template, :);
    
    matched_amps = NaN(height(ttt), 1);
    
    for i = 1:height(ttt)
        time_diffs = abs(filtered_waves_by_temp.detection_time - ttt.datetime(i));
        [min_diff, closest_idx] = min(time_diffs);
    
        if min_diff <= max_time_diff
            matched_amps(i) = filtered_waves_by_temp.amps(closest_idx);
        end
    end

    ttt.amplitude = matched_amps;

    figure();
    hold on;

    for b = 1:length(bin_centers)
        idx = bin_indices == b;
    
        velocities = ttt.velocity(idx);
        amplitudes = ttt.amplitude(idx);
    
        valid = ~isnan(velocities) & ~isnan(amplitudes);
        velocities = velocities(valid);
        amplitudes = amplitudes(valid);
    
        jitter = (rand(size(velocities)) - 0.5) * 5;
        scatter(bin_centers(b) + jitter, velocities, amplitudes * 1e-3, amplitudes, 'filled');
    end

    % Plot reference lines
    yline(3.23, 'b--', 'LineWidth', 2);
    yline(2.22, 'b--', 'LineWidth', 2);

    xlabel('Depth (m)');
    ylabel('Velocity (m/s)');
    title('Velocity binned by depth — All velocities for a given wave');
    % subtitle(sprintf('Template %d — Latitude 21.0012 - 21.1012', this_template));
    xlim([-850 -300]);
    ylim([0 5]);
    colorbar;
    caxis([0 6.5e5]);
    colormap(jet);
    grid on;
    hold off;
end


%% PART 3: Bin 5-point average velocities by 5-point average depth

figure();
hold on;

% Step 1: Smooth with moving average
window_size = 5;
avg_velocity = movmean(ttt.velocity, window_size, 'Endpoints', 'discard');
avg_depth = movmean(depths_at_data_lon, window_size, 'Endpoints', 'discard');

% Step 2: Bin by average depth
bin_edges = -1000:25:0;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[~, ~, bin_indices] = histcounts(avg_depth, bin_edges);

% Step 3: Compute average and std per bin
avg_velocity_by_bin = arrayfun(@(b) mean(avg_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));
std_velocity_by_bin = arrayfun(@(b) std(avg_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));

errorbar(bin_centers, avg_velocity_by_bin, std_velocity_by_bin, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
scatter(bin_centers, avg_velocity_by_bin, 80, 'filled');

yline(3.23, 'b--', 'LineWidth', 2)
yline(2.22, 'b--', 'LineWidth', 2)

text(-1100, 1.2, ['slope: ' num2str(p(1), '%.5f')], 'Color', 'r', 'FontSize', fs-4)

xlim([-800, -300])
ylim([0, 5])
grid on;

xlabel('Depth (m)', 'FontSize', fs)
ylabel('5-point Average Velocity (m/s)', 'FontSize', fs)
title('Smoothed Avg. Velocity Binned by Avg. Depth — Latitudinal Band Over OBS', 'FontSize', fs)
subtitle('Latitude 21.0012 - 21.1012', 'FontSize', fs-2)

%% PART 4: Correlating avg. velocities with amplitudes
% close all;

% Now apply 5-point moving averages
n = height(ttt);
valid_len = n - window_size + 1;

% avg_velocity = NaN(valid_len, 1);
% avg_amplitude = NaN(valid_len, 1);
% avg_depth = NaN(valid_len, 1);

for i = 1:valid_len
    v_window = ttt.velocity(i:i+window_size-1);
    a_window = ttt.amplitude(i:i+window_size-1);
    d_window = depths_at_data_lon(i:i+window_size-1);

    avg_velocity(i) = nanmean(v_window);
    avg_amplitude(i) = nanmean(a_window);
    avg_depth(i) = nanmean(d_window);
end

% Ensure all are column vectors
avg_velocity = avg_velocity(:);
avg_amplitude = avg_amplitude(:);
avg_depth = avg_depth(:);

valid_idx = ~isnan(avg_velocity) & ~isnan(avg_amplitude) & ~isnan(avg_depth);
avg_velocity = avg_velocity(valid_idx);
avg_amplitude = avg_amplitude(valid_idx);
avg_depth = avg_depth(valid_idx);

% Bin by depth
bin_edges = -1000:25:0;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[~, ~, bin_indices] = histcounts(avg_depth, bin_edges);

figure;
hold on;
for b = 1:length(bin_centers)
    idx = bin_indices == b;

    v = avg_velocity(idx);
    a = avg_amplitude(idx);

    % Jitter for depth
    jitter = (rand(size(v)) - 0.5) * 5;

    scatter(bin_centers(b) + jitter, v, log10(a * 1e50), log10(a), 'filled');  % size & color by amp
end

yline(3.23, 'b--', 'LineWidth', 2)
yline(2.22, 'b--', 'LineWidth', 2)

xlabel('Depth (m)');
ylabel('Smoothed Velocity (m/s)');
title('Smoothed Velocity Binned by Depth — Colored by Smoothed Amplitude');
subtitle('5-point windows, Lat 21.0012–21.1012');
xlim([-850 -300]);
ylim([0 5]);
colorbar;
% caxis([0 6.5e5])
colormap(jet);
grid on;
hold off;
