% Amanda Syamsul
% March 7th, 2025
% Analysis on Internal Wave propagation proximal to OBS

% close all; clear all; clc

data_dir = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(data_dir);

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/data/');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');
addpath('/Users/amandasyamsul/Documents/MATLAB/DrosteEffect-BrewerMap-3.2.5.0')
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

filtered_waves = readtable('filtered_waves.csv');
all_detections = readtable('all_detections4.csv');
luzon = readtable('luzon.csv');

all_data_2020 = readtable('OIW_2020_data.csv');
OBS_latband_data_2020 = readtable('IW_OBS_latband_2020.csv');
avglat = readtable('avg_OBS_latband_data_2020.csv');
OBS_data_2020 = readtable('IW_OBS_2020.csv'); % square around OBS

QC_OBS_data_2020 = readtable('QC_IW_OBS_2020.csv');
QC_OBS_latband_data_2020 = readtable('QC_IW_OBS_latband_2020.csv');

OBS_latband_data_2020_adjusted = readtable('OBS_latband_2020_adjusted.csv');
fs = 13;

% Defining datasets

% Choose which datasets to use for this analysis
square_around_OBS = OBS_data_2020;
latband_across_OBS = OBS_latband_data_2020_adjusted;
dataset = latband_across_OBS;


plot_year = 2020;
obs_coords = [21.00116 117.40267];

% Define geographical bounds (longitude and latitude)
lon1 = 116.2; lon2 = 117.7;  % Longitude range of the study area
lat1 = 20.0; lat2 = 21.5;    % Latitude range of the study area

% Define x-y coordinate system bounds (based on satellite image pixels)
x_min = 0; x_max = 500;      % X-coordinates in image space
y_min = -400; y_max = 0;     % Y-coordinates in image space

% Computte the range of coordinates
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

% Plot average depth versus longitude
% figure();
% plot(lon, avg_depth_along_lon, 'LineWidth', 2);
% xline(obs_coords(2), 'g-', 'LineWidth',2)
% xline(117.0792445, 'k-', 'LineWidth',2)
% legend('depth', 'OBS', 'Dongsha')
% xlabel('Longitude', 'FontSize', fs)
% ylabel('Average depth (m)', 'FontSize', fs)
% title('Average Depth Along Longitude', 'FontSize', fs, 'FontAngle', 'italic');
% grid on;
%% T and A over time
close all;
obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

figure()

ax1 = subplot(2,1,1); hold on

www  = unique(all_detections.Type);     % still a cell array of char
cols = lines(numel(www));               % colors for each type
shapes = {'o','s','d','^','v','>'};     % circle, square, diamond, triangle up, down, right

for k = 1:numel(www)
    idx = strcmp(all_detections.Type, www{k});   % indices for this type
    
    % plot with specific marker and color
    plot(all_detections.DetectionTime(idx), minutes(all_detections.period(idx)), ...
        shapes{mod(k-1, numel(shapes))+1}, ...
        'MarkerFaceColor', cols(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', www{k});
end

grid on
ylabel('Period (minutes)')
legend('Location','northeast')

ax2 = subplot(2,1,2); hold on
for k = 1:numel(www)
    idx = strcmp(all_detections.Type, www{k});
    
    plot(all_detections.DetectionTime(idx), all_detections.amps(idx), ...
        shapes{mod(k-1, numel(shapes))+1}, ...
        'MarkerFaceColor', cols(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', www{k});
end

grid on
ylabel('Amplitude (psi)')
legend('Location','northeast')

linkaxes([ax1 ax2], 'x')


%% speed vs T


%% T vs A
close all;

data_to_bin = all_detections.period;
ampli = all_detections.amps;

% Define bin edges
startT = min(data_to_bin);
endT   = max(data_to_bin);
edgesT = startT:minutes(1):endT;

[N,~,binIdx] = histcounts(data_to_bin, edgesT);

amps_avg = accumarray(binIdx((binIdx>0)), ampli((binIdx>0)), [], @mean);


binDates = edgesT(1:end-1)';  % start of each bin

Ts = minutes(edgesT(1:end-1));


figure()
loglog(Ts, amps_avg, 'b-', 'LineWidth',2)
xlabel('Period (minutes)')
ylabel('Amplitude (psi)')
grid on
hold on
p = polyfit(log10(Ts), log10(amps_avg), 1);
xfit = linspace(min(Ts), max(Ts), 100);
yfit = 10.^(polyval(p, log10(xfit)));
loglog(xfit, yfit, 'r--','LineWidth',2)


% unbinned version
Tvals = minutes(data_to_bin);
Avals = ampli;
% Fit power law: log A ~ alpha log T + const
% p_unbin = polyfit(log10(Tvals), log10(Avals), 1);
loglog(Tvals, Avals,'k.')
legend('binned T vs A', 'slope = -1', 'unbinned data')
hold on
% xfit = logspace(log10(min(Tvals)), log10(max(Tvals)), 200);
% yfit = 10.^(polyval(p_unbin, log10(xfit)));
% loglog(xfit, yfit, 'r--','LineWidth',2)
% title(sprintf('Slope = %.2f', p_unbin(1)))

%% Plot fortnightly cycle (with three edge sets)
close all;
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

% OTPS (unchanged)
t_otps = (all_detections.DetectionTime(1):hours(1):all_detections.DetectionTime(end));
z_otps = tmd_predict(tide_model, lat, lon, t_otps, 'h');

% Your 3 edge vectors (as given)
startDate1 = min(all_detections.DetectionTime) - days(2);
endDate1   = datetime(2020, 2, 27);
edges1 = startDate1:days(13.66):all_detections.DetectionTime(end);

% Quick visual check of edges vs data
figure()
ax(1) = subplot(4,1,1); hold on
plot(t_otps, z_otps)
xline(edges1, 'k')
grid on
ylabel('Tidal height (m)')

ax(2) = subplot(4,1,2); hold on
scatter(all_detections.DetectionTime, minutes(all_detections.period), 12, 'filled')
xline(edges1, 'k')
% xline(edges2, 'r')
% xline(edges3, 'b')
grid on
ylabel('Period (mins)')


% Compute fortnightly averages using ALL 3 edge sets
clc;

data_to_bin = minutes(all_detections.period);

% Bin once using the combined edges (left-inclusive to avoid double-counting)
[~,~,binIdx] = histcounts(all_detections.DetectionTime, edges_all);

valid = (binIdx > 0) & ~isnan(data_to_bin);
% Pre-size so we keep bins even if empty (as NaN)
fortnightly_avg = accumarray(binIdx(valid), data_to_bin(valid), [numel(edges_all)-1, 1], @mean, NaN);

binDates = edges_all(1:end-1)';  % start of each bin

ax(3) = subplot(4,1,3); hold on
scatter(binDates, fortnightly_avg, 28, 'filled')
plot(binDates, fortnightly_avg, '-')
grid on
ylabel('Average period (mins)')

linkaxes(ax, 'x')

%% Boxplot for month
% Extract year-month from DetectionTime
dt = all_detections.DetectionTime;

% Format as "MMM yyyy" so Nov 2019 and Nov 2020 stay separate
month_labels = cellstr(datestr(dt, 'mmm yyyy'));

% --- Boxplot by month for Period ---
figure()

subplot(2,1,1)
boxplot(minutes(all_detections.period), month_labels, ...
    'Colors','k','Symbol','k+','Whisker',1.5);
grid on
title('Period by Month')
ylabel('Period (minutes)')
xlabel('Month')
xtickangle(45)  % rotate labels for readability

% --- Boxplot by month for Amplitude ---
subplot(2,1,2)
boxplot(all_detections.amps, month_labels, ...
    'Colors','k','Symbol','k+','Whisker',1.5);
grid on
title('Amplitude by Month')
ylabel('Amplitude (psi)')
xlabel('Month')
xtickangle(45)


%% Boxplot for month velocity
% Extract year-month from DetectionTime
dt = avglat.datetime;

% Format as "MMM yyyy" so Nov 2019 and Nov 2020 stay separate
month_labels = cellstr(datestr(dt, 'mmm yyyy'));

% --- Boxplot by month for Period ---
figure()

subplot(2,1,1)
boxplot(avglat.average_velocity, month_labels, ...
    'Colors','k','Symbol','k+','Whisker',1.5);
grid on
title('Average velocity by Month')
ylabel('Speed (m/s)')
xlabel('Month')
xtickangle(45)  % rotate labels for readability

% --- Boxplot by month for Amplitude ---
subplot(2,1,2)
boxplot(avglat.average_back_azimuth, month_labels, ...
    'Colors','k','Symbol','k+','Whisker',1.5);
grid on
title('Back Azimuth by Month')
ylabel('Back azimuth (psi)')
xlabel('Month')
xtickangle(45)
%% T vs A
figure()

typeA_idx = all_detections.Type == "A";
typeA2_idx = all_detections.Type == "A2";
typeB_idx = all_detections.Type == "B";

% subplot(3,1,1)
scatter(all_detections.DetectionTime(typeA_idx), all_detections.period(typeA_idx), 'filled');
grid on
ylabel('period (mins)')
hold on
% xlabel('period (minutes)')
title('period over time')
% subplot(3,1,2)
scatter(all_detections.DetectionTime(typeA2_idx), all_detections.period(typeA2_idx), 'filled', 'k');
grid on
ylabel('period (mins)')
legend('A', 'A2')
% xlabel('period (minutes)')
% title('period over time (type A2)')
% subplot(3,1,3)
% scatter(all_detections.DetectionTime(typeB_idx), all_detections.period(typeB_idx), 'filled');
% grid on
% ylabel('period (mins)')
% % xlabel('period (minutes)')
% title('period over time (type B)')
%% PART 1: Plot propagation speed binned by depth

% 1.1 Plot average propagation speed binned by depth

fig = figure();
subplot(1,2,1)

% If only plotting specific times, define the time range
start_date = datetime(plot_year, 5, 1); 
end_date = datetime(plot_year, 8, 31); 
range = (dataset.datetime >= start_date) & (dataset.datetime <= end_date);
data_lon = dataset.longitude(range);
data_velocity = dataset.velocity(range);

% Match longitudes to corresponding average depths
idx = knnsearch(lon(:), data_lon(:)); 
depths_at_data_lon = avg_depth_along_lon(idx); % Get average depths

% Bin velocities by depth
bin_edges = -1000:25:0; % depth bins
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[N, edges, bin_indices] = histcounts(depths_at_data_lon, bin_edges);

% Average velocity per depth bin
avg_velocity = arrayfun(@(b) mean(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));
std_velocity = arrayfun(@(b) std(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));

x = bin_centers;
y = avg_velocity;

% Get rate of slowing in shallow waters (after 900 m depth) & exclude anomalous points
valid_idx = (x>-900) & ~isnan(y) & ~isinf(y);% & (y>1) & (y<3.5);
x_clean = x(valid_idx);
y_clean = y(valid_idx);

% Perform polyfit on cleaned data
p = polyfit(x_clean, y_clean, 1);

xline(-900, 'k--', 'LineWidth',2)
hold on;
errorbar(bin_centers, avg_velocity, std_velocity, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
scatter(bin_centers, avg_velocity, 80, 'filled');
add_reference_lines()

% Perform polynomial fit
x_fit = linspace(min(bin_centers), max(bin_centers), 100);
y_fit = polyval(p, x_fit);  % Evaluate polynomial at x_fit

% Plot the fit
plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Linear Fit');

% title(num2str(years(i)), 'FontSize', fs, 'FontAngle', 'italic');
xlim([-800, -300])
ylim([0,5])
grid on;
title('Avg. speed binned by depth — Latitudinal band over OBS')
subtitle('Latitude 21.0012 - 21.1012')
xlabel('Depth (m)');
ylabel('propagation speed (m/s)');
% subtitle('Latitude 21.0012 - 21.1012, Longitude 117.3777 - 117.4777')

all=axes(fig,'visible','off'); 
all.XLabel.Visible='on';
all.YLabel.Visible='on';
xlabel(all, 'Depth (m)', 'FontSize', fs-2)
ylabel(all, 'Avg. propagation speed (m/s)', 'FontSize', fs-2)

% 1.2 Plot all velocities binned by depth
ttt = dataset(range, :);
ttt_speed = ttt.velocity;

subplot(1,2,2)
hold on;
for b = 1:length(bin_centers)
    % Get all propagation speed values in the bin
    velocities_in_bin = ttt_speed(bin_indices == b);
    
    % Scatter plot with jitter
    jitter = (rand(size(velocities_in_bin)) - 0.5) * 5; % Adjust jitter magnitude as needed
    scatter(bin_centers(b) + jitter, velocities_in_bin, 10, 'filled');
    hold on
end

add_reference_lines()

xlabel('Depth (m)');
ylabel('propagation speed (m/s)');
title('All velocities binned by depth — Latitudinal band over OBS')
subtitle('Latitude 21.0012 - 21.1012')
ylim([0 5])
grid on;
hold off;

%% PART 2: Correlating propagation speed with amplitude & period
close all;
clc;

ttt.amplitude = NaN(height(ttt), 1);
ttt.period = NaN(height(ttt), 1);
unique_templates = unique(all_detections.template);
max_time_diff = hours(4); 
matched_amps = NaN(height(ttt), 1);
matched_periods = NaN(height(ttt), 1);

figure()
subplot(1,2,1)
for t = 1:length(unique_templates)

amps_for_this_template = NaN(height(ttt), 1);
periods_for_this_template = NaN(height(ttt), 1);

% for t=16
    this_template = unique_templates(t);
    filtered_waves_by_temp = all_detections(all_detections.template == this_template, :);

    % matched_amps = NaN(height(ttt), 1);
    
    for i = 1:height(ttt)
        time_diffs = abs(filtered_waves_by_temp.DetectionTime - ttt.datetime(i));
        [min_diff, closest_idx] = min(time_diffs);
    
        if min_diff <= max_time_diff
            matched_amps(i) = filtered_waves_by_temp.amps(closest_idx);
            matched_periods(i) = minutes(filtered_waves_by_temp.period(closest_idx));

            amps_for_this_template(i) = filtered_waves_by_temp.amps(closest_idx);
            periods_for_this_template(i) = minutes(filtered_waves_by_temp.period(closest_idx));

            % matched_amps is a cumulative table of all templates, while 
             % amps_for_this_template is reset every iteration
        
        end
    end

    hold on;

    for b = 1:length(bin_centers)
        idx = bin_indices == b;
    
        velocities = ttt_speed(idx);
        amplitudes = amps_for_this_template(idx);
    
        valid = ~isnan(velocities) & ~isnan(amplitudes);
        velocities = velocities(valid);
        amplitudes = amplitudes(valid);
    
        jitter = (rand(size(velocities)) - 0.5) * 5;
        scatter(bin_centers(b) + jitter, velocities, amplitudes * 1e4/4, amplitudes, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    

    add_reference_lines()
    % xlabel('Depth (m)');
    % ylabel('propagation speed (m/s)');
    % title('Propagation speed binned by depth');
    % subtitle(sprintf('Template %d — Latitude 21.0012 - 21.1012', this_template));
    xlim([-800 -200]);
    ylim([0 5]);
    % cb = colorbar('southoutside');                 
    % cb.Label.String = 'Amplitude (psi)';
    caxis([0 0.07]);
    % colormap(flipud(autumn));
    colormap(brewermap([],"YlOrRd"))
    grid on;
    set(gca, 'XDir','reverse')
    hold off;
end

ttt.amplitude = matched_amps;
ttt.period = matched_periods;

%% PART 3a: Correlating 5-point avg propagation speeds with amplitudes (Choose between 3a and 3b)
% close all;

window_size=5;
% Now apply 5-point moving averages
n = height(ttt);
valid_len = n - window_size + 1;

% Method 1
for i = 1:valid_len
    v_window = ttt_speed(i:i+window_size-1);
    a_window = ttt.amplitude(i:i+window_size-1);
    d_window = depths_at_data_lon(i:i+window_size-1);

    avg_velocity(i) = nanmean(v_window);
    avg_amplitude(i) = nanmean(a_window);
    avg_depth(i) = nanmean(d_window);
end

% Method 2
% avg_velocity = movmean(ttt_speed, window_size, 'Endpoints','discard');
% avg_amplitude = movmean(ttt.amplitude, window_size, 'Endpoints','discard');
% avg_depth = movmean(depths_at_data_lon, window_size, 'Endpoints','discard');
% 
% avg_velocity = avg_velocity(:);
% avg_amplitude = avg_amplitude(:);
% avg_depth = avg_depth(:);
% 
% valid_idx = ~isnan(avg_velocity) & ~isnan(avg_amplitude) & ~isnan(avg_depth);
% avg_velocity = avg_velocity(valid_idx);
% avg_amplitude = avg_amplitude(valid_idx);
% avg_depth = avg_depth(valid_idx);

% Bin by depth
bin_edges = -1000:25:0;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[~, ~, bin_indices] = histcounts(avg_depth, bin_edges);

subplot(1,2,2)
plot_moving_averages(bin_centers,bin_indices,avg_velocity,avg_amplitude)

%% PART 3B: Day-by-day moving average
window_size=5;
% Step 1: Get just the date part (no time)
dates_only = dateshift(ttt.datetime, 'start', 'day');

% Step 2: Get unique dates and group indices
[unique_dates, ~, day_group_ids] = unique(dates_only);

% Initialize storage
avg_velocity    = [];
avg_amplitude   = [];
avg_depth       = [];

% Loop over each day
for d = 1:length(unique_dates)
    idx_today = find(day_group_ids == d);
    
    % Skip if too few points for window
    if length(idx_today) < window_size
        continue
    end
    
    % Extract data for this day
    v_day = ttt_speed(idx_today);
    a_day = ttt.amplitude(idx_today);
    d_day = depths_at_data_lon(idx_today);
    
    valid_len = length(idx_today) - window_size + 1;
    
    avg_v_day = NaN(valid_len,1);
    avg_a_day = NaN(valid_len,1);
    avg_d_day = NaN(valid_len,1);

    for i = 1:valid_len
        v_window = v_day(i:i+window_size-1);
        a_window = a_day(i:i+window_size-1);
        d_window = d_day(i:i+window_size-1);
        
        avg_v_day(i) = nanmean(v_window);
        avg_a_day(i) = nanmean(a_window);
        avg_d_day(i) = nanmean(d_window);
    end

    % Append to full arrays
    avg_velocity    = [avg_velocity; avg_v_day];
    avg_amplitude   = [avg_amplitude; avg_a_day];
    avg_depth       = [avg_depth; avg_d_day];
end

% Bin by depth
bin_edges = -1000:25:0;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[~, ~, bin_indices] = histcounts(avg_depth, bin_edges);

subplot(1,2,2)
plot_moving_averages(bin_centers,bin_indices,avg_velocity,avg_amplitude)
set(gca, 'XDir','reverse')

%%

function plot_moving_averages(bin_centers,bin_indices,avg_velocity,avg_amplitude)
    hold on;
    for b = 6:length(bin_centers)
        idx = bin_indices == b;
    
        v = avg_velocity(idx);
        a = avg_amplitude(idx);
    
        % Jitter for depth
        jitter = (rand(size(v)) - 0.5) * 5;
    
        % scatter(bin_centers(b-window_size) + jitter, v, a * 1e3, a, 'filled', 'MarkerFaceAlpha', 0.7);
        scatter(bin_centers(b) + jitter, v, a * 1e4/4, a, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    add_reference_lines()
    
    % xlabel('Depth (m)');
    % ylabel('propagation speed (m/s)');
    % title('5-pt avg. propagation speed binned by Depth');
    % subtitle('Lat 21.0012–21.1012');
    xlim([-800 -200]);
    ylim([0.5 3.5]);
    % cb = colorbar('southoutside');             
    cb.Label.String = 'Amplitude (psi)'; 
    caxis([0 0.07]);
    colormap(brewermap([],"YlOrRd"))
    grid on;
    % hold off;
end

function add_reference_lines()
    yline(3.23, 'b--', 'LineWidth',2);
    yline(2.22, 'b--', 'LineWidth',2);
    % text(-850, 3.33, 'mean speed in deep basin (Ramp et al., 2010)', 'color', 'b')
    % text(-850, 2.33, 'mean speed over cont. slope (Ramp et al., 2010)', 'color', 'b')
end

function bubblelegend(sizeVals, labels, color, location)
% BUBBLELEGEND Create a size legend for scatter/bubble plots.
%
% bubblelegend(sizeVals, labels, color, location)
%   sizeVals - vector of marker sizes (same units as scatter SizeData)
%   labels   - cell array of text labels (same length as sizeVals)
%   color    - RGB triplet or color name (e.g., 'k', [0 0 0])
%   location - legend location string ('northwest','northeast', etc.)
%
% Example:
%   bubblelegend([20 40 60], {'Small','Medium','Large'}, 'k', 'northeast')

    hold on
    ax = gca;
    xL = xlim(ax);
    yL = ylim(ax);

    % Position bubbles in top right corner by default
    xPos = xL(2) + 0.02 * range(xL);
    yPosStart = yL(2) - 0.05 * range(yL);

    % Plot each bubble & label
    for i = 1:length(sizeVals)
        yPos = yPosStart - (i-1) * 0.07 * range(yL);
        scatter(xPos, yPos, sizeVals(i), color, 'filled');
        text(xPos + 0.03 * range(xL), yPos, labels{i}, ...
            'VerticalAlignment','middle');
    end

    % Expand axis so legend is visible
    ax.XLim = [xL(1) xL(2) + 0.2*range(xL)];
end

