% Amanda Syamsul
% March 7th, 2025
% Analysis on Internal Wave propagation proximal to OBS

close all; clear all; clc

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
all_detections = readtable('all_detections5.csv');
luzon = readtable('luzon.csv');

all_data_2020 = readtable('OIW_2020_data.csv');
OBS_latband_data_2020 = readtable('IW_OBS_latband_2020.csv');
avglat = readtable('avg_OBS_latband_data_2020.csv');
OBS_data_2020 = readtable('IW_OBS_2020.csv'); % square around OBS

% QC_OBS_data_2020 = readtable('QC_IW_OBS_2020.csv');
% QC_OBS_latband_data_2020 = readtable('QC_IW_OBS_latband_2020.csv');

OBS_latband_data_2020_adjusted = readtable('OBS_latband_2020_adjusted.csv');
fs = 13;

% Choose which datasets to use for this analysis
square_around_OBS = OBS_data_2020;
latband_across_OBS = OBS_latband_data_2020_adjusted;

dataset = latband_across_OBS;
plot_year = 2020;
obs_coords = [21.00116 117.40267];

% Define geographical bounds (longitude and latitude)
lon1 = 116.2; lon2 = 117.7;  
lat1 = 20.0; lat2 = 21.5;   

% Define x-y coordinate system bounds (based on satellite image pixels)
x_min = 0; x_max = 500;      % X-coord in image space
y_min = -400; y_max = 0;     % Y-coord in image space

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
avg_depth_along_lon = mean(bathymetry, 1, 'omitnan'); 
avg_depth_along_lon_latband = mean(bathymetry_latband, 1, 'omitnan'); 

% alldet_small = all_detections(:, {'crest_time','crests', 'trough_time', 'troughs', 'period'});  
% writetable(alldet_small, 'OBP_measurements.csv');

%% Manuscript Figure 2b: Plot average depth versus longitude
close all;

figure();
plot(lon, avg_depth_along_lon_latband, 'LineWidth', 2);
hold on
plot(lon, avg_depth_along_lon, 'LineWidth', 2);
xline(obs_coords(2), 'LineWidth',2)
xline(117.0792445, 'k--', 'LineWidth',2)
legend('Latitude 21.0012\circN to 21.1012\circN','Latitude 20.0\circN to 21.5\circN', 'OBP', 'Cutoff longitude')
xlabel('Longitude', 'FontSize', fs)
ylabel('Average depth (m)', 'FontSize', fs)
% title('Average Depth Along Longitude', 'FontSize', fs, 'FontAngle', 'italic');
grid on;

%% Manuscript Figure 4 a,b: T and A over time
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
% ylabel('Period (minutes)')
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
% ylabel('Amplitude (psi)')
legend('Location','northeast')

linkaxes([ax1 ax2], 'x')

%% Manuscript Figure 5: T vs A
close all;clc;

data_to_bin = all_detections.period;
ampli = all_detections.amps;

% Define bin edges
startT = min(data_to_bin);
endT   = max(data_to_bin);
edgesT = startT:minutes(1):endT;

[N,~,binIdx] = histcounts(data_to_bin, edgesT);

amps_avg = accumarray(binIdx((binIdx>0)), ampli((binIdx>0)), [], @mean);
amps_std = accumarray(binIdx(binIdx>0), ampli(binIdx>0), [], @std);
amps_count = accumarray(binIdx(binIdx>0), 1, [], @sum);
amps_sem = amps_std ./ sqrt(amps_count);  % standard error

binDates = edgesT(1:end-1)';  % start of each bin

Ts = minutes(edgesT(1:end-1));

figure()
% unbinned version
Tvals = minutes(data_to_bin);
Avals = ampli;
% Fit power law: log A ~ alpha log T + const
% p_unbin = polyfit(log10(Tvals), log10(Avals), 1);
loglog(Tvals, Avals, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 18);
hold on
p_unbin = polyfit(log10(Tvals), log10(Avals), 1)

hold on
% loglog(Ts, amps_avg, 'bo', 'MarkerFaceColor','b', 'MarkerSize',10)
errorbar(Ts, amps_avg, amps_sem, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 8, ...
         'CapSize', 6, 'LineWidth', 2)  % error bars
xlabel('Period (minutes)')
ylabel('Amplitude (Pa)')
grid on
p = polyfit(log10(Ts), log10(amps_avg), 1)
xfit = linspace(min(Ts), max(Ts), 100);
yfit = 10.^(polyval(p, log10(xfit)));

loglog(xfit, yfit, 'b--','LineWidth',3)

% Reference line with slope -1
hold on
C = 10^(mean(log10(amps_avg)) + mean(log10(Ts)));  % scaling so line passes through center
Aref = C * xfit.^(-1);  % slope = -1 in log-log space

loglog(xfit, Aref, 'r--', 'LineWidth', 3)
text(mean(xfit), C*mean(xfit)^(-1)*1.2, 'slope = -1', 'Color', 'r', 'FontSize', 12)

legend('Unbinned T vs A', 'Binned T vs A','Slope = -1.34')

x_unbin = log10(Tvals(:));
y_unbin = log10(Avals(:));

lm_unbin = fitlm(x_unbin, y_unbin);   % linear model
disp(lm_unbin)

x_bin = log10(Ts(:));
y_bin = log10(amps_avg(:));

lm_bin = fitlm(x_bin, y_bin);
disp(lm_bin)

text(5, 1e2, ['unbinned p-value= ' num2str(lm_unbin.Coefficients.pValue(2))])
text(5, 60, ['binned p-value= ' num2str(lm_bin.Coefficients.pValue(2))])

%% Manuscript Figure 6 a,b,c: Plot fortnightly cycle (with three edge sets)
close all;
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

% OTPS (unchanged)
t_otps = (all_detections.DetectionTime(1):hours(1):all_detections.DetectionTime(end));
z_otps = tmd_predict(tide_model, lat, lon, t_otps, 'h');

startDate1 = min(all_detections.DetectionTime) - days(2);
endDate1   = datetime(2020, 2, 27);
edges1 = startDate1:days(13.66):all_detections.DetectionTime(end);

figure()
ax(1) = subplot(4,1,1); hold on
plot(t_otps, z_otps, 'Color',	[0.6350, 0.0780, 0.1840])
xline(edges1, 'k')
ylabel('Tidal height (m)')

ax(2) = subplot(4,1,2); hold on
scatter(all_detections.DetectionTime, minutes(all_detections.period), 12, 'k', 'filled')
xline(edges1, 'k')
ylabel('Period (mins)')

data_to_bin = minutes(all_detections.period);

% Bin once using the combined edges (left-inclusive to avoid double-counting)
[~,~,binIdx] = histcounts(all_detections.DetectionTime, edges1);
valid = (binIdx > 0) & ~isnan(data_to_bin);
fortnightly_avg = accumarray(binIdx(valid), data_to_bin(valid), [numel(edges1)-1, 1], @mean, NaN);

binDates = edges1(1:end-1)';  % start of each bin

ax(3) = subplot(4,1,3); hold on
scatter(binDates, fortnightly_avg, 28, 'k', 'filled')
plot(binDates, fortnightly_avg, '-')
grid on
ylabel('Average period (mins)')

linkaxes(ax, 'x')


%% PART 1: Plot propagation speed binned by depth
close all;

fig = figure();
subplot(1,2,1)
hold on;
for b = 1:length(bin_centers)
    % Get all propagation speed values in the bin
    velocities_in_bin = ttt_speed(bin_indices == b);
    
    jitter = (rand(size(velocities_in_bin)) - 0.5) * 5; % Adjust jitter magnitude as needed
    scatter(bin_centers(b) + jitter, velocities_in_bin, 10, 'filled');
    hold on
end
xline(-900, 'k--', 'LineWidth',2)
xlabel('Depth (m)');
ylabel('propagation speed (m/s)');
set(gca, "XDir", "reverse")
title('Propagation speeds binned by depth')
subtitle('Latitude 21.0012 - 21.1012')
ylim([0 5])
grid on;
hold off;

subplot(1,2,2)

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

hold on;
errorbar(bin_centers, avg_velocity, std_velocity, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
scatter(bin_centers, avg_velocity, 80, 'filled');
xline(-900, 'k--', 'LineWidth',2)
% title(num2str(years(i)), 'FontSize', fs, 'FontAngle', 'italic');
xlim([-1000, -400])
ylim([0,5])
grid on;
title('Avg. propagation speeds binned by depth')
subtitle('Latitude 21.0012 - 21.1012')
xlabel('Depth (m)');
ylabel('propagation speed (m/s)');

set(gca, "XDir", "reverse")

%% Manuscript Figure S6: Correlating propagation speed with amplitude & period
close all;
clc;

ttt = dataset(range, :);
ttt_speed = ttt.velocity;

ttt.amplitude = NaN(height(ttt), 1);
ttt.period = NaN(height(ttt), 1);
unique_templates = unique(all_detections.template);
max_time_diff = hours(4); 
matched_amps = NaN(height(ttt), 1);
matched_periods = NaN(height(ttt), 1);

figure()
subplot(1,2,1)
add_reference_lines()
for t = 1:length(unique_templates)
% for t=16
    
amps_for_this_template = NaN(height(ttt), 1);
periods_for_this_template = NaN(height(ttt), 1);

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
    
    hold on

    for b = 1:length(bin_centers)
        idx = bin_indices == b;
    
        velocities = ttt_speed(idx);
        amplitudes = amps_for_this_template(idx);
    
        valid = ~isnan(velocities) & ~isnan(amplitudes);
        velocities = velocities(valid);
        amplitudes = amplitudes(valid);
    
        jitter = (rand(size(velocities)) - 0.5) * 5;
        scatter(bin_centers(b) + jitter, velocities, amplitudes/3, amplitudes, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    xlim([-1000 -500]);
    ylim([0 5]);
    % cb = colorbar('southoutside');                 
    % cb.Label.String = 'Amplitude (psi)';
    caxis([100 500]);
    colormap(brewermap([],"YlOrRd"))
    grid on;
    set(gca, 'XDir','reverse')
    hold off;
end

ttt.amplitude = matched_amps;
ttt.period = matched_periods;

% PART 3B: Day-by-day moving average
window_size=5;
% Step 1: Get just the date part (no time)
dates_only = dateshift(ttt.datetime, 'start', 'day');

% Step 2: Get unique dates and group indices
[unique_dates, ~, day_group_ids] = unique(dates_only);

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

    avg_velocity    = [avg_velocity; avg_v_day];
    avg_amplitude   = [avg_amplitude; avg_a_day];
    avg_depth       = [avg_depth; avg_d_day];
end

% Bin by depth
bin_edges = -1000:25:0;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[~, ~, bin_indices] = histcounts(avg_depth, bin_edges);

subplot(1,2,2)
add_reference_lines()
hold on
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
        scatter(bin_centers(b) + jitter, v, a/3, a, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    % xlabel('Depth (m)');
    % ylabel('propagation speed (m/s)');
    % title('5-pt avg. propagation speed binned by Depth');
    % subtitle('Lat 21.0012–21.1012');
    xlim([-1000 -500]);
    ylim([1 3]);
    % cb = colorbar('southoutside');             
    cb.Label.String = 'Amplitude (Pa)'; 
    caxis([100 500]);
    colormap(brewermap([],"YlOrRd"))
    grid on;
    % hold off;
end

function add_reference_lines()
% Mean speed on continental slope from Ramp (2010)
    y0 = 2.22;
    dy = 0.18;
    
    ymin = y0 - dy;
    ymax = y0 + dy;
    
    xmin = -2491;
    xmax = -350;
    
    patch([xmin xmax xmax xmin], [ymin ymin ymax ymax], [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
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

