% Amanda Syamsul
% April 11th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc;

data_dir = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(data_dir);
%%
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/QC_may-aug2020/');

%% setup filename
% !ls -1 0*svg > Filenames

% Read in CSV files
luzon = readtable('luzon.csv');

all_data_2015 = readtable('OIW_2015_data.csv');
avg_data_2015 = readtable('avg_OIW_2015_data.csv');

all_data_2016 = readtable('OIW_2016_data.csv');
avg_data_2016 = readtable('avg_OIW_2016_data.csv');

all_data_2017 = readtable('OIW_2017_data.csv');
avg_data_2017 = readtable('avg_OIW_2017_data.csv');

all_data_2018 = readtable('OIW_2018_data.csv');
avg_data_2018 = readtable('avg_OIW_2018_data.csv');

all_data_2019 = readtable('OIW_2019_data.csv');
avg_data_2019 = readtable('avg_OIW_2019_data.csv');

all_data_2020 = readtable('OIW_2020_data.csv');
avg_data_2020 = readtable('avg_OIW_2020_data.csv');

all_data_2021 = readtable('OIW_2021_data.csv');
avg_data_2021 = readtable('avg_OIW_2021_data.csv');

all_data_2022 = readtable('OIW_2022_data.csv');
avg_data_2022 = readtable('avg_OIW_2022_data.csv');

all_data_2023 = readtable('OIW_2023_data.csv');
avg_data_2023 = readtable('avg_OIW_2023_data.csv');

all_data_2024 = readtable('OIW_2024_data.csv');
avg_data_2024 = readtable('avg_OIW_2024_data.csv');

fs = 13;
%% Data calculation
% close all;
% 
% [velocity,velocity_std_dev,velocity_std_err,back_azimuth,t_curr_array, mid_x]=calculate_velocity('QC_2020', true, true, 0.9);
% 
% % Saving data to table & SVG 
% 
% col_names = {'datetime', 'velocity','velocity_std_dev', 'back_azimuth', 'x_coord'};
% T = table(t_curr_array', velocity', velocity_std_dev', back_azimuth', mid_x', 'VariableNames', col_names);
% writetable(T,'QC_IW_OBS_latband_2020_30pt.csv')

% Get table of averages
% avg_table = get_averages(T);
% writetable(avg_table, 'avg_OIW_2018_data.csv');

%% Defining datasets

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

% Years for which data is processed
years = 2015:2024;

for i = 1:length(years)
    year = years(i);
    
    % Extract data for the given year using eval
    data_all{i} = eval(sprintf('all_data_%d', year));
    data_avg{i} = eval(sprintf('avg_data_%d', year));

    % Convert datetime to numeric days since January 1st of the year
    days_all{i} = datenum(data_all{i}.datetime) - datenum(year, 1, 1);
    days_avg{i} = datenum(data_avg{i}.datetime) - datenum(year, 1, 1);

    % Compute time of day in decimal hours
    timeofday{i} = hour(data_all{i}.datetime) + minute(data_all{i}.datetime) / 60;

    % Convert x-coordinates to longitudes
    longitude{i} = lon1 + (data_all{i}.x_coord - x_min) / image_width * (lon2 - lon1);
    latitudes{i} = lat1 + (all_ys - y_min) / image_height * (lat2 - lat1);

    % Add longitude as a new column in data_all{i}
    data_all{i}.longitude = longitude{i};
end

% Conversion equation:
% x_coord = x_min + ((longitude - lon1) / (lon2 - lon1)) * image_width;
% y_coord = y_min + ((latitude - lat1) / (lat2 - lat1)) * image_height;

bathymetry1 = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "elevation");
bathymetry = bathymetry1';
lon = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "lon");
lat = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "lat");

% Compute average depth along each longitude point
avg_depth_along_lon = mean(bathymetry, 1, 'omitnan'); % Average along latitude dimension

for i = 1:length(years)
    sorted_xes = sort(unique(data_all{i}.x_coord) );
    minxes(i) = sorted_xes(2); % ignore i=1 because we don't want x=0
end

min_xcoord = min(minxes);
minlon = lon1 + (min_xcoord - x_min) / image_width * (lon2 - lon1);

%% Bathymetry figures
close all;

figure()
contourf(lon, lat, bathymetry, 50, 'LineStyle', 'none');
hold on
plot(obs_coords(2), obs_coords(1), 'p', 'MarkerFaceColor', 'green', 'MarkerSize', 15);

% OBS bounds box
space = 0.05;
obs_x1 = obs_coords(2) - space/2;
obs_x2 = obs_coords(2) + space*(3/2);
obs_y1 = obs_coords(1);
obs_y2 = obs_coords(1) + space*2;

x = [obs_x1, obs_x2, obs_x2, obs_x1, obs_x1];
y = [obs_y1, obs_y1, obs_y2, obs_y2, obs_y1];
plot(x, y, 'k-', 'LineWidth', 3);

xline(117.0792445 ,'k--', 'LineWidth',2)
colormap('turbo');
colorbar;
grid on;
legend('bathymetry', 'OBS', 'OBS bounds', 'data cutoff')
legend show;
title('Bathymetric Contour Map','FontSize', fs, 'FontAngle', 'italic');
xlabel('Longitude', 'FontSize', fs)
ylabel('Latitude', 'FontSize', fs)


% Plot average depth versus longitude
figure();
plot(lon, avg_depth_along_lon, 'LineWidth', 2);
xline(obs_coords(2), 'g-', 'LineWidth',2)
xline(117.0792445, 'k-', 'LineWidth',2)
legend('depth', 'OBS', 'Dongsha')
xlabel('Longitude', 'FontSize', fs)
ylabel('Average depth (m)', 'FontSize', fs)
title('Average Depth Along Longitude', 'FontSize', fs, 'FontAngle', 'italic');
grid on;
%% Plot average velocity binned by depth
close all;
cmap = autumn(length(years));

fig = figure();
num_years = length(years);

for i = 1:num_years
    % Subplot grid (adjust rows/cols depending on number of years)
    subplot(2, ceil(num_years/2), i)
    
    % Time filter
    start_date = datetime(years(i), 5, 1); 
    end_date = datetime(years(i), 8, 31); 
    range = (data_all{i}.datetime >= start_date) & (data_all{i}.datetime <= end_date);

    data_lon = longitude{i}(range);
    data_velocity = data_all{i}.velocity(range);
    
    % Match longitudes to bathymetry
    idx = knnsearch(lon(:), data_lon(:));
    depths_at_data_lon = avg_depth_along_lon(idx);

    % Bin data by depth
    bin_edges = -2000:50:0;
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    [~, ~, bin_indices] = histcounts(depths_at_data_lon, bin_edges);

    % Compute stats
    avg_velocity = arrayfun(@(b) mean(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));
    std_velocity = arrayfun(@(b) std(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));

    % Fit and clean
    x = bin_centers;
    y = avg_velocity;
    valid_idx = (x > -900) & ~isnan(y) & (y > 1) & (y < 3.5);
    x_clean = x(valid_idx);
    y_clean = y(valid_idx);

    p = polyfit(x_clean, y_clean, 1);
    slopes(i) = p(1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);

    % Plot
    hold on;
    xline(-900, 'k--', 'LineWidth', 2)
    errorbar(bin_centers, avg_velocity, std_velocity, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
    scatter(bin_centers, avg_velocity, 80, cmap(i, :), 'filled', 'MarkerEdgeColor', 'k');
    plot(x_fit, y_fit, 'r-', 'LineWidth', 1.5);
    yline(3.23, 'b--', 'LineWidth', 2)
    yline(2.22, 'b--', 'LineWidth', 2)

    title(num2str(years(i)), 'FontSize', fs, 'FontAngle', 'italic');
    text(-1100, 1.2, ['slope: ' num2str(p(1), '%.5f')], 'color', 'r', 'FontSize', fs-4)
    xlim([-1200, -400])
    ylim([0, 3.5])
    grid on;

    % Bootstrapping
    f = y_clean;
    b = x_clean(:);
    a = ones(length(f), 1);
    Nbs = 1000;

    c_store1 = zeros(Nbs, 1);
    c_store2 = zeros(Nbs, 1);

    for k = 1:Nbs
        I = randi(length(f), 1, length(f));
        f_twiddle = f(I);
        b_twiddle = b(I);
        G_twiddle = [ones(length(I), 1), b_twiddle];
        c_twiddle = (G_twiddle' * G_twiddle) \ (G_twiddle' * f_twiddle(:));
        c_store1(k) = c_twiddle(2);
        c_store2(k) = c_twiddle(1);
    end

    stdev_slope(i) = std(c_store1);
    stdev_int(i) = std(c_store2);
end

all=axes(fig,'visible','off'); 
all.XLabel.Visible='on';
all.YLabel.Visible='on';
xlabel(all, 'Depth (m)', 'FontSize', fs-2)
ylabel(all, 'Average velocity (m/s)', 'FontSize', fs-2)

%% Slope of velocity as a function of depth
figure()
scatter(years, slopes, 50, 'filled');
hold on
errorbar(years, slopes, stdev_slope, 'k.', 'HandleVisibility', 'off');
grid on
axis ij
xlim([2014,2025])
title('Rate of Change of Velocity with Depth Over Time', 'FontSize', fs, 'FontAngle', 'italic')
ylabel('slope (m/s per m)', 'FontSize', fs-2)
xlabel('time', 'FontSize', fs-2)

%% Back Azimuth (2015-2024)

fig2 = figure()
colormap('autumn')
for i = 1:length(years)
    subplot(2, 5, i)
    % Error bar plot
    errorbar(data_avg{i}.datetime, data_avg{i}.median_back_azimuth, ...
        data_avg{i}.back_azimuth_std_dev, 'w.', 'CapSize', 0);
    hold on;
    
    % Scatter plot
    scatter(data_avg{i}.datetime, data_avg{i}.median_back_azimuth, ...
        55, days_avg{i}, 'filled');
    
    % Filter out NaN values for polynomial fit
    x = days_avg{i}; 
    y = data_avg{i}.median_back_azimuth; 
    t = data_avg{i}.datetime; 
    valid_indices = ~isnan(y); 
    x_valid = x(valid_indices);
    y_valid = y(valid_indices);
    t_valid = t(valid_indices); 

    % Perform polynomial fit
    p = polyfit(x_valid, y_valid, 3);

    % Generate fitted values over datetime range
    t_fit = linspace(min(t_valid), max(t_valid), 100); 
    x_fit = interp1(t_valid, x_valid, t_fit, 'linear', 'extrap'); 
    y_fit = polyval(p, x_fit); 

    % Plot the fit
    plot(t_fit, y_fit, 'w', 'LineWidth', 1.5);
    
    ylim([85 116]);
    axis ij;
    title(num2str(years(i)), 'FontSize', fs, 'FontAngle', 'italic');
    grid on;

    % Set axes background to dark gray
    set(gca, 'Color', [0.1 0.1 0.1]);
    set(gca, 'GridColor', 'w');
    
    time_x = [datetime(years(i), 5, 1, 3, 0, 0), datetime(years(i), 8, 31, 10, 0, 0)];
    xlim(time_x);

    caxis([120 240]);
end

% sgtitle('Daily Medians of Back Azimuth in South China Sea')
c = colorbar('Position', [0.92 0.11 0.02 0.815]);
% c.Label.String = 'Day of Year';
c.Ticks = [120, 150, 180, 210, 240];  % Monthly ticks
c.TickLabels = {'May', 'June', 'July', 'Aug', 'Sept'};
c.FontSize = fs;

all=axes(fig2,'visible','off'); 
all.XLabel.Visible='on';
all.YLabel.Visible='on';
xlabel(all, 'Time', 'FontSize', fs)
ylabel(all, 'Back azimuth (m/s)', 'FontSize', fs)
sgtitle('Daily Medians of Back Azimuth in South China Sea', 'FontSize', fs+4, 'FontAngle', 'italic')


%% Velocity vs Back Azimuth (day of year)

close all;
figure()
colormap('autumn')

for i = 1:length(years)
    subplot(2, 5, i)

    scatter( data_avg{i}.average_velocity,  data_avg{i}.average_back_azimuth, ...
        55, days_avg{i}, 'filled')
    % errorbar(data_avg{i}.average_velocity,  data_avg{i}.average_back_azimuth, ...
    %      data_avg{i}.average_velocity,  data_avg{i}.average_back_azimuth, 'k.', 'CapSize', 0)
    ylim([80 120])
    xlabel('velocity');
    ylabel('back azimuth');
    axis ij
    title(num2str(years(i)), FontSize=14, FontAngle = 'italic');
    grid on

    % set(gcf, 'Color', 'k');
    set(gca, 'Color', [0.1 0.1 0.1]);
    % set(gca, 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w');
    set(gca, 'GridColor', 'w');
end

% sgtitle('Daily medians of internal wave velocity against back azimuth in South China Sea')

caxis([120 240])


c = colorbar('Position', [0.92 0.11 0.02 0.815]); 
c.Label.String = 'Day of Year';
c.Ticks = [120, 150, 180, 210, 240];  % Monthly ticks
c.TickLabels = {'May', 'June', 'July', 'Aug', 'Sept'};
% c.Label.Color = 'w'; 
% c.Color = 'w';      

%% Velocity vs Back Azimuth (time of day)

close all;
figure()
colormap('hot')

for i = 1:length(years)
    subplot(2, 5, i)
    timeofday_thisyear=timeofday{i};
    III=(timeofday_thisyear>4)&(timeofday_thisyear<6);

    scatter( data_all{i}.velocity(III),  data_all{i}.back_azimuth(III), 55, days_all{i}(III), 'filled')
    ylim([80 120])
    xlabel('velocity');
    ylabel('back azimuth');
    axis ij
    title(num2str(years(i)), FontSize=14, FontAngle = 'italic');
    grid on
end

sgtitle('Daily medians of internal wave velocity against back azimuth in South China Sea')

% Set the color axis to ensure full coverage of tick labels
% caxis([120 240])

c = colorbar('Position', [0.92 0.11 0.02 0.815]); 
c.Label.String = 'time of day';
c.Ticks = [4, 5, 6, 7, 8, 9];  % Monthly ticks
c.TickLabels = {'4 pm', '5 pm', '6 pm', '7 pm', '8 pm', '9 pm'};

%% Individual year Back azimuth circle plots
close all;

figure()
plot_azimuth2(all_data_2020, 'Direction of wave source',  'May to September 2020')

figure()
plot_azimuth2(all_data_2019, 'Direction of wave source',  'May to September 2019')

%% Combined BA circle plots
close all;

figure()
% Convert back azimuth to radians for polar plot
azimuths_rad1 = deg2rad(all_data_2019.back_azimuth);
azimuths_rad2 = deg2rad(all_data_2020.back_azimuth);
azimuths_rad3 = deg2rad(all_data_2021.back_azimuth);

% Scale normalized time to radius range
radius1 = all_data_2019.velocity;
radius2 = all_data_2020.velocity;
radius3 = all_data_2021.velocity;

p = polaraxes; 
hold on;

% Adjust the polar plot to set 0 degrees to North and clockwise direction
p.ThetaDir = 'clockwise';
p.ThetaZeroLocation = 'top';

% Create polar scatter plot with velocity as radius using scatter
scatter(p, azimuths_rad1, radius1, 30, 'filled', 'DisplayName', '2019');
hold on
scatter(p, azimuths_rad2, radius2, 30, 'filled', 'DisplayName', '2020');
scatter(p, azimuths_rad3, radius3, 30, 'filled', 'DisplayName', '2021');

% Legend to show radius represents velocity
% h = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
% h2 = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
% 
% title(h, 'Radius = Velocity (m/s)');
% title(h2, 'Radius = Velocity (m/s)');

title('Direction of wave source',  'April-Sept 2019-2021', FontSize=14, FontAngle = 'italic');
legend show;
hold off;

%%
function plot_azimuth2(data,label, subtitle)
    % Convert back azimuth to radians for polar plot
    azimuths_rad = deg2rad(data.back_azimuth);
    
    % Scale normalized time to desired radius range
    radius = data.velocity;
    
    p = polaraxes; 
    hold on;

    % Adjust the polar plot to set 0 degrees to North and clockwise direction
    p.ThetaDir = 'clockwise';
    p.ThetaZeroLocation = 'top';

    scatter(p, azimuths_rad, radius, 30, 'filled');
    
    h = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
    title(h, 'Radius = Velocity (m/s)');
    title(label, subtitle, FontSize=14, FontAngle = 'italic');
    legend show;
    hold off;
end


function plot_azimuth(data, label, subtitle)
    scatter(data.datetime, data.back_azimuth, 20, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
    grid on;
    title(label, subtitle, FontSize=14, FontAngle = 'italic');
    xlabel('Time');
    ylabel('Back Azimuth (degrees)');
    hold off;
end
