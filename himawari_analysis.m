% Amanda Syamsul
% April 11th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc;

data_dir = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(data_dir);


addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/data');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/QC_may-aug2020/');

% setup filename
% !ls -1 0*svg > Filenames

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

fs = 14;
% Data calculation
% close all;
% 
% [velocity,velocity_std_dev,velocity_std_err,back_azimuth,t_curr_array, mid_x]=calculate_velocity('all_2024_data', true, false, 0.9);

% % Saving data to table & SVG 
% 
% col_names = {'datetime', 'velocity','velocity_std_dev', 'back_azimuth', 'x_coord'};
% T = table(t_curr_array', velocity', velocity_std_dev', back_azimuth', mid_x', 'VariableNames', col_names);
% writetable(T,'IW_OBS_latband_2024.csv')

% % Get table of averages
% avg_table = get_averages(T);
% writetable(avg_table, 'avg_OIW_2018_data.csv');

% Defining datasets

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
plot(obs_coords(2),obs_coords(1), 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','blue','Marker', '^','MarkerSize',10);

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
colorbar('southoutside');
grid on;
legend('bathymetry', 'OBS', 'OBS bounds', 'data cutoff')
legend show;
title('Bathymetric Contour Map (21.5N to 20.0S; 116.2W 117.7E)','FontSize', fs);
% xlabel('Longitude', 'FontSize', fs)
% ylabel('Latitude', 'FontSize', fs)


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

%%
close all;

% Target latitude
target_lat = 20.95;

% Find index of closest latitude
[~, lat_idx] = min(abs(lat - target_lat));

% Extract bathymetry along that latitude
bathy_profile = bathymetry1(:, lat_idx);

% Plot
figure()
plot(lon, bathy_profile, 'k-', 'LineWidth', 1.5)
hold on
xline(117.36)
xline(117.427)
slope = (-97+615)/(117.36-117.427);
% Choose coordinates where you want the text to appear
xpos = mean(lon);   % midpoint of longitude range
ypos = mean(bathy_profile);  % somewhere around average depth

% Add text
text(xpos, ypos, ['slope = ' num2str(slope, '%.3f')])

ylabel('Depth (m)')
title(['Bathymetry along latitude ', num2str(lat(lat_idx))])
grid on

%% Plot average velocity binned by depth
close all;
cmap = autumn(length(years));

fig = figure();
num_years = length(years);

for i = 1:num_years
% for i = 9
    % Subplot grid (adjust rows/cols depending on number of years)
    % subplot(2, ceil(num_years/2), i)
    subplot(5, 2, i)
    
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
    valid_idx = (x > -900) & ~isnan(y) & (y > 0) & (y < 3.5);
    x_clean = x(valid_idx);
    y_clean = y(valid_idx);

    p = polyfit(x_clean, y_clean, 1);
    slopes(i) = p(1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);

    % Plot
    hold on;
    add_reference_lines()
    plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5);
    xline(-900, 'k--', 'LineWidth', 2)
    errorbar(bin_centers, avg_velocity, std_velocity, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
    scatter(bin_centers, avg_velocity, 80, cmap(i, :), 'filled', 'MarkerEdgeColor', 'k');

    title(num2str(years(i)), 'FontSize', fs, 'FontName','.AppleSystemUIFont');
    text(-500, 0.7, ['slope: ' num2str(p(1), '%.5f')], 'color', 'k', 'FontSize', fs-2)
    xlim([-1200, -400])
    ylim([0, 3.5])
    set(gca, "XDir", "reverse")
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
% xlabel(all, 'Depth (m)', 'FontSize', fs-2)
% ylabel(all, 'Average propagation speed (m/s)', 'FontSize', fs-2)

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

close all;

% ---- Style defaults ----
fs   = 10;                      % base font size
lw   = 0.8;                     % errorbar line width
ms   = 16;                      % scatter marker size
yl   = [85 116];                % y-limits for back azimuth
mo_ticks = 5:8;                 % May–Aug tick months

% ---- Figure + tiling ----
fig2 = figure('Color','w','Units','inches','Position',[0 0 7.0 9.0]); %#ok<NASGU>
tlo  = tiledlayout(5,2,'TileSpacing','compact','Padding','compact');

colormap('autumn')

for i = 1:numel(years)
    ax = nexttile; hold(ax,'on');

    % Fetch year-i data
    T   = data_avg{i}.datetime(:);
    baz = data_avg{i}.median_back_azimuth(:);
    sd  = data_avg{i}.back_azimuth_std_dev(:);
    std_dev(i) = mean(sd);

    % Keep valid + time-sort
    valid = isfinite(baz) & isfinite(sd) & isfinite(datenum(T));
    T = T(valid); baz = baz(valid); sd = sd(valid);
    [T,ord] = sort(T); baz = baz(ord); sd = sd(ord);

    % Error bars (light gray, no caps)
    if ~isempty(T)
        errorbar(ax, T, baz, sd, ...
            'LineStyle','none','Color',[0.65 0.65 0.65], ...
            'CapSize',0,'LineWidth',lw);
    end

    % Points (black, slightly opaque to reduce clutter)
    if ~isempty(T)
        scatter(ax, T, baz, ms, 'k', 'filled', ...
            'MarkerFaceAlpha',0.9, 'MarkerEdgeColor','none');
    end

    % Axes formatting
    ylim(ax, yl);
    set(ax, 'YDir','reverse');          % same intent as axis ij
    xlim(ax, [datetime(years(i),5,1) datetime(years(i),8,31)]);

    % Month ticks & labels (May–Aug)
    xt = datetime(years(i), mo_ticks, 1);
    set(ax, 'XTick', xt, 'XTickLabel', cellstr(datestr(xt,'mmm')));

    % Grid & styling
    grid(ax,'on'); ax.MinorGridAlpha = 0.15; ax.GridAlpha = 0.25; ax.LineWidth = 0.75;
    title(ax, num2str(years(i)), 'FontSize', fs, 'FontWeight','bold');
    set(ax, 'FontSize', fs, 'FontName', '.AppleSystemUIFont');
    yticks(ax, 85:5:115);
end

% Shared labels
xlabel(tlo, 'Date', 'FontSize', fs+4, 'FontWeight','bold');
ylabel(tlo, 'Back azimuth (°)', 'FontSize', fs+4, 'FontWeight','bold');

disp(sprintf('Mean std dev across years: %.2f', nanmean(std_dev)));

% Optional export
% exportgraphics(gcf, 'back_azimuth_2015_2024.png', 'Resolution', 600);
% exportgraphics(gcf, 'back_azimuth_2015_2024.pdf');

%% Propagation speed (2015-2024)
close all;

fig2 = figure()
colormap('autumn')
for i = 1:length(years)
    subplot(5, 2, i)

    scatter(data_avg{i}.datetime, data_avg{i}.average_velocity, 'filled');
    hold on
    % Filter out NaN values for polynomial fit
    x = days_avg{i}; 
    y = data_avg{i}.average_velocity; 
    t = data_avg{i}.datetime; 
    valid_indices = ~isnan(y) & ~isinf(y); 
    x_valid = x(valid_indices);
    y_valid = y(valid_indices);
    t_valid = t(valid_indices); 

    % Perform polynomial fit
    p = polyfit(x_valid, y_valid, 2);

    % Generate fitted values over datetime range
    t_fit = linspace(min(t_valid), max(t_valid), 100); 
    x_fit = interp1(t_valid, x_valid, t_fit, 'linear', 'extrap'); 
    y_fit = polyval(p, x_fit); 

    % Plot the fit
    plot(t_fit, y_fit, 'k', 'LineWidth', 1.5);
    
    ylim([0 3]);
    title(num2str(years(i)), 'FontSize', fs, 'FontName', '.AppleSystemUIFont');
    grid on;
    time_x = [datetime(years(i), 5, 1, 3, 0, 0), datetime(years(i), 8, 31, 10, 0, 0)];
    xlim(time_x);

end

all=axes(fig2,'visible','off'); 
all.XLabel.Visible='on';
all.YLabel.Visible='on';

%% Boxplot of propagation speeds by month (all years combined)
close all;

all_dates = [];
all_vels  = [];

% Collect all years into one array
for i = 1:length(years)
    t = data_avg{i}.datetime;
    v = data_avg{i}.average_velocity;

    % Remove NaN/Inf values
    valid_idx = ~isnan(v) & ~isinf(v);
    all_dates = [all_dates; t(valid_idx)];
    all_vels  = [all_vels; v(valid_idx)];
end

% Extract months
all_months = month(all_dates);

% Boxplot
figure()
boxplot(all_vels, all_months, ...
    'Colors','k','Symbol','k+','Whisker',1.5);

ylim([0 3])
xlabel('Month')
ylabel('Propagation speed (m/s)')
title('Propagation speed by Month (2015–2024)')
grid on

%% Velocity vs Back Azimuth (day of year)

% close all;
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

% close all;
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
% close all;

figure()
plot_azimuth2(all_data_2020, 'Direction of wave source',  'May to September 2020')

figure()
plot_azimuth2(all_data_2019, 'Direction of wave source',  'May to September 2019')

%% Combined BA circle plots
close all;

figure()
p = polaraxes; 


for i = 1:length(years)
% for i = 3
    % Convert back azimuth to radians for polar plot
    azimuths_rad = deg2rad(data_all{i}.back_azimuth);
    
    % Scale normalized time to radius range
    radius = data_all{i}.velocity;
    
    % Create polar scatter plot with velocity as radius using scatter
    scatter(p, azimuths_rad, radius, 30, 'filled', 'DisplayName', num2str(years(i)), 'MarkerFaceAlpha', 1-(0.05*i));
    hold on
    
    % title('Direction of wave source',  'April-Sept 2019-2021', FontSize=14, FontAngle = 'italic');
    legend show;
    % hold off;
end

% Adjust the polar plot to set 0 degrees to North and clockwise direction
p.ThetaDir = 'clockwise';
p.ThetaZeroLocation = 'top';
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

function add_reference_lines()
    yline(3.23, 'b--', 'LineWidth',2);
    yline(2.22, 'b--', 'LineWidth',2);
    % text(-850, 3.33, 'mean speed in deep basin (Ramp et al., 2010)', 'color', 'b')
    % text(-850, 2.33, 'mean speed over cont. slope (Ramp et al., 2010)', 'color', 'b')
end
