% Amanda Syamsul
% April 11th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');

% setup filename file
% !ls -1 0*svg > Filenames

%% Read in CSV files

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

det = readtable('detections_reviewed.csv');

% may_away = readtable('may_away.csv');
% may_near = readtable('may_near.csv');
% may_10km_60 = readtable('may_10km_60.csv');
% may_5km_60 = readtable('may_5km_60.csv');
% may_10km_60c = readtable('may_10km_60c.csv');
% may_5km_60c = readtable('may_5km_60c.csv');


%% Data calculation

[velocity,velocity_std_dev,velocity_std_err,back_azimuth,t_curr_array]=calculate_velocity('all_2015_data', false, true, 10000, 0.9);

% Saving data to table & SVG 

% Define column names
col_names = {'datetime', 'velocity','velocity_std_dev', 'back_azimuth'};
T = table(t_curr_array', velocity', velocity_std_dev', back_azimuth', 'VariableNames', col_names);

writetable(T,'OIW_2015_data.csv')

% Make table of averages

% Convert the time column to date-only format
T.date = dateshift(T.datetime, 'start', 'day');

% Group the data by date
[G, dateGroups] = findgroups(T.date);

% Calculate daily statistics while ignoring NaN values
daily_avg_velocity = splitapply(@(x) nanmean(x), T.velocity, G);
daily_median_velocity = splitapply(@(x) nanmedian(x), T.velocity, G);
daily_stdev_velocity = splitapply(@(x) nanstd(x), T.velocity, G);
daily_avg_backazimuth = splitapply(@(x) nanmean(x), T.back_azimuth, G);
daily_median_backazimuth = splitapply(@(x) nanmedian(x), T.back_azimuth, G);
daily_stdev_backazimuth = splitapply(@(x) nanstd(x), T.back_azimuth, G);
daily_max_backazimuth = splitapply(@(x) nanmax(x), T.back_azimuth, G);
daily_min_backazimuth = splitapply(@(x) nanmin(x), T.back_azimuth, G);

% Create a new table with the aggregated values
averages = table(dateGroups, daily_avg_velocity, daily_median_velocity, daily_stdev_velocity, daily_avg_backazimuth, daily_median_backazimuth, daily_stdev_backazimuth, daily_max_backazimuth, daily_min_backazimuth, ...
    'VariableNames', {'datetime', 'average_velocity','median_velocity', 'velocity_std_dev', 'average_back_azimuth', 'median_back_azimuth','back_azimuth_std_dev', 'max_ba', 'min_ba'});

writetable(averages,'avg_OIW_2015_data.csv')


%%

avgdayofyear15=(datenum(avg_data_2015.datetime)-datenum(2015,1,1));
avgdayofyear16=(datenum(avg_data_2016.datetime)-datenum(2016,1,1));
avgdayofyear17=(datenum(avg_data_2017.datetime)-datenum(2017,1,1));
avgdayofyear18=(datenum(avg_data_2018.datetime)-datenum(2018,1,1));
avgdayofyear19=(datenum(avg_data_2019.datetime)-datenum(2019,1,1));
avgdayofyear20=(datenum(avg_data_2020.datetime)-datenum(2020,1,1));
avgdayofyear21=(datenum(avg_data_2021.datetime)-datenum(2021,1,1));
avgdayofyear22=(datenum(avg_data_2022.datetime)-datenum(2022,1,1));
avgdayofyear23=(datenum(avg_data_2023.datetime)-datenum(2023,1,1));
avgdayofyear24=(datenum(avg_data_2024.datetime)-datenum(2024,1,1));

dayofyear15=(datenum(all_data_2015.datetime)-datenum(2015,1,1));
dayofyear16=(datenum(all_data_2016.datetime)-datenum(2016,1,1));
dayofyear17=(datenum(all_data_2017.datetime)-datenum(2017,1,1));
dayofyear18=(datenum(all_data_2018.datetime)-datenum(2018,1,1));
dayofyear19=(datenum(all_data_2019.datetime)-datenum(2019,1,1));
dayofyear20=(datenum(all_data_2020.datetime)-datenum(2020,1,1));
dayofyear21=(datenum(all_data_2021.datetime)-datenum(2021,1,1));
dayofyear22=(datenum(all_data_2022.datetime)-datenum(2022,1,1));
dayofyear23=(datenum(all_data_2023.datetime)-datenum(2023,1,1));
dayofyear24=(datenum(all_data_2024.datetime)-datenum(2024,1,1));

start_time = 15; % 3 PM in 24-hour format
end_time = 21 + 50/60; % 9:50 PM in decimal hours

data_all = {all_data_2015, all_data_2016, all_data_2017, all_data_2018, all_data_2019, all_data_2020, all_data_2021, all_data_2022, all_data_2023, all_data_2024};
days_all = {dayofyear15, dayofyear16, dayofyear17, dayofyear18, dayofyear19, dayofyear20, dayofyear21, dayofyear22, dayofyear23, dayofyear24};

data_avg = {avg_data_2015, avg_data_2016, avg_data_2017, avg_data_2018, avg_data_2019, avg_data_2020, avg_data_2021, avg_data_2022, avg_data_2023, avg_data_2024};
days_avg = {avgdayofyear15, avgdayofyear16, avgdayofyear17, avgdayofyear18, avgdayofyear19, avgdayofyear20, avgdayofyear21, avgdayofyear22, avgdayofyear23, avgdayofyear24};

years = 2015:2024;

timeofday = cell(1, length(years));

for i = 1:length(years)
    % Calculate time of day in decimal hours
    time_of_day = hour(data_all{i}.datetime) + minute(data_all{i}.datetime) / 60;
    
    % Assign each year's data to a cell in the cell array
    timeofday{i} = time_of_day;
end


%% Back Azimuth (2015-2024)

close all;
figure()
colormap('autumn')

years = 2015:2024;

for i = 1:length(years)
    subplot(2, 5, i)
    errorbar( data_avg{i}.datetime,  data_avg{i}.median_back_azimuth, ...
    data_avg{i}.back_azimuth_std_dev, 'w.', 'CapSize', 0)
    hold on
    scatter( data_avg{i}.datetime,  data_avg{i}.median_back_azimuth, ...
    55, days_avg{i}, 'filled')
    ylim([80 120])
    xlabel('datetime');
    ylabel('back azimuth');
    axis ij
    title(num2str(years(i)), FontSize=14, FontAngle = 'italic');
    grid on

    % Set axes background to dark gray
    set(gca, 'Color', [0.1 0.1 0.1]);
    set(gca, 'GridColor', 'w');
    time_x = [datetime(years(i), 5, 1, 3, 0, 0), datetime(years(i), 8, 31, 10, 0, 0)];
    xlim([time_x])

    caxis([120 240])
     
end

sgtitle('Daily medians of back azimuth in South China Sea')
c = colorbar('Position', [0.92 0.11 0.02 0.815]);
c.Label.String = 'Day of Year';
c.Ticks = [120, 150, 180, 210, 240];  % Monthly ticks
c.TickLabels = {'May', 'June', 'July', 'Aug', 'Sept'};


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
    % % Set figure background to black
    % set(gcf, 'Color', 'k');
    % 
    % Set axes background to dark gray
    set(gca, 'Color', [0.1 0.1 0.1]);
    % 
    % Set grid and text colors to white for visibility
    % set(gca, 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w');
    set(gca, 'GridColor', 'w');
end

% sgtitle('Daily medians of internal wave velocity against back azimuth in South China Sea')


% Set the color axis to ensure full coverage of tick labels
caxis([120 240])

% Add a single colorbar at the far right
c = colorbar('Position', [0.92 0.11 0.02 0.815]); 
c.Label.String = 'Day of Year';
c.Ticks = [120, 150, 180, 210, 240];  % Monthly ticks
c.TickLabels = {'May', 'June', 'July', 'Aug', 'Sept'};
% c.Label.Color = 'w';  % Set colorbar label text to white
% c.Color = 'w';        % Set colorbar tick labels to white
% 


%% Velocity vs Back Azimuth (time of day)

close all;
figure()
colormap('hot')

for i = 1:length(years)
    subplot(3, 3, i)
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


%% Running average (Back azimuth over time)
close all; 

num_points = 7;
h = ones(1, num_points)/ num_points;

figure()

for i = 1:length(years)

    subplot(3, 3, i)

    x = data_avg{i}.median_back_azimuth;
    t = data_avg{i}.datetime;

    y = conv(h, x);
    y = y(1: length(t));
    % y = y( ( ceil(num_points/2) : length(y)- floor(num_points/2) ) );
    std = data_avg{i}.back_azimuth_std_dev;

    plot(t, x, 'b', 'LineWidth', 1)
    hold on
    errorbar(t, x, std, std, 'k.', 'CapSize', 0)
    plot(t, y, 'ro-', 'LineWidth', 1)
    xlim([t(num_points), t(end)])
    axis ij
    grid on
    legend('median', 'std. deviation', '7-pt running average')
    xlabel('time')
    ylabel('back azimuth')
    title(num2str(years(i)), FontSize=14, FontAngle = 'italic');

end

%% All data on one same axis (colored by day of year)

figure(6)

for i = 1:length(years)
    scatter(data_all{i}.velocity, data_all{i}.back_azimuth, 40, days_all{i} , 'Filled')
    hold on
end

c = colorbar('Position', [0.92 0.11 0.02 0.815]); 
c.Label.String = 'Day of Year';
c.Ticks = [120, 150, 180, 210, 240];  % Monthly ticks
c.TickLabels = {'May', 'June', 'July', 'Aug', 'Sept'};

% Uncomment if coloring by time of day
% c = colorbar('Position', [0.92 0.11 0.02 0.815]); 
% c.Label.String = 'time of day';
% c.Ticks = [4, 5, 6, 7, 8, 9];  % Monthly ticks
% c.TickLabels = {'4 pm', '5 pm', '6 pm', '7 pm', '8 pm', '9 pm'};

title('Oceanic Internal Wave Velocity vs Back Azimuth 2015-2023')
grid on
ylim([80 120])
xlabel('velocity');
ylabel('back azimuth');

%% Plot velocities with different constraints
close all;

% figure()
% hold on;
% plot_velocity(may, 'total distance')
% plot_velocity(may_away, 'waves > 42 km E away from atoll')
% plot_velocity(may_near, 'waves < 42 km E from atoll')
% xlim(time_x);
% title('Internal wave velocities in May 2020');

%% Direction of wave source plots

close all; 

figure(10)

time_x = [datetime(2019, 4, 15, 3, 0, 0), datetime(2019, 9, 15, 10, 0, 0)];
subplot(4, 2, 1)
plot_azimuth(all_data_2019, 'Direction of Wave Source', '2019')
hold on
plot(avg_data_2019.datetime, avg_data_2019.average_back_azimuth, 'LineWidth',2)
xlim(time_x);
ylim([80,120])

time_x = [datetime(2020, 4, 15, 3, 0, 0), datetime(2020, 9, 15, 10, 0, 0)];
subplot(4, 2, 2)
plot_azimuth(all_data_2020, 'Direction of Wave Source', '2020')
hold on
plot(avg_data_2020.datetime, avg_data_2020.average_back_azimuth, 'LineWidth',2)
xlim(time_x);
ylim([80,120])

time_x = [datetime(2021, 4, 15, 3, 0, 0), datetime(2021, 9, 15, 10, 0, 0)];
subplot(4, 2, 3)
plot_azimuth(all_data_2021, 'Direction of Wave Source', '2021')
hold on
plot(avg_data_2021.datetime, avg_data_2021.average_back_azimuth, 'LineWidth',2)
xlim(time_x);
ylim([80,120])

time_x = [datetime(2022, 4, 15, 3, 0, 0), datetime(2022, 9, 15, 10, 0, 0)];
subplot(4, 2, 4)
plot_azimuth(all_data_2022, 'Direction of Wave Source', '2022')
hold on
plot(avg_data_2022.datetime, avg_data_2022.average_back_azimuth, 'LineWidth',2)
xlim(time_x);
ylim([80,120])

time_x = [datetime(2023, 4, 15, 3, 0, 0), datetime(2023, 9, 15, 10, 0, 0)];
subplot(4, 2, 5)
plot_azimuth(all_data_2023, 'Direction of Wave Source', '2022')
hold on
% plot(avg_data_2022.datetime, avg_data_2022.average_back_azimuth, 'LineWidth',2)
xlim(time_x);
ylim([80,120])


%% Individual year BA circle plots
close all;

figure()
plot_azimuth2(all_data_2020, 'Direction of wave source',  'May to September 2020')

figure()
plot_azimuth2(all_data_2019, 'Direction of wave source',  'May to September 2019')

%% Combined BA circle plot
close all;

figure()
% Convert back azimuth to radians for polar plot
azimuths_rad1 = deg2rad(all_data_2019.back_azimuth);
azimuths_rad2 = deg2rad(all_data_2020.back_azimuth);
azimuths_rad3 = deg2rad(all_data_2021.back_azimuth);

% Scale normalized time to desired radius range
radius1 = all_data_2019.velocity;
radius2 = all_data_2020.velocity;
radius3 = all_data_2021.velocity;

p = polaraxes; % Create a polar axes context
hold on;

% Adjust the polar plot to set 0 degrees to North and clockwise direction
p.ThetaDir = 'clockwise';
p.ThetaZeroLocation = 'top';

% Create polar scatter plot with velocity as radius using scatter
scatter(p, azimuths_rad1, radius1, 30, 'filled', 'DisplayName', '2019');
hold on
scatter(p, azimuths_rad2, radius2, 30, 'filled', 'DisplayName', '2020');
scatter(p, azimuths_rad3, radius3, 30, 'filled', 'DisplayName', '2021');

% Add a custom legend to indicate that the radius represents velocity
% h = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
% h2 = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
% 
% title(h, 'Radius = Velocity (m/s)');
% title(h2, 'Radius = Velocity (m/s)');

% Enhance the plot
title('Direction of wave source',  'April-Sept 2019-2021', FontSize=14, FontAngle = 'italic');
legend show;
hold off;

%% Only run this if OBS data was downloaded
% time_x = [datetime(2020, 5, 1, 3, 0, 0), datetime(2020, 7, 31, 10, 0, 0)];
% time_x = [datetime(2020, 7, 15, 3, 0, 0), datetime(2020, 8, 15, 10, 0, 0)];
% 
% figure;
% subplot(2, 1, 1);
% plot_azimuth(all_data, 'total distance')
% xlim(time_x);
% grid on;
% subplot(2, 1, 2);
% plot(hh1_mt_filtered, hh1_filtered_detrended, 'DisplayName', 'HH1 Filtered');
% plot(hh2_mt_filtered, hh2_filtered_detrended, 'DisplayName', 'HH2 Filtered');
% plot(hhz_mt_filtered, hhz_filtered_detrended, 'DisplayName', 'HHZ Filtered');
% title('Filtered Data');
% xlabel('Time');
% ylabel('Amplitude');
% xlim(time_x);
% legend;
% grid on;

%% Daily velocity plots

% maydays = [5, 6, 7, 9, 10, 12, 26, 27, 29];
% junedays = [4 11 19 20 21 22 24 25 27 28];
% julydays = [6 11 12 13 27];

julydays = [13 27];

figure(4)

for p = 1:length(julydays)
    time_x = [datetime(2020, 7, julydays(p), 1, 0, 0), datetime(2020, 7, julydays(p), 11, 0, 0)];
    
    % Plot 1: Tidal velocity
    subplot(length(julydays)*2, 1, 2*p-1); % Odd index for tidal velocity
    plot(luzon.l_dates, luzon.l_vel)
    xlim(time_x - hours(44.5));
    % title(['Tidal velocities in the Luzon Strait in May 2020 (44.5 hour lag) - Day ', num2str(maydays(p))]);
    % xlabel('time');
    % ylabel('velocity (m/s)');
    
    % Plot 2: Internal wave velocity
    subplot(length(julydays)*2, 1, 2*p); % Even index for internal wave velocity
    plot_velocity(all_data, 'July 2020')
    xlim(time_x);
    % title(['Internal wave velocities in May 2020 - Day ', num2str(maydays(p))]);
end

%%


% % Adding type to 2020 data
% 
% all_data_2020.type = repmat("Unknown", height(all_data_2020), 1);
% 
% % Loop through each row in all_data_2020
% for i = 1:height(all_data_2020)
%     time_window = abs(det.DetectionTime - all_data_2020.datetime(i)) <= hours(2);
% 
%     % If a match is found, use the first matched 'Type'
%     if any(time_window)
%         all_data_2020.type(i) = det.Type(find(time_window, 1));
%     end
% end
% 
% dayofyear20=(datenum(all_data_2020.datetime)-datenum(2020,1,1));
% 
% figure()
% hold on
% 
% % Separate data by type
% typeA = all_data_2020.type == "A";
% typeA2 = all_data_2020.type == "A2";
% typeUnknown = all_data_2020.type == "Unknown"; 
% 
% % Plot each type with a different marker shape
% scatter(all_data_2020.velocity(typeA), all_data_2020.back_azimuth(typeA), 50, dayofyear20(typeA), 'o', 'filled');
% scatter(all_data_2020.velocity(typeA2), all_data_2020.back_azimuth(typeA2), 70, dayofyear20(typeA2), 's');
% % scatter(all_data_2020.velocity(typeUnknown), all_data_2020.back_azimuth(typeUnknown), 30, dayofyear20(typeUnknown), '^', 'filled');
% 
% % Set axis limits, labels, and colorbar
% ylim([80 120])
% colorbar()
% xlabel('velocity');
% ylabel('back azimuth');
% axis ij
% title('Velocity & direction of wave source',  'May-Sept 2020', FontSize=14, FontAngle = 'italic');
% legend('A', 'A2', 'Unknown')
% 
% c = colorbar();
% c.Label.String = 'Day of Year';
% c.Ticks = [120, 150, 180, 210, 240]; % Roughly monthly ticks
% c.TickLabels = {'May', 'June', 'July', 'Aug', 'Sept'};
% 
% hold off

%%
function plot_azimuth2(data,label, subtitle)
    % Convert back azimuth to radians for polar plot
    azimuths_rad = deg2rad(data.back_azimuth);
    
    % Scale normalized time to desired radius range, e.g., [0, 2]
    radius = data.velocity;
    
    p = polaraxes; % Create a polar axes context
    hold on;

    % Adjust the polar plot to set 0 degrees to North and clockwise direction
    p.ThetaDir = 'clockwise';
    p.ThetaZeroLocation = 'top';

    % Create polar scatter plot with velocity as radius using scatter
    scatter(p, azimuths_rad, radius, 30, 'filled');
    
    % Add a custom legend to indicate that the radius represents velocity
    h = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
    title(h, 'Radius = Velocity (m/s)');

    % Enhance the plot
    title(label, subtitle, FontSize=14, FontAngle = 'italic');
    legend show;
    hold off;
end


function plot_azimuth(data, label, subtitle)
   
    % Plot mean back azimuth with error bars
    % errorbar(datetime, back_azimuth, std_deback_azimuth, 'r-', 'LineWidth', 1, 'DisplayName', label);
    
    % Plot individual back azimuth values
    scatter(data.datetime, data.back_azimuth, 20, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
    
    grid on;
    
    % Label the axes
    title(label, subtitle, FontSize=14, FontAngle = 'italic');
    
    xlabel('Time');
    ylabel('Back Azimuth (degrees)');
    
    % Show legend
    % legend show;
    
    hold off;
end

function plot_avg_velocity(data, caption)

    % Plot velocity vs time with error bars
    % errorbar(data.datetime, data.average_velocity, data.velocity_std_dev, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'standard deviation');
    % scatter(data.datetime, data.average_velocity, 'ro', 'Filled')
    hold on
    plot(data.datetime, data.average_velocity,'-o','LineWidth', 1,'DisplayName', caption);
    xlabel('Time');
    ylabel('Velocity (m/s)');
    grid on;
    set(gca, 'FontSize', 12);
end

function plot_velocity(data, caption)
    % Plot velocity vs time with error bars
    errorbar(data.datetime, data.velocity, data.velocity_std_err, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'standard error');
    hold on
    plot(data.datetime, data.velocity,'-','LineWidth', 1,'DisplayName', caption);
    xlabel('Time');
    ylabel('Velocity (m/s)');
    legend show;
    grid on;
    set(gca, 'FontSize', 12);
end

function [velocity,velocity_std_dev,velocity_std_err,back_azimuth,t_curr_array]=calculate_velocity(data_file, verbose, only_obs, constraint, prop)

    data = readtable(data_file, 'ReadVariableNames', false); 
    
    data_array = table2array(data);

    for f=1:length(data_array)-1
        file1 = data_array{f};
        file2 = data_array{f+1};
    
        try
            prev = loadsvg(file1,1,1);
        catch
            continue; % Skip to the next iteration if file1 causes an error
        end
        
        try
            curr = loadsvg(file2,1,1);
        catch   
            continue; % Skip to the next iteration if file2 causes an error
        end

        obs = loadsvg('OBS.svg',1,1);
    
        % Failsafe in case there is a smaller vector in the svg that should
        % be ignored
        if (length(prev)>1)
            for icheck=1:length(prev)
                Ncheck(icheck)=length(prev{icheck});
            end
            [~,Ilong]=max(Ncheck);
            prev=prev(Ilong); % pick the longest vector
        end
    
        if (length(curr)>1)
            for icheck=1:length(curr)
                Ncheck(icheck)=length(curr{icheck});
            end
            [~,Ilong]=max(Ncheck);
            curr=curr(Ilong); % pick the longest vector
        end
    
        close all; 
        
        dt_steps = readmatrix('dt-picking.csv');
        
        % Extract last three digits from file names
        prev_num = str2double(file1(end-14:end-11));
        curr_num = str2double(file2(end-14:end-11));

        % Find the corresponding dt values from the CSV data
        prev_dt = find(dt_steps(:,1) == prev_num);
        curr_dt = find(dt_steps(:,1) == curr_num);
        dt = (curr_dt - prev_dt) * 600;

        % Check if dt is less than 0 (between two diff. days), and if so, skip to the next iteration
        if dt < 0
            continue;
        end
        
        if verbose
          disp(['Time difference (dt) between SVG files: ' num2str(dt) ' seconds']);
        end
      
        % Set datetime format for previous wave
        % Extract month, day, and time parts from the filename
        year_prev = str2double(file1(1:4));
        month_prev = str2double(file1(5:6));
        day_prev = str2double(file1(7:8));
        time_prev = file1(9:12);
        
        % Determine hour and minute from the time part
        if length(time_prev) == 4
            hour_prev = str2double(time_prev(1:2));  % First two digits are the hour
            min_prev = str2double(time_prev(3:4)); % Last two digits are the minutes
        else
            hour_prev = str2double(time_prev(1:2)); % First two digits are the hour
            min_prev = str2double(time_prev(3:4)); % Last two digits are the minutes
        end
        
        % Construct the datetime object assuming time is in AM UTC
        t_prev = datetime(year_prev, month_prev, day_prev, hour_prev, min_prev, 0, 'TimeZone', 'UTC');

        % Set datetime format for previous wave
        % Extract month, day, and time parts from the filename
        year_curr = str2double(file2(1:4));
        month_curr = str2double(file2(5:6));
        day_curr = str2double(file2(7:8));
        time_curr = file2(9:12);
        
        % Determine hour and minute from the time part
        if length(time_curr) == 4
            hour_curr = str2double(time_curr(1:2));  % First two digits are the hour
            min_curr = str2double(time_curr(3:4)); % Last two digits are the minutes
        else
            hour_curr = str2double(time_curr(1:2)); % First two digits are the hour
            min_curr = str2double(time_curr(3:4)); % Last two digits are the minutes
        end
        
        % Construct the datetime object assuming time is in AM UTC
        t_curr = datetime(year_curr, month_curr, day_curr, hour_curr, min_curr, 0, 'TimeZone', 'UTC');        
        
        % Define latitude range covered by the image in meters
        % Using scale bar in satellite images, 25 km = 70.03 pixels
        meters_per_pix = 25000/70.03;
        
        % close all; clc
        
        original_x_prev=prev{1,1}(:,1);
        original_y_prev=prev{1,1}(:,2);
        original_x_curr=curr{1,1}(:,1);
        original_y_curr=curr{1,1}(:,2);


        %% Interpolate data
    
        % Previous wave:
        % Combine x and y into a matrix and sort rows by the first column (x)
        data = [original_x_prev, original_y_prev];
        data = sortrows(data, 1);
    
        % Find unique x values and average corresponding y values
        [y_unique, ia, ic] = unique(data(:,2));
        x_avg = accumarray(ic, data(:,1), [], @mean);
    
        % Interpolation query points
        yq = linspace(min(y_unique), max(y_unique), numel(y_unique)*10);  % A more reliable query range
    
        % Interpolate
        vq1 = interp1(y_unique, x_avg, yq);
    
        % Combine and sort x values
        x_combined = [data(:,1); vq1'];
        y_combined = [data(:,2); yq'];
    
        % Sort the combined arrays by x values
        [y_sorted, sort_index] = sort(y_combined);
        x_sorted = x_combined(sort_index);
    
        x_prev = x_sorted;
        y_prev = y_sorted;
    
        % Curent wave:
    
        % Combine x and y into a matrix and sort rows by the first column (x)
        data = [original_x_curr, original_y_curr];
        data = sortrows(data, 1);
    
        % Find unique x values and average corresponding y values
        [y_unique, ia, ic] = unique(data(:,2));
        x_avg = accumarray(ic, data(:,1), [], @mean);
    
        % Interpolation query points
        yq = linspace(min(y_unique), max(y_unique), numel(y_unique)*10);  % A more reliable query range
    
        % Interpolate
        vq1 = interp1(y_unique, x_avg, yq);
    
        % Combine and sort x values
        x_combined = [data(:,1); vq1'];
        y_combined = [data(:,2); yq'];
    
        % Sort the combined arrays by x values
        [y_sorted, sort_index] = sort(y_combined);
        x_sorted = x_combined(sort_index);
    
        x_curr = x_sorted;
        y_curr = y_sorted;

        % Cropping the longer vector to match wave lengths
        
        % Calculate y-lengths (ranges)
        yLengthPrev = max(y_prev) - min(y_prev);
        yLengthCurr = max(y_curr) - min(y_curr);
        
        % Compare y-lengths and choose the vector with the shorter y
        if yLengthCurr > yLengthPrev
            short_x = x_prev;
            short_y = y_prev;
            chosen_name = 'Previous';
            long_x = x_curr;
            long_y = y_curr;
        else
            short_x = x_curr;
            short_y = y_curr;
            chosen_name = 'Previous';
            long_x = x_prev;
            long_y = y_prev;
        end
        
        % Find the minimum and maximum y-values of the shorter vector
        ymin = min(short_y);
        ymax = max(short_y);
        
        % Find the indices in the longer vector where the y-values are within the range of the shorter vector
        indices = long_y >= ymin & long_y <= ymax;
        
        % Crop the longer vector to match the range and y-max/y-min of the shorter vector
        cropped_long_y = long_y(indices);
        cropped_long_x = long_x(indices);

        % Compare y-lengths and replace longer wave vector with cropped version
        if yLengthCurr > yLengthPrev
            y_curr = cropped_long_y;
            x_curr = cropped_long_x;
        else
            y_prev = cropped_long_y;
            x_prev = cropped_long_x;
        end
                
        if verbose
            % Visualization for verification
            figure;
            plot(x_prev, y_prev, 'b-', 'DisplayName', 'Previous Wave');
            hold on;
            plot(x_curr, y_curr, 'g-', 'DisplayName', 'Current Wave');
            plot(cropped_long_x, cropped_long_y, 'r--', 'DisplayName', sprintf('Adjusted %s Wave', chosen_name));
            plot(392,159, 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);
            legend show;
            xlabel('X Axis');
            ylabel('Y Axis');
            title(sprintf('Comparison of Adjusted and Original Vectors'));
        end


        %% Use to only calculate near OBS

        if only_obs
            %% Define the boundary based on constraint and meters_per_pix
            bound = constraint / meters_per_pix;

            %% Check the proportion of x_prev and x_curr within the bounds of OBS
            % in_bounds_prev = (x_prev >= (392-bound)) & (x_prev <= (392+bound));
            % in_bounds_curr = (x_curr >= (392-bound)) & (x_curr <= (392+bound));

            % Cut off waves at atoll
            in_bounds_prev = (x_prev >= 258);
            in_bounds_curr = (x_curr >= 258);

            % Calculate the proportion of points within the bounds
            proportion_in_bounds_prev = sum(in_bounds_prev) / length(x_prev);
            proportion_in_bounds_curr = sum(in_bounds_curr) / length(x_curr);

            % Check if at least prop% of the points are within the bounds for both waves
            if proportion_in_bounds_prev < prop || proportion_in_bounds_curr < prop
                continue; % Skip iteration if less than prop% of the points are within the bounds
            end

            % COMMENT OUT code below if USING entire y-axis
            % % Filter the arrays to keep only values above Dongsha (y=251)
            % obs_filter_prev = y_prev <= 251; % less than because points are flipped
            % obs_filter_curr = y_curr <= 251;
            % 
            % % Apply the filters
            % x_prev = x_prev(obs_filter_prev);
            % y_prev = y_prev(obs_filter_prev);
            % x_curr = x_curr(obs_filter_curr);
            % y_curr = y_curr(obs_filter_curr);

            % Skip iteration if any of the arrays are empty after filtering
            if isempty(x_prev) || isempty(y_prev) || isempty(x_curr) || isempty(y_curr)
                continue;
            end
            %

            if verbose
                figure;
                plot(x_prev, y_prev, '-*', 'DisplayName', 'Previous Wave');
                hold on;
                plot(x_curr, y_curr, '-*', 'DisplayName', 'Current Wave');
                plot(392, 159, 'p', 'DisplayName', 'OBS Station', 'MarkerFaceColor', 'red', 'MarkerSize', 15);
                plot(199, 251, 'p', 'DisplayName', 'Dongsha', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);
                legend show;
            end
        end

        % Failsafe for OBS -- skip if one wave is too short
        if length(x_prev) < 0.3*length(x_curr)
            continue
        end

        if length(x_curr) < 0.3*length(x_prev)
            continue
        end

        %% Finding measurement points:
    
        % Define the fractions of the line where you need the indices
        fractions = [1/10, 2/10, 3/10, 4/10, 1/2, 6/10, 7/10, 8/10, 9/10];
    
        % Preallocate the matrix to store results
        calc_points = zeros(length(fractions), 2);

        % Loop over each fraction
        for i = 1:length(fractions)
            % Calculate the indices for the current and previous lines
            index_prev = round(length(y_prev) * fractions(i));
            index_curr = round(length(y_curr) * fractions(i));

            
            % Test code to see if this solves issue of rounding to 0
            if index_prev==0
                continue
            end

            if index_curr==0
                continue
            end
    
            % Store the results in the matrix
            calc_points(i, :) = [index_prev index_curr];
        end
            
        for i=1:length(fractions)
            % Calculate the slope at varying points on the previous wave
            gap=10;
    
            % edit gap for waves that have too few points!
            if (calc_points(i,1) - gap)<1
                gap=round(calc_points(i,1)/2);
            end
    
            if (calc_points(i,1) - gap)==0
                gap=0;
            end
            
            % Failsafe if vector is still too short -- skip to next iteration
            if (calc_points(i,1) - gap)==0
                continue
            end

            x1a = x_prev(calc_points(i,1) - gap);
            y1a = y_prev(calc_points(i,1) - gap);
            
            gap = 10;
            if (calc_points(i,1) + gap) > length(x_prev)
                gap = 1; % this has to be fixed!!!!
                if (calc_points(i,1) + gap) > length(x_prev)
                    continue;
                end
            end

            x2a = x_prev(calc_points(i,1) + gap);
            y2a = y_prev(calc_points(i,1) + gap);
            slope_prev(i) = (y2a - y1a) / (x2a - x1a);
            
            % Calculate the slope at varying points on the current wave
            gap=10;
            if (calc_points(i,2) - gap)<1
                gap=round(calc_points(i,2)/2);
            end
    
            if (calc_points(i,2) - gap)==0
                gap=0;
            end
    
            x1a= x_curr(calc_points(i,2) - gap);
            y1a= y_curr(calc_points(i,2) - gap);
            
            gap=10;

            if (calc_points(i,2) + gap)>length(x_curr)
                gap=1; % this has to be fixed!!!!
            end
            
            % Failsafe for OBS
            if (calc_points(i,2) + gap)>length(x_curr)
                continue
            end

            x2a = x_curr(calc_points(i,2) + gap);
            y2a = y_curr(calc_points(i,2) + gap);
            slope_curr(i) = (y2a - y1a) / (x2a - x1a);
        
            % Measure perp. dist from slope mid point to closest point on wave of next image
        
            % Define the coordinates of the points
            
            m1(i) = slope_prev(i);
            m2(i) = slope_curr(i);
        
            x1(i) = x_prev(calc_points(i,1));
            y1(i) = y_prev(calc_points(i,1));
        
            x2(i) = x_curr(calc_points(i,2));
            y2(i) = y_curr(calc_points(i,2));
    
            xs_num(i) = y2(i) - (m2(i)*x2(i)) -(x1(i)/m1(i)) -y1(i);
            xs_den(i) = -(1/m1(i)) -m2(i);
            xs(i) = xs_num(i)/xs_den(i);
            ys(i) = m2(i) * (xs(i) - x2(i)) + y2(i);
            ms(i) = (ys(i) - y1(i)) / (xs(i) - x1(i));
            
            xvals = 200:500;
            y_perp(:,i) = -(1/m1(i)) * (xvals - x1(i)) + y1(i);
    
            point1(i,:) = [x1(i) y1(i)];
            point2(i,:) = [xs(i) ys(i)];
        
            % Calculate velocity using all points
            if any(isnan(point1(i,:))) || any(isnan(point2(i,:)))
                continue
            end
                
            distance_pix(i) = norm(point1(i,:) - point2(i,:));
            distance_m(i) = distance_pix(i) * meters_per_pix;
            v(i) = distance_m(i)/dt;
        end
        
        try
            % m* (ms) is the slope of the line connecting both waves

            % Calculate the slope angle in degrees
            mean_ms = mean(ms);

            % Assuming mean_ms is a slope value, convert it to an angle in degrees
            ms_degree = atand(mean_ms);
            
            % Convert the slope angle to azimuth
            % Assume north is 0 degrees, east is 90 degrees, south is 180 degrees, west is 270 degrees
            if ms_degree >= 0
                back_az = 90 + ms_degree;
            else
                back_az = 90 - abs(ms_degree);
            end

            back_azimuth(f+1) = mean(back_az);
            velocity(f+1) = mean(v);
            velocity_std_dev(f+1)=std(v);

            % standard error
            n = length(v);
            mean_velocity = mean(v);
            velocity_std_err(f+1) = velocity_std_dev(f+1) / sqrt(n);

        catch
            continue
        end

        % Store the datetime object
        t_curr_array(f+1) = t_curr;

        if verbose
            disp(t_curr);
        end
    
        %% logic testing
        if verbose
            figure(10)
            hold on
            
            plot(x_prev,y_prev,'b'); leg7 = "previous wave";
            plot(x_curr,y_curr,'k'); leg8 = "current wave";
            
            for p = 1:length(fractions)

                if calc_points(p, 1) == 0
                    continue
                end
                
                % if length(fractions) > length(p)
                %     continue
                % end

                plot(x_prev(calc_points(p, 1)), y_prev(calc_points(p,1)), 'bo')
                plot(x_curr(calc_points(p, 2)), y_curr(calc_points(p,2)), 'ko')
                plot(xs(p),ys(p),'m*'); leg4 = "perp. intersect";
                plot(x1(p),y1(p),'b*'); leg4 = "";
                plot(x2(p),y2(p),'k*'); leg4 = "perp. intersect";
                y_perp(:,p) = -(1/m1(p)) * (xvals - x1(p)) + y1(p);
                plot(xvals,y_perp(:,p), 'r--'); leg3 = "perp. line";
                plot(392,159, 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);
                plot(199, 251, 'p', 'DisplayName', 'Dongsha', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);

            end
            
        end

        if verbose
            disp(['Average velocity between SVG files: ' num2str(velocity) 'm/s'])
            disp(['next!'])
        end
    end
end

function compare(manual_velocity,auto_velocity,auto_velocity_error)
    manual = readtable(manual_velocity);
    lim = height(manual);
    avg_manual = mean(manual.v(2:end));
    avg_auto = mean(auto_velocity(2:end));
    diff_avg = abs(avg_manual-avg_auto);

    figure(20)
    hold on
    plot(manual.seconds(1:lim),manual.v(1:lim), 'DisplayName', 'manual measurements');
    errorbar(manual.seconds(1:lim), auto_velocity(1:lim), auto_velocity_error(1:lim),'DisplayName', 'automatic measurements');
    xlabel('time (seconds)');
    ylabel('manually measured velocity (m/s)');
    title('comparison between manually and automatically measured velocities');
    disp(['diff. in avg. velocity = ', num2str(diff_avg), ' m/s'])
    legend show;
    hold off
end

function filtered_array = filter_continuous_ones(array)
    
    % Find the start and end indices of each contiguous segment of 1s
    d = diff([0; array; 0]);
    start_indices = find(d == 1);
    end_indices = find(d == -1) - 1;
    
    % Calculate the length of each segment
    segment_lengths = end_indices - start_indices + 1;
    
    % Find the index of the longest segment
    [~, longest_segment_index] = max(segment_lengths);
    
    % Get the start and end indices of the longest segment
    longest_segment_start = start_indices(longest_segment_index);
    longest_segment_end = end_indices(longest_segment_index);
    
    % Initialize the output array with zeros
    filtered_array = zeros(size(array));
    
    % Set the elements of the longest segment to 1
    filtered_array(longest_segment_start:longest_segment_end) = 1;

    filtered_array = logical(filtered_array);
end
