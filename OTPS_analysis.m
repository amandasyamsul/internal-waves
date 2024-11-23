% Amanda Syamsul
% Nov 21st 2024
% Analysis on correlation between back azimuth and phase

close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS');

%% OTPS data
z = readtable('z_2015_2024.csv');
u = readtable('u_2015_2024.csv');

%% IW data
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

data_all = {all_data_2015, all_data_2016, all_data_2017, all_data_2018, all_data_2019, all_data_2020, all_data_2021, all_data_2022, all_data_2023, all_data_2024};
data_avg = {avg_data_2015, avg_data_2016, avg_data_2017, avg_data_2018, avg_data_2019, avg_data_2020, avg_data_2021, avg_data_2022, avg_data_2023, avg_data_2024};


% [l_time, l_vel, d_time, d_height] = tide_2019_2020();
% 
% l_dates = datetime(l_time, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% d_dates = datetime(d_time, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');

% l_vel = l_vel ./ 100; % Convert cm/s to m/s
% luzon = table(l_dates, l_vel);
% dongsha = table(d_dates, d_height);
%
% writetable(luzon,'luzon.csv')
% writetable(dongsha, 'dongsha.csv')

%% File processing -- only needs to be run once

% % Step 1: Read the .out file
% filename = 'z_2015_2024.out';  % Input .out file
% output_csv = 'z_2015_2024.csv'; % Desired output .csv file
% 
% % Open the file
% fid = fopen(filename, 'r');
% 
% % Step 2: Skip header lines
% header_lines = 6; % Adjust based on your file's structure
% for i = 1:header_lines
%     fgetl(fid); % Skip each header line
% end
% 
% % Step 3: Read data from the rest of the file
% data = textscan(fid, '%f %f %s %s %f %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
% 
% % Close the file
% fclose(fid);
% 
% % Step 4: Combine date and time into a single datetime column
% datetime_col = strcat(data{3}, {' '}, data{4}); % Combine date and time
% datetime_col = datetime(datetime_col, 'InputFormat', 'MM.dd.yyyy HH:mm:ss');
% 
% % Step 5: Create a table for output
% T = table(data{1}, data{2}, datetime_col, data{5}, data{6}, ...
%     'VariableNames', {'Lat', 'Lon', 'Datetime', 'Z_m', 'Depth_m'});
% 
% % Step 6: Write to a .csv file
% writetable(T, output_csv);
% 
% disp(['Data converted and saved as ', output_csv]);
% %%
% 
% % Step 1: Read the .out file
% filename = 'u_2015_2024.out';  % Input .out file
% output_csv = 'u_2015_2024.csv'; % Desired output .csv file
% 
% % Open the file
% fid = fopen(filename, 'r');
% 
% % Step 2: Skip header lines
% header_lines = 6; % Adjust based on your file's structure
% for i = 1:header_lines
%     fgetl(fid); % Skip each header line
% end
% 
% % Step 3: Read data from the rest of the file
% data = textscan(fid, '%f %f %s %s %f %f %f %f %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
% 
% % Close the file
% fclose(fid);
% 
% % Step 4: Combine date and time into a single datetime column
% datetime_col = strcat(data{3}, {' '}, data{4}); % Combine date and time
% datetime_col = datetime(datetime_col, 'InputFormat', 'MM.dd.yyyy HH:mm:ss');
% 
% % Step 5: Create a table for output
% T = table(data{1}, data{2}, datetime_col, data{5}, data{6}, data{7}, data{8}, data{9}, ...
%     'VariableNames', {'Lat', 'Lon', 'Datetime', 'U_m2_per_s', 'V_m2_per_s', 'u_cm_per_s', 'v_cm_per_s', 'Depth_m'});
% 
% % Step 6: Write to a .csv file
% writetable(T, output_csv);
% 
% disp(['Data converted and saved as ', output_csv]);

%% IW velocities corrected with E-W tidal velocities

years = 2015:2024;

for i = 1:length(years)

    time_x = [datetime(years(i),5,10,0,0,0), datetime(years(i),8,31,0,0,0)];
    
    [is_match, idx_in_u, idx_in_data] = intersect(u.Datetime, data_all{i}.datetime);
    
    % Create filtered tables for both datasets using the matched indices
    filtered_table = u(idx_in_u, :);
    filtered_table2 = data_all{i}(idx_in_data, :);
    
    normd_v = filtered_table2.velocity - (filtered_table.u_cm_per_s / 10e5);
    
    figure(i);
    ax1 = subplot(2, 1, 1)
    plot(u.Datetime, we_velocity./10e5, 'LineWidth', 1.5);
    ylabel('E-W Velocity (m/s)');
    title('E-W Velocity Over Time');
    grid on;
    xlim([time_x])
    
    ax2 = subplot(2, 1, 2)
    
    scatter(data_all{i}.datetime, data_all{i}.velocity, 'filled')
    hold on
    scatter(filtered_table2.datetime, normd_v, 'ro')
    legend('IW velocities', 'IW velocities corrected by E-W tidal velocities')
    xlim([time_x])
    ylabel('velocity (m/s)')
    title('Internal Wave velocities')
    grid on
    
    linkaxes([ax1 ax2], 'x');
end

%% Luzon Strait tidal plots for 2020

close all; 
time_x = [datetime(2022, 5, 1, 3, 0, 0), datetime(2020, 8, 31, 10, 0, 0)];

close all
figure(1)
% Plot 2: Tidal velocity
subplot(2, 1, 1); % Second subplot in a 2x2 grid
plot(luzon.l_dates, luzon.l_vel)
xlim(time_x - hours(44.5));
title('Tidal velocities in the Luzon Strait in 2020 (44.5 hour lag)');
xlabel('time');
ylabel('velocity (m/s)');

% Plot 1: Internal wave velocity
subplot(2, 1, 2); % First subplot in a 2x2 grid
plot_velocity(all_data_2020,'velocity')
% plot_avg_velocity(may_10km_60c, '10 km E-W of OBS - waves 60% in bounds above Dongsha')
xlim(time_x);
title('Internal wave velocities in 2020');


close all;

time_x = [datetime(2020, 7, 1, 3, 0, 0), datetime(2020, 7, 31, 11, 0, 0)];
figure()
subplot(2,1,1)
plot(luzon.l_dates, luzon.l_vel)
xlim(time_x - hours(44.5));
xline([avg_data.datetime]- hours(44.5), 'r-')
title('Tidal velocities in the Luzon Strait', FontSize=14, FontAngle = 'italic');
ylabel('Velocity (m/s)')

subplot(2,1,2)
plot_avg_velocity(avg_data, 'Average velocity')
xlim(time_x);
xline([avg_data.datetime], 'r-')
title('Daily average internal wave velocities', FontSize=14, FontAngle = 'italic');


