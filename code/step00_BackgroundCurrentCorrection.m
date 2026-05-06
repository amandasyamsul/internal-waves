clear all; close all; clc
% Amanda Syamsul — April 30th 2025
% Correcting for background currents from HYCOM at lat 20.5 (all longitudes)
% Relevant outputs: Fig. S2

project_root = '/Users/amandasyamsul/Documents/MATLAB/OIW';

% Change folder path to choose current-depth at 25m OR 100m
data_folder = fullfile(project_root, '/hycom/currents/25m');
cd(data_folder);

addpath(fullfile(project_root, 'code'));

avg_data_2020 = readtable('avg_OIW_2020_data.csv');
all_data_2020 = readtable('OIW_2020_data.csv');
OBS_latband_data_2020 = readtable('IW_OBS_latband_2020.csv');

% The filenames for 25m are both uv100m. This is a just a naming error.
nc_files = dir(fullfile(data_folder, 'HYCOM_uv100m_*.nc'));
num_files = length(nc_files);

all_u = [];
all_v = [];
all_dates = []; 

for i = 1:num_files
    filename = nc_files(i).name;
    fpath = fullfile(data_folder, filename);
    disp(['Processing: ', filename]);

    % Check file size
    file_info = dir(fpath);
    if file_info.bytes < 5000
        warning(['Skipping corrupt or empty file: ', filename]);
        continue;
    end

    u = ncread(filename, 'water_u'); % eastward velocity
    v = ncread(filename, 'water_v'); % northward velocity
    lon = ncread(filename, 'lon');

    all_u = [all_u; u'];
    all_v = [all_v; v']; 
    all_dates = [all_dates; datetime(extractBetween(filename, 'HYCOM_uv100m_', '.nc'), 'InputFormat', 'yyyyMMdd')];

end

back_az_dates = avg_data_2020.datetime;
back_az = avg_data_2020.average_back_azimuth;

directional_currents = NaN(size(back_az));

for i = 1:height(avg_data_2020)
    date_idx = find(all_dates == back_az_dates(i));

    if isempty(date_idx)
        continue 
    end

    u_val = all_u(date_idx, :);  
    v_val = all_v(date_idx, :);

    theta_rad = deg2rad(back_az(i) + 180); % to convert to azimuth (propagation direction)
    all_dir_curr = u_val .* sin(theta_rad) + v_val .* cos(theta_rad);
    directional_currents(i) = mean(all_dir_curr);
end

%% Visualizing currents
figure()

ax(1) = subplot(3,1,1)

plot(all_dates, mean(all_u, 2), 'b-', 'LineWidth', 2)
hold on
plot(all_dates, mean(all_v, 2), 'r-', 'LineWidth', 2)
grid on
% ylim([-0.15 0.15])
ylabel('Background velocity (m/s)')
legend('Eastward velocity (u)', 'Northward velocity (v)')
title('Currents at 25m, lat 20.5, averaged across longitudes of satellite image')

ax(2) = subplot(3,1,2)
plot(back_az_dates, directional_currents, 'o-', 'LineWidth', 2)
grid on
ylabel('Background velocities in wave direction (m/s)')

ax(3) = subplot(3,1,3)
plot(back_az_dates, avg_data_2020.average_velocity, 'o-', 'LineWidth', 2)
hold on
plot(back_az_dates, avg_data_2020.average_velocity- directional_currents, 'o-', 'LineWidth', 2)

legend('Speeds from satellite imagery', 'Speeds adjusted for background velocities (m/s)')
ylabel('Propagation speeds (m/s)')

grid on;

linkaxes(ax, 'x')
%% Fig. S2

figure()

plot(back_az_dates, avg_data_2020.average_velocity, 'o-', 'LineWidth', 2)
hold on
plot(back_az_dates, avg_data_2020.average_velocity-directional_currents, 'g^', 'LineWidth', 2, 'MarkerFaceColor','g')
grid on
legend('Propagation speeds', 'Propagation speeds adjusted for background currents at 25m')
ylabel('Propagation speeds (ms^{-1})')

