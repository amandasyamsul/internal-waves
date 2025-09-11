% Amanda Syamsul
% June 17th 2025
% Correcting IW speeds for currents using OTPS

close all; clear all; clc

data_dir = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(data_dir);

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/data');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

%% Defining datasets
all_data_2020 = readtable('OIW_2020_data.csv');
OBS_latband_data_2020 = readtable('IW_OBS_latband_2020.csv');
OBS_data_2020 = readtable('IW_OBS_2020.csv'); % square around OBS

obs_coords = [21.00116 117.40267];
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
square_around_OBS = OBS_data_2020;
latband_across_OBS = OBS_latband_data_2020;

% Define geographical bounds (longitude and latitude)
lon1 = 116.2; lon2 = 117.7;  % Longitude range of the study area
lat1 = 20.0; lat2 = 21.5;    % Latitude range of the study area

% Define x-y coordinate system bounds (based on satellite image pixels)
x_min = 0; x_max = 500;      % X-coordinates in image space
y_min = -400; y_max = 0;     % Y-coordinates in image space

% Compute the range of coordinates
image_width = x_max - x_min;    % Total width of the image in x-coordinates
image_height = y_max - y_min;   % Total height of the image in y-coordinates

lon_obs = lon1 + (square_around_OBS.x_coord - x_min) / image_width * (lon2 - lon1);
lon_obs_latband = lon1 + (latband_across_OBS.x_coord - x_min) / image_width * (lon2 - lon1);
square_around_OBS.longitude = lon_obs;
latband_across_OBS.longitude = lon_obs_latband;

%% Adjust speeds for current velocities
clear adjusted_speed;

% dataset = latband_across_OBS;
% 
% % If only plotting specific times, define the time range
% start_date = datetime(2020, 5, 1); 
% end_date = datetime(2020, 8, 31); 
% range = (dataset.datetime >= start_date) & (dataset.datetime <= end_date);
% dataset = dataset(range,:);
% idx_with_ba = ~isnan(dataset.back_azimuth);
% dataset = dataset(idx_with_ba,:);
% 
% tic;
% wb = waitbar(0, 'Adjusting speeds...');
% 
% for i = 1:height(dataset)
%     time = dataset.datetime(i);
%     u(i) = tmd_predict(tide_model, lat, dataset.longitude(i), time, 'u'); % + is east
%     v(i) = tmd_predict(tide_model, lat, dataset.longitude(i), time, 'v'); % + is north
% 
%     theta_rad = deg2rad(dataset.back_azimuth(i) + 180); % +180 to get direction of propagation instead of generation
%     current_velocity(i) = u(i) * sin(theta_rad) + v(i) * cos(theta_rad); % in direction of wave
% 
%     % Adjust propagation speed
%     adjusted_speed(i) = dataset.velocity(i) - current_velocity(i);
% 
%     waitbar(i / height(dataset), wb);
% end
% 
% close(wb);
% toc;
% 
% dataset.u = u';
% dataset.v = v';
% dataset.current_velocity = current_velocity';
% dataset.adjusted_speed = adjusted_speed';

% writetable(dataset, 'OBS_latband_2020_adjusted.csv');

%%
close all;

adjusted_data = readtable('OBS_latband_2020_adjusted.csv');

figure()
ax1 = subplot(2,1,1)
scatter(adjusted_data.datetime, adjusted_data.velocity, 'filled')
grid on
ylabel('propagation speed (m/s)')

ax2 = subplot(2,1,2)
scatter(adjusted_data.datetime, adjusted_data.u, 'filled', 'r')
hold on
scatter(adjusted_data.datetime, adjusted_data.v, 'filled', 'b')
grid on
ylabel('zonal velocity (m/s)')

linkaxes([ax1 ax2], 'x')