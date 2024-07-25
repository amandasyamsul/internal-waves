% Amanda Syamsul
% July 19th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');

% setup filename file
% !ls -1 0*svg > Filenames

%% Read in data
may_10km_60c = readtable('may_10km_60c.csv');
june_10km_60c = readtable('june_10km_60c.csv');
july_10km_60c = readtable('july_10km_60c.csv');
detection = readtable('detections_reviewed.csv');
hdh = load('bb01_dec_1Hz_all_HDH.mat');
hh1 = load('bb01_dec_1Hz_all_HH1.mat');
hh2 = load('bb01_dec_1Hz_all_HH2.mat');
hhz = load('bb01_dec_1Hz_all_HHZ.mat');

%%
m = 7; % month
dataset = july_10km_60c;
%% Convert the times from datenum to datetime
hh1_times = datetime(hh1.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
hh2_times = datetime(hh2.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
hhz_times = datetime(hhz.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');

detection_times = detection.DetectionTime;
% Filter detection times for May 2020
month_ind = (year(detection_times) == 2020) & (month(detection_times) == m);
month_times = detection_times(month_ind);

% Filter OBS times for May 2020
hh1_month_ind = (year(hh1_times) == 2020) & (month(hh1_times) == m);
hh1_month_times = hh1_times(hh1_month_ind);

hh2_month_ind = (year(hh2_times) == 2020) & (month(hh2_times) == m);
hh2_month_times = hh2_times(hh2_month_ind);

hhz_month_ind = (year(hhz_times) == 2020) & (month(hhz_times) == m);
hhz_month_times = hhz_times(hhz_month_ind);

% Extract the corresponding data points from hh1, hh2, and hhz
hh1_mt = hh1_times(hh1_month_ind);
hh1_md = hh1.d(hh1_month_ind);

hh2_mt = hh2_times(hh2_month_ind);
hh2_md = hh2.d(hh2_month_ind);

hhz_mt = hhz_times(hhz_month_ind);
hhz_md = hhz.d(hhz_month_ind);

%% Design the Butterworth bandpass filter

% Define the sampling frequency (Fs) and cutoff frequency
Fs = 1; % Sampling frequency (1 Hz in this case)
fnyq = Fs / 2;
Fc1 = 1 / 1000;
Fc2 = 1 / 400;

% Bandpass Butterworth filter design
[b, a] = butter(2, [Fc1 Fc2] / fnyq, 'bandpass'); % 2nd order Butterworth filter

% Define the time ranges to exclude around the artifacts in May
artifact_times = [datetime('5-May-2020 06:20:00'),
    datetime('12-May-2020 06:20:00'), 
    datetime('19-May-2020 06:20:00'), 
    datetime('26-May-2020 06:20:00')];
exclude_durations = minutes(20);

% Define the time ranges to exclude around the artifacts in June
artifact_times = [datetime('2-June-2020 06:20:00'), 
    datetime('9-June-2020 06:20:00'), 
    datetime('16-June-2020 06:20:00'),
    datetime('23-June-2020 06:20:00'),
    datetime('30-June-2020 06:20:00')];

%Define the time ranges to exclude around the artifacts in July
artifact_times = [datetime('21-July-2020 06:30:00'),
    datetime('28-July-2020 06:30:00')];

exclude_durations = minutes(20);

% Create a logical index to include valid data points
valid_indices_hh1 = true(size(hh1_mt));
valid_indices_hh2 = true(size(hh2_mt));
valid_indices_hhz = true(size(hhz_mt));

% Loop through each artifact time and update the valid indices
for i = 1:length(artifact_times)
    exclude_start = artifact_times(i);
    exclude_end = artifact_times(i) + exclude_durations;
    
    valid_indices_hh1 = valid_indices_hh1 & ~(hh1_mt >= exclude_start & hh1_mt <= exclude_end);
    valid_indices_hh2 = valid_indices_hh2 & ~(hh2_mt >= exclude_start & hh2_mt <= exclude_end);
    valid_indices_hhz = valid_indices_hhz & ~(hhz_mt >= exclude_start & hhz_mt <= exclude_end);
end

% Exclude half of July
exclude_start = artifact_times(i);
exclude_end = artifact_times(i) + exclude_durations;

valid_indices_hh1 = valid_indices_hh1 & ~(hh1_mt >= exclude_start & hh1_mt <= exclude_end);
valid_indices_hh2 = valid_indices_hh2 & ~(hh2_mt >= exclude_start & hh2_mt <= exclude_end);
valid_indices_hhz = valid_indices_hhz & ~(hhz_mt >= exclude_start & hhz_mt <= exclude_end);

% Filter out the artifact time ranges from the data
hh1_mt_filtered = hh1_mt(valid_indices_hh1);
hh1_md_filtered = hh1_md(valid_indices_hh1);

hh2_mt_filtered = hh2_mt(valid_indices_hh2);
hh2_md_filtered = hh2_md(valid_indices_hh2);

hhz_mt_filtered = hhz_mt(valid_indices_hhz);
hhz_md_filtered = hhz_md(valid_indices_hhz);

% Apply the bandpass filter to the data
hh1_filtered = filtfilt(b, a, hh1_md_filtered);
hh2_filtered = filtfilt(b, a, hh2_md_filtered);
hhz_filtered = filtfilt(b, a, hhz_md_filtered);

% Detrend and de-mean the filtered data
hh1_filtered_detrended = detrend(hh1_filtered) - mean(hh1_filtered);
hh2_filtered_detrended = detrend(hh2_filtered) - mean(hh2_filtered);
hhz_filtered_detrended = detrend(hhz_filtered) - mean(hhz_filtered);

%% Plotting seismic data for specific waves

close all;
row = 16;
% time_x = [top_table.ClosestDetectionTime(row) - hours(1), top_table.ClosestDetectionTime(row) + hours(1)];
time_x = [datetime(2020,m,15,0,0,0), datetime(2020,m,30,24,0,0)];

figure(4);

% Original data plot
subplot(2, 1, 1);
hold on;
plot(hh1_mt, hh1_md, 'DisplayName', 'HH1 Original');
plot(hh2_mt, hh2_md, 'DisplayName', 'HH2 Original');
plot(hhz_mt, hhz_md, 'DisplayName', 'HHZ Original');
title('Original Data');
xlabel('Time');
ylabel('Amplitude');
xlim(time_x);
legend;
grid on;

% Filtered data plot
subplot(2, 1, 2);
hold on;
plot(hh1_mt_filtered, hh1_filtered_detrended, 'DisplayName', 'HH1 Filtered');
plot(hh2_mt_filtered, hh2_filtered_detrended, 'DisplayName', 'HH2 Filtered');
plot(hhz_mt_filtered, hhz_filtered_detrended, 'DisplayName', 'HHZ Filtered');
title('Filtered Data');
xlabel('Time');
ylabel('Amplitude');
xlim(time_x);
legend;
grid on;

%%
non_zero_indices = dataset.velocity > 0;
non_zero_may = dataset(non_zero_indices,:);

% Assuming month_times is your table and contains datetime values in a column, e.g., 'Timestamps'
window = 2.78; % Define the window

% Loop through each time in month_times
for times = 1:height(month_times)
    % Compute the start and end of the velocity window
    velocity_window_start(times) = month_times(times) - hours(window / 2);
    velocity_window_end(times) = month_times(times) + hours(window / 2);
end

col_names = {'window_start', 'detection_time', 'window_end'};
detection_window = table(month_times, velocity_window_start', velocity_window_end', 'VariableNames', col_names);

%%

% Initialize an array to store the average velocities for each window
average_velocities = zeros(height(detection_window), 1);

% Loop through each row in detection_window
for i = 1:height(detection_window)
    % Get the current detection window
    window_start = detection_window.window_start(i);
    window_end = detection_window.window_end(i);
    
    % Find the indices of velocities that fall within the current detection window
    in_window_indices = (non_zero_may.time > window_start) & (non_zero_may.time < window_end);
    
    % Extract the velocities that are within the current detection window
    velocities_in_window = non_zero_may.velocity(in_window_indices);
    
    % Calculate the average velocity within the current detection window
    if ~isempty(velocities_in_window)
        average_velocities(i) = mean(velocities_in_window);
    else
        average_velocities(i) = NaN; % Handle case with no data points in the window
    end
end

% Add the average velocities as a new column to detection_window
detection_window.average_velocity = average_velocities;


%%
% Close all figures
close all;

% Preallocate array to store the maximum amplitudes
hhz_amps = zeros(height(detection_window), 1);

for i = 1:height(detection_window)
    % Get the current detection window
    window_start = detection_window.window_start(i);
    window_end = detection_window.window_end(i);
    
    % Find the indices of the time points within the current detection window
    v_ind = (hhz_mt_filtered >= window_start) & (hhz_mt_filtered <= window_end);
    
    % Extract the filtered and detrended data within the current detection window
    hhz_v = hhz_filtered_detrended(v_ind);
    
    % Calculate the maximum amplitude in the filtered data within the current window
    if ~isempty(hhz_v)
        hhz_amp = max(abs(hhz_v));
    else
        hhz_amp = NaN; % Handle case with no data points in the window
    end
    
    % Store the result
    hhz_amps(i) = hhz_amp;
end

% Add the hhz_amps as a new column to detection_window
detection_window.hhz_amps = hhz_amps;

%%
% Plot the data
figure;
scatter(detection_window.average_velocity, detection_window.hhz_amps, 50, 'r', 'filled', 'DisplayName', 'Average velocity for each detection');

xlabel('Velocity (m/s)', 'FontSize', 12);
ylabel('HHZ Amplitude', 'FontSize', 12);
title('Internal wave velocity in June 2020 vs Max. Abs. HHZ Amplitude');

legend('show');
grid on;
box on;
hold off;