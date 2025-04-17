% Amanda Syamsul
% July 19th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');

% setup filename file
% !ls -1 0*svg > Filenames

detection = readtable('detections_reviewed.csv');
hdh = load('bb01_dec_1Hz_all_HDH.mat');
% hh1 = load('bb01_dec_1Hz_all_HH1.mat');
% hh2 = load('bb01_dec_1Hz_all_HH2.mat');
% hhz = load('bb01_dec_1Hz_all_HHZ.mat');

%%
m = [5:8]; % month

% Convert the times from datenum to datetime
hdh_times = datetime(hdh.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% hh1_times = datetime(hh1.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% hh2_times = datetime(hh2.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% hhz_times = datetime(hhz.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');

col_names = {'det', 'type', 'template'};
detection_times = table(detection.DetectionTime, detection.Type, detection.TemplateHDH, 'VariableNames', col_names);

month_ind = (year(detection_times.det) == 2020) & ismember(month(detection_times.det), m);
month_times = detection_times.det(month_ind);


month_ind = (year(detection_times.det) == 2020) & ismember(month(detection_times.det), m);
month_times = detection_times.det(month_ind);

% Filter OBS times

hdh_month_ind = (year(hdh_times) == 2020) & ismember(month(hdh_times), m);
hdh_month_times = hdh_times(hdh_month_ind);

% hh1_month_ind = (year(hh1_times) == 2020) & ismember(month(hh1_times), m);
% hh1_month_times = hh1_times(hh1_month_ind);
% 
% hh2_month_ind = (year(hh2_times) == 2020) & ismember(month(hh2_times), m);
% hh2_month_times = hh2_times(hh2_month_ind);
% 
% hhz_month_ind = (year(hhz_times) == 2020) & ismember(month(hhz_times), m);
% hhz_month_times = hhz_times(hhz_month_ind);

% Extract the corresponding data points from hdh, hh1, hh2, and hhz
hdh_mt = hdh_times(hdh_month_ind);
hdh_md = hdh.d(hdh_month_ind);

% hh1_mt = hh1_times(hh1_month_ind);
% hh1_md = hh1.d(hh1_month_ind);
% 
% hh2_mt = hh2_times(hh2_month_ind);
% hh2_md = hh2.d(hh2_month_ind);
% 
% hhz_mt = hhz_times(hhz_month_ind);
% hhz_md = hhz.d(hhz_month_ind);


%% Design the Butterworth bandpass filter

% Define the sampling frequency (Fs) and cutoff frequency
Fs = 1; % Sampling frequency (1 Hz in this case)
fnyq = Fs / 2;
Fc1 = 1 / 500;
Fc2 = 1 / 400;

% Bandpass Butterworth filter design
[b, a] = butter(2, [Fc1 Fc2] / fnyq, 'bandpass'); % 2nd order Butterworth filter

% Define the time ranges to exclude around the artifacts
artifact_times = [datetime('5-May-2020 06:20:00'),
    datetime('12-May-2020 06:20:00'), 
    datetime('19-May-2020 06:20:00'), 
    datetime('26-May-2020 06:20:00'),
    datetime('2-June-2020 06:20:00'), 
    datetime('9-June-2020 06:20:00'), 
    datetime('16-June-2020 06:20:00'),
    datetime('23-June-2020 06:20:00'),
    datetime('30-June-2020 06:20:00'),
    datetime('21-July-2020 06:30:00'),
    datetime('28-July-2020 06:30:00'),
    datetime('4-August-2020 06:30:00'),
    datetime('11-August-2020 06:30:00'),
    datetime('15-August-2020 13:40:00'), % needs bigger exclusion
    datetime('18-August-2020 06:30:00'),
    datetime('25-August-2020 06:30:00')
    ];

exclude_durations = minutes(20);

% July has a data gap from 7/2-7/13, resume from 14-July-2020 14:30:00

% Create a logical index to include valid data points
% valid_indices_hh1 = true(size(hh1_mt));
% valid_indices_hh2 = true(size(hh2_mt));
% valid_indices_hhz = true(size(hhz_mt));
valid_indices_hdh = true(size(hdh_mt));

% Loop through each artifact time and update the valid indices
for i = 1:length(artifact_times)
    exclude_start = artifact_times(i);
    exclude_end = artifact_times(i) + exclude_durations;
    
    % valid_indices_hh1 = valid_indices_hh1 & ~(hh1_mt >= exclude_start & hh1_mt <= exclude_end);
    % valid_indices_hh2 = valid_indices_hh2 & ~(hh2_mt >= exclude_start & hh2_mt <= exclude_end);
    % valid_indices_hhz = valid_indices_hhz & ~(hhz_mt >= exclude_start & hhz_mt <= exclude_end);
    valid_indices_hdh = valid_indices_hdh & ~(hdh_mt >= exclude_start & hdh_mt <= exclude_end);
end

% Define the gap period
gap_start = datetime('1-July-2020 00:00:00');
gap_end = datetime('14-July-2020 14:30:00');

% Update valid indices to exclude the gap period
% valid_indices_hh1 = valid_indices_hh1 & ~(hh1_mt >= gap_start & hh1_mt <= gap_end);
% valid_indices_hh2 = valid_indices_hh2 & ~(hh2_mt >= gap_start & hh2_mt <= gap_end);
% valid_indices_hhz = valid_indices_hhz & ~(hhz_mt >= gap_start & hhz_mt <= gap_end);
valid_indices_hdh = valid_indices_hdh & ~(hdh_mt >= gap_start & hdh_mt <= gap_end);

% Filter out the artifact time ranges from the data
% hh1_mt_filtered = hh1_mt(valid_indices_hh1);
% hh1_md_filtered = hh1_md(valid_indices_hh1);
% 
% hh2_mt_filtered = hh2_mt(valid_indices_hh2);
% hh2_md_filtered = hh2_md(valid_indices_hh2);
% 
% hhz_mt_filtered = hhz_mt(valid_indices_hhz);
% hhz_md_filtered = hhz_md(valid_indices_hhz);

hdh_mt_filtered = hdh_mt(valid_indices_hdh);
hdh_md_filtered = hdh_md(valid_indices_hdh);

% Apply the bandpass filter to the data
% hh1_filtered = filtfilt(b, a, hh1_md_filtered);
% hh2_filtered = filtfilt(b, a, hh2_md_filtered);
% hhz_filtered = filtfilt(b, a, hhz_md_filtered);

% Detrend and de-mean the filtered data
% hh1_filtered_detrended = detrend(hh1_filtered) - mean(hh1_filtered);
% hh2_filtered_detrended = detrend(hh2_filtered) - mean(hh2_filtered);
% hhz_filtered_detrended = detrend(hhz_filtered) - mean(hhz_filtered);

%% Plotting seismic data for specific waves

% close all;
% time_x = [datetime(2020,7,17,2,30,0), datetime(2020,7,17,3,30,0)];
% 
% figure();
% 
% % Original data plot
% subplot(2, 1, 1);
% hold on;
% plot(hh1_mt, hh1_md, 'DisplayName', 'HH1 Original');
% plot(hh2_mt, hh2_md, 'DisplayName', 'HH2 Original');
% plot(hhz_mt, hhz_md, 'DisplayName', 'HHZ Original');
% title('Original Data');
% xlabel('Time');
% ylabel('Amplitude');
% xlim(time_x);
% legend;
% grid on;
% 
% % Filtered data plot
% subplot(2, 1, 2);
% hold on;
% plot(hh1_mt_filtered, hh1_filtered_detrended, 'DisplayName', 'HH1 Filtered');
% plot(hh2_mt_filtered, hh2_filtered_detrended, 'DisplayName', 'HH2 Filtered');
% plot(hhz_mt_filtered, hhz_filtered_detrended, 'DisplayName', 'HHZ Filtered');
% title('Filtered Data');
% xlabel('Time');
% ylabel('Amplitude');
% xlim(time_x);
% legend;
% grid on;

%% Calculating time taken to cross OBS

month_templates = detection_times.template(month_ind);

dataset = readtable('QC_IW_OBS_latband_2020.csv');
non_zero_idx = dataset.velocity > 0;
filt_data = dataset(non_zero_idx,:);

t_cross_obs = 20000 / (mean(filt_data.velocity)); % 20000 m because 10 km E and W of OBS

% Assuming month_times is your table and contains datetime values in a column
window = (t_cross_obs / 3600); % convert sec to hours

%% Loop through each time in month_times
for times = 1:height(month_times)
    % Compute the start and end of the velocity window
    velocity_window_start(times) = month_times(times) - hours(window / 2);
    velocity_window_end(times) = month_times(times) + hours(window / 2);
end

col_names = {'window_start', 'detection_time', 'window_end'};
detection_window = table(velocity_window_start', month_times, velocity_window_end', 'VariableNames', col_names);

my_data = hdh_md_filtered;
my_time = hdh_mt_filtered;

for i = 1:height(detection_window)
    % Get the current detection window
    window_start = detection_window.window_start(i);
    window_end = detection_window.window_end(i);
    
    % Find the indices of the time points within the current detection window
    wave_passes_idx = (my_time >= window_start) & (my_time <= window_end);
    t_wave = my_time(wave_passes_idx);

    % Extract the filtered and detrended data within the current detection window
    wave_passes = my_data(wave_passes_idx); % no bandpass filt or detrend/demean
    wave_passes = wave_passes - mean(wave_passes); % demean data
    
    % Calculate the maximum amplitude in the filtered data within the current window
    if ~isempty(wave_passes)
        [amp_max, amp_ind] = max(abs(wave_passes));
        [crest, crest_ind] = max(wave_passes);
        [trough, trough_ind] = min(wave_passes);
    end

    % Store the results
    amps(i) = amp_max;
    amp_time(i) = t_wave(amp_ind);
    crests(i) = crest;
    crest_time(i) = t_wave(crest_ind);
    troughs(i) = trough;
    trough_time(i) = t_wave(trough_ind);
end

%% 
detection_window.amps = amps';
detection_window.amp_time = amp_time';
detection_window.crests = crests';
detection_window.crest_time = crest_time';
detection_window.troughs = troughs';
detection_window.trough_time = trough_time';
detection_window.template = month_templates;

%%

% for i = 1:height(detection_window)
for i = 107
    % close all;

    time_x = [velocity_window_start(i), velocity_window_end(i)];
    time_indices = (my_time >= time_x(1)) & (my_time <= time_x(2));

    % Extract and demean within the time range
    data_subset = my_data(time_indices);
    demeaned_data = data_subset - mean(data_subset);

    % Extract corresponding time values
    time_subset = my_time(time_indices);

    figure()
    plot(time_subset, demeaned_data)
    hold on

    % Plot crests and troughs
    scatter(detection_window.crest_time(i), detection_window.crests(i), 'filled', 'r')
    scatter(detection_window.trough_time(i), detection_window.troughs(i), 'filled', 'g')

    legend('HHZ demeaned', 'Absolute max', 'Absolute min')
    title('Filtered Data');
    xlabel('Time');
    ylabel('Amplitude');
    grid on
    xlim(time_x);
    grid on;
end
%%

bad_indices = [18, 21, 23, 25, 26,27, 29, 30, 38, 63, 69, 116,124,137,141,143,171,174]; % List of indices to remove
filtered_waves = detection_window;
filtered_waves(bad_indices, :) = []; % Removes the row(s)

%%

% close all;

figure()
plot(filtered_waves.detection_time, filtered_waves.amps)
hold on
scatter(filtered_waves.detection_time, filtered_waves.crests, 'filled')
scatter(filtered_waves.detection_time, -filtered_waves.troughs, 'filled')
legend('max abs amplitude', 'crest amplitude','absolute trough amplitude')
ylim([0 7e5])
grid on
