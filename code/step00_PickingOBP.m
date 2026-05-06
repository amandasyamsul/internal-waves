% Amanda Syamsul
% July 19th 2024
% Picking period and max amplitude from OBP

close all; clear all; clc

project_root = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(project_root);

% Add project folders
addpath(fullfile(project_root, 'code'));
addpath(fullfile(project_root, 'data'));
addpath(fullfile(project_root, 'seismic'));

detection = readtable('detections_reviewed.csv');
all_detections = readtable('all_detections_new.csv');
hdh = load('bb01_dec_1Hz_all_HDH.mat');
hhz = load('bb01_dec_1Hz_all_HH1.mat');
hh2 = load('bb01_dec_1Hz_all_HH2.mat');
hh1 = load('bb01_dec_1Hz_all_HHZ.mat'); % switch HH1 and HHZ (they were labeled wrong)

month_times = detection.DetectionTime;

% Timeframe for 2020 M7.5 Kuril Earthquake 
kuril = [datetime(2020, 3, 25, 2, 56, 16) datetime(2020, 3, 25, 3, 7, 16)];

% OBS coordinates
obs_coords = [21.00116, 117.40267];

% Convert the times from datenum to datetime (HDH, HH1, HH2, and HHZ all use same times)
t = datetime(hdh.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');

% Use this block when running full analysis
hdh_mt = t; 
hdh_md = hdh.d; 
hh1_md = hh1.d;
hh2_md = hh2.d;
hhz_md = hhz.d;

hdh_md_psi = hdh_md/1676128.86; % correction factor to convert to psi
hdh_md = hdh_md_psi * 6894.76; % convert psi to Pa

%% Define the time ranges to exclude around the artifacts
artifact_times=[datetime('5-May-2020 06:20:00'),
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
                datetime('15-August-2020 13:40:00'),
                datetime('18-August-2020 06:30:00'),
                datetime('25-August-2020 06:30:00') ];

exclude_durations = minutes(20);
artifact_end = artifact_times + exclude_durations;

% July has a data gap from 7/2-7/13, resume from 14-July-2020 14:30:00
gap_start = datetime('1-July-2020 00:00:00');
gap_end = datetime('14-July-2020 14:30:00');

% Create a logical index to only include valid data points
hdh_valid_idx = true(size(hdh_mt));

% Loop through each calibration pulse (artifact) time and update the valid indices
for i = 1:length(artifact_times)
    hdh_valid_idx = hdh_valid_idx & ~(hdh_mt >= artifact_times(i) & hdh_mt <= artifact_end(i));
end

% Update valid indices to exclude the gap period
hdh_valid_idx = hdh_valid_idx & ~(hdh_mt >= gap_start & hdh_mt <= gap_end);

% Filter out the artifact time ranges from the data
t_valid_hdh = hdh_mt(hdh_valid_idx);
d_valid_hdh = hdh_md(hdh_valid_idx);
d_valid_hhz = hhz_md(hdh_valid_idx);
d_valid_hh1 = hh1_md(hdh_valid_idx);
d_valid_hh2 = hh2_md(hdh_valid_idx);

% Interpolate data
t_interp = (hdh_mt(1):seconds(1):hdh_mt(end));
d_interp = interp1(t_valid_hdh, d_valid_hdh, t_interp);

d_interp_hhz = interp1(t_valid_hdh, d_valid_hhz, t_interp);
d_interp_hh1 = interp1(t_valid_hdh, d_valid_hh1, t_interp);
d_interp_hh2 = interp1(t_valid_hdh, d_valid_hh2, t_interp);

%% Raw vs Filtered OBP data for specific days
% close all;
clc;

nn = 320;
con_window=sin(pi*(1:nn)./nn);
con_window=con_window./sum(con_window);

window = 1.5337;

smoothed_data = conv(d_interp, con_window, 'same');

kuril = [datetime(2020, 3, 25, 2, 56, 16) datetime(2020, 3, 25, 3, 7, 16)];
wave_idx = 70;
figure()
plot(t_interp, d_interp-mean(d_interp), 'Color', [0.8 0.8 0.8]);
hold on
plot(t_interp, smoothed_data-mean(smoothed_data), 'k', 'LineWidth',1.5);
grid on
legend('OBP', 'OBP Smoothed')

xlim([all_detections.DetectionTime(wave_idx)-hours(window / 4) all_detections.DetectionTime(wave_idx)+hours(window / 2)] )

% Example of timeframe that was detected but is not an internal wave
% xlim([datetime(2019, 12, 12)-hours(window/4) datetime(2019, 12, 12)+hours(window/2) ]) 
%% Bandpass Butterworth filter design
% clc
% Fc1 = 1 / 10000;
% Fc2 = 1 / 100;
% [b, a] = butter(2, [Fc1 Fc2] / fnyq, 'bandpass'); % 2nd order Butterworth filter
% hdh_filtered_demeaned = filtfilt(b, a, hdh_detrended);
% hhz_filtered_demeaned = filtfilt(b, a, hhz_detrended);
% hh1_filtered_demeaned = filtfilt(b, a, hh1_detrended);
% hh2_filtered_demeaned = filtfilt(b, a, hh2_detrended);

%% Calculating time taken to cross OBS

% month_templates = detection.TemplateHDH(month_idx);
month_templates = detection.TemplateHDH;

dataset = readtable('IW_OBS_latband_2020.csv');
non_zero_idx = dataset.velocity > 0;
filt_data = dataset(non_zero_idx,:);

t_cross_obs = 10000 / (mean(filt_data.velocity)); % 1000 m because 5 km E and W of OBS

% Assuming month_times is your table and contains datetime values in a column
window = (t_cross_obs / 3600); % convert sec to hours

% Loop through each time in month_times

clc;

for times = 1:height(month_times)
    % Compute the start and end of the velocity window
    velocity_window_start(times) = month_times(times) - hours(window / 4);
    velocity_window_end(times) = month_times(times) + hours(window / 2);
end

col_names = {'window_start', 'detection_time', 'window_end'};
detection_window = table(velocity_window_start', month_times, velocity_window_end', 'VariableNames', col_names);

nn = 320;
con_window=sin(pi*(1:nn)./nn);
con_window=con_window./sum(con_window);

demeaned_data = d_interp - mean(d_interp);

my_data = conv(demeaned_data,con_window,'same');
my_time = t_interp;

for i = 1:height(detection_window)
    % Get the current detection window
    window_start = detection_window.window_start(i);
    window_end = detection_window.window_end(i);
    
    % Find the indices of the time points within the current detection window
    wave_passes_idx = (my_time >= window_start) & (my_time <= window_end);
    t_wave = my_time(wave_passes_idx);

    % Extract the filtered and detrended data within the current detection window
    wave_passes = my_data(wave_passes_idx); % no bandpass filt or detrend/demean
    
    % Calculate period, max amp, and min amp
    if ~isempty(wave_passes)
        [crest, crest_idx] = max(wave_passes);
        before_peak = wave_passes(1:crest_idx);
        [trough, trough_idx] = min(before_peak);
        % [amp_max, amp_idx] = max(abs(wave_passes));
        
        % Sign of signal
        wave_sign = sign(wave_passes);
        
        % Find where signal changes sign (zero crossings)
        zero_crossings = find(diff(wave_sign) ~= 0);
        
        % % T1: last zero before trough
        % T1_candidates = zero_crossings(zero_crossings < trough_idx);
        % if ~isempty(T1_candidates)
        %     T1_idx = T1_candidates(end);
        %     T1 = t_wave(T1_idx);
        % else
        %     T1 = NaT;
        % end
        % 
        % % T2: first zero after crest
        % T2_candidates = zero_crossings(zero_crossings > crest_idx);
        % if ~isempty(T2_candidates)
        %     T2_idx = T2_candidates(1);
        %     T2 = t_wave(T2_idx);
        % else
        %     T2 = NaT;
        % end

    end

    % Store the results
    % amp_time(i) = t_wave(amp_idx);
    crests(i) = crest;
    crest_time(i) = t_wave(crest_idx);
    troughs(i) = trough;
    trough_time(i) = t_wave(trough_idx);
    % amps(i) = crest-tough;
    % T1s(i) = T1';
    % T2s(i) = T2';
end

% detection.amp_time = amp_time';
detection.crests = crests';
detection.crest_time = crest_time';
detection.troughs = troughs';
detection.trough_time = trough_time';
detection.template = month_templates;
% detection.T1s = T1s';
% detection.T2s = T2s';
detection.period = (crest_time-trough_time)';
detection.amps = (crests-troughs)'

%% Check each wave manually
clc; close all;
% for f = 1:height(detection)
   
for f = [3, 61, 252, 376]

    time_x = [velocity_window_start(f), velocity_window_end(f)];
    time_indices = (my_time >= time_x(1)) & (my_time <= time_x(2));

    % Extract and demean within the time range
    data_subset = my_data(time_indices);

    % Extract corresponding time values
    time_subset = my_time(time_indices);

    figure()
    plot(t_interp, demeaned_data, 'Color', [0.8 0.8 0.8]); % unfiltered data in background
    hold on
    plot(time_subset, data_subset, 'LineWidth',2, 'Color','k')

    % Plot crests and troughs
    scatter(detection.crest_time(f), detection.crests(f), 200, 'filled', 'b^')
    scatter(detection.trough_time(f), detection.troughs(f),200,  'filled', 'rv')

    yline(0, 'LineWidth',1);
    legend('OBP', 'OBP Smoothed', 'Crest', 'Trough')
    grid on
    xlim(time_x);
    f
    % detection.Type(f)
    grid on;
end

% Examples in figure 4d are f = 3, 61, 252, 376

%% Manually remove signals that may not be internal waves

% bad_indices2 = [41 46 54 62 63 64 65 67 75 78 86 87 97 99 135 164 178 183 197 209 259 260 284 330 401 420 443 444 455] ;
% filtered_detections = detection;
% filtered_detections(bad_indices2, :) = []; % Removes the rows
% writetable(filtered_detections, 'all_detections_new.csv');

% all_detections5.csv: amp measured as max(abs(wave)) 
% all_detections_new.csv: amp measured as crest minus trough