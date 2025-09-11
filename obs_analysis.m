% Amanda Syamsul
% July 19th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc

dr = '/Users/amandasyamsul/Documents/MATLAB/OIW/code';
cd(dr);

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/data/');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');

% setup filename file
% !ls -1 0*svg > Filenames

detection = readtable('detections_reviewed.csv');
all_detections = readtable('all_detections3.csv');
hdh = load('bb01_dec_1Hz_all_HDH.mat');
hhz = load('bb01_dec_1Hz_all_HHZ.mat');
hh1 = load('bb01_dec_1Hz_all_HH1.mat');
hh2 = load('bb01_dec_1Hz_all_HH2.mat');

m = [5:8]; % month

% OBS coordinates
obs_coords = [21.00116, 117.40267];

% Convert the times from datenum to datetime
hdh_times = datetime(hdh.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
hh1_times = datetime(hh1.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
hh2_times = datetime(hh2.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
hhz_times = datetime(hhz.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');

% col_names = {'det', 'type', 'template'};
% detection_times = table(detection.DetectionTime, detection.Type, detection.TemplateHDH, 'VariableNames', col_names);

% month_idx = (year(detection.DetectionTime) == 2020) & ismember(month(detection.DetectionTime), m);
% month_times = detection.DetectionTime(month_idx);
month_times = detection.DetectionTime;

% Filter data to match desired timeframe
hdh_month_idx = (year(hdh_times) == 2020) & ismember(month(hdh_times), m);
hdh_month_times = hdh_times(hdh_month_idx);

hh1_month_idx = (year(hh1_times) == 2020) & ismember(month(hh1_times), m);
hh1_month_times = hh1_times(hh1_month_idx);

hh2_month_idx = (year(hh2_times) == 2020) & ismember(month(hh2_times), m);
hh2_month_times = hh2_times(hh2_month_idx);

hhz_month_idx = (year(hhz_times) == 2020) & ismember(month(hhz_times), m);
hhz_month_times = hhz_times(hhz_month_idx);

% hdh_mt = hdh_times(hdh_month_idx);
hdh_mt = hdh_times;

% hdh_md = hdh.d(hdh_month_idx);
hdh_md = hdh.d;
hhz_md = hhz.d;
hh1_md = hh1.d;
hh2_md = hh2.d;

hdh_md = hdh_md/1676128.86; % correction factor to convert pressure to psi
%%
% close all;
% 
% figure()
% ax(1) = subplot(2,1,1)
% plot(hdh_mt, hdh_md, 'DisplayName', 'HDH');
% ax(2) = subplot(2,1,2)
% plot(hh1_mt, hh1_md, 'DisplayName', 'HH1');
% hold on
% plot(hhz_mt, hhz_md, 'DisplayName', 'HHZ');
% plot(hh2_mt, hh2_md, 'DisplayName', 'HH2');
% legend()
% linkaxes([ax], 'x')
%% Define the time ranges to exclude around the artifacts
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

% Create a logical index to only include valid data points
valid_indices_hdh = true(size(hdh_mt));

% Loop through each calibration pulse time and update the valid indices
for i = 1:length(artifact_times)
    exclude_start = artifact_times(i);
    exclude_end = artifact_times(i) + exclude_durations;
    valid_indices_hdh = valid_indices_hdh & ~(hdh_mt >= exclude_start & hdh_mt <= exclude_end);
    % valid_indices_hh1 = valid_indices_hh1 & ~(hh1_mt >= exclude_start & hh1_mt <= exclude_end);
    % valid_indices_hh2 = valid_indices_hh2 & ~(hh2_mt >= exclude_start & hh2_mt <= exclude_end);
    % valid_indices_hhz = valid_indices_hhz & ~(hhz_mt >= exclude_start & hhz_mt <= exclude_end);
end

% Define the gap period
gap_start = datetime('1-July-2020 00:00:00');
gap_end = datetime('14-July-2020 14:30:00');

% Update valid indices to exclude the gap period
valid_indices_hdh = valid_indices_hdh & ~(hdh_mt >= gap_start & hdh_mt <= gap_end);

% Filter out the artifact time ranges from the data
t_valid_hdh = hdh_mt(valid_indices_hdh);
d_valid_hdh = hdh_md(valid_indices_hdh);

d_valid_hhz = hhz_md(valid_indices_hdh);
d_valid_hh1 = hh1_md(valid_indices_hdh);
d_valid_hh2 = hh2_md(valid_indices_hdh);

t_interp = (hdh_mt(1):seconds(1):hdh_mt(end));
d_interp = interp1(t_valid_hdh, d_valid_hdh, t_interp);

d_interp_hhz = interp1(t_valid_hdh, d_valid_hhz, t_interp);
d_interp_hh1 = interp1(t_valid_hdh, d_valid_hh1, t_interp);
d_interp_hh2 = interp1(t_valid_hdh, d_valid_hh2, t_interp);

% Apply the bandpass filter to the data, demean and detrend the filtered data
hdh_demeaned = d_interp - mean(d_interp);
hdh_detrended = detrend(hdh_demeaned);

hhz_demeaned = d_interp_hhz - mean(d_interp_hhz);
hhz_detrended = detrend(hhz_demeaned);

hh1_demeaned = d_interp_hh1 - mean(d_interp_hh1);
hh1_detrended = detrend(hh1_demeaned);

hh2_demeaned = d_interp_hh2 - mean(d_interp_hh2);
hh2_detrended = detrend(hh2_demeaned);


% W=hamming(1,length(hdh_detrended));
% hdh_tapered = hdh_detrended.*W';

% Design the Butterworth bandpass filter

% Define the sampling frequency (Fs) and cutoff frequency
Fs = 1; % Sampling frequency (1 Hz in this case)
fnyq = Fs / 2;

% Bandpass Butterworth filter design
% Fc1 = 1 / (3600*26);
% Fc2 = 1 / (3600*9);
% [b, a] = butter(2, [Fc1 Fc2] / fnyq, 'bandpass'); % 2nd order Butterworth filter
% hdh_filtered_twice = filter(b, a, hdh_tapered);

% low pass
% Fc = 1/(3600*8);
% [b, a] = butter(2, Fc / fnyq, 'low'); % 2nd order Butterworth filter

%% Bandpass Butterworth filter design
clc
Fc1 = 1 / 100;
Fc2 = 1 / 20;
[b, a] = butter(2, [Fc1 Fc2] / fnyq, 'bandpass'); % 2nd order Butterworth filter
hdh_filtered_demeaned = filtfilt(b, a, hdh_detrended);
hhz_filtered_demeaned = filtfilt(b, a, hhz_detrended);
hh1_filtered_demeaned = filtfilt(b, a, hh1_detrended);
hh2_filtered_demeaned = filtfilt(b, a, hh2_detrended);
%% Plotting seismic data for specific waves

close all; clc
% time_x = [datetime(2020,3,25,2, 57, 0), datetime(2020,3,25, 2, 58, 0)];
time_x = [datetime(2020,5,9,4, 0, 0), datetime(2020,5,9, 5, 30, 0)];
t_cut_idx = (t_valid_hdh > time_x(1)) & (t_valid_hdh < time_x(2));
t_cut = t_valid_hdh(t_cut_idx);

figure();

nn = 320;
con_window=sin(pi*(1:nn)./nn);
con_window=con_window./sum(con_window);
%old_con_window=ones(1,nn)/nn;

hdh_conv = conv(d_valid_hdh,con_window,'same');
hhz_conv = conv(d_valid_hhz,con_window,'same');
% hh1_conv = conv(d_valid_hh1,con_window,'same');
% hh2_conv = conv(d_valid_hh2,con_window,'same');

obs_correction = 421192.92*754.3; % from Heather's email

% Filtered data plot
ax(1) = subplot(2,1,1);
% hdh_det = detrend((hdh_conv(t_cut_idx)-mean(hdh_conv(t_cut_idx)) ) / obs_correction );
% plot(t_cut, hdh_det * 6894.76, 'DisplayName', 'HDH convolved (Pa)');
% hold on
% xlim(time_x);
% legend;
% grid on;
% plot(t_interp, cumsum(hdh_filtered_demeaned) );
% hold on
% xlim(time_x);
% legend('HDH demeaned, detrended, filtered')
% grid on;

% ax(2) = subplot(2,1, 2);
% hhz_velocity = detrend((hhz_conv(t_cut_idx)-mean(hhz_conv(t_cut_idx)) ) / obs_correction );
% plot(t_cut, hhz_velocity, 'DisplayName', 'HHZ convolved') % we're not integrating over time because dt = 1 Hz
% plot(t_valid_hdh, d_valid_hhz)
% plot(t_interp, cumsum(hhz_filtered_demeaned) );
% hold on
% xlim(time_x);
% legend('HHZ demeaned, detrended, filtered')
% grid on;

ax(1) = subplot(2,1,1);
hold on;
hdh_plot = detrend((hdh_conv(t_cut_idx)-mean(hdh_conv(t_cut_idx)) ) );
plot(t_cut, hdh_plot )
legend('HDH convolved')
xlim(time_x);
legend;
grid on;

ax(2) = subplot(2,1, 2);
hhz_velocity = detrend((hhz_conv(t_cut_idx)-mean(hhz_conv(t_cut_idx)) ) / obs_correction );
plot(t_cut, cumsum(hhz_velocity) )
legend('HHZ convolved & integrated')
xlim(time_x);
legend;
grid on;

linkaxes([ax], 'x')

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
%old_con_window=ones(1,nn)/nn;

my_data = conv(d_interp,con_window,'same');
% my_data = d_interp;
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
    wave_passes = wave_passes - mean(wave_passes); % demean data
    
    % Calculate period, max amp, and min amp
    if ~isempty(wave_passes)
        [amp_max, amp_idx] = max(abs(wave_passes));
        [crest, crest_idx] = max(wave_passes);
        before_peak = wave_passes(1:crest_idx);
        [trough, trough_idx] = min(before_peak);
        
        % Sign of signal
        wave_sign = sign(wave_passes);
        
        % Find where signal changes sign (zero crossings)
        zero_crossings = find(diff(wave_sign) ~= 0);
        
        % T1: last zero before trough
        T1_candidates = zero_crossings(zero_crossings < trough_idx);
        if ~isempty(T1_candidates)
            T1_idx = T1_candidates(end);
            T1 = t_wave(T1_idx);
        else
            T1 = NaT;
        end
        
        % T2: first zero after crest
        T2_candidates = zero_crossings(zero_crossings > crest_idx);
        if ~isempty(T2_candidates)
            T2_idx = T2_candidates(1);
            T2 = t_wave(T2_idx);
        else
            T2 = NaT;
        end

    end

    % Store the results
    amps(i) = amp_max;
    amp_time(i) = t_wave(amp_idx);
    crests(i) = crest;
    crest_time(i) = t_wave(crest_idx);
    troughs(i) = trough;
    trough_time(i) = t_wave(trough_idx);
    T1s(i) = T1';
    T2s(i) = T2';
end
%%
detection.amps = amps';
detection.amp_time = amp_time';
detection.crests = crests';
detection.crest_time = crest_time';
detection.troughs = troughs';
detection.trough_time = trough_time';
detection.template = month_templates;
detection.T1s = T1s';
detection.T2s = T2s';
detection.period = (crest_time-trough_time)';
%% Check each wave manually

% for f = 1:height(detection)
   
for f = 376
    % close all;

    time_x = [velocity_window_start(f), velocity_window_end(f)];
    % time_x = [sorted_detections.DetectionTime(f)-hours(window / 4),sorted_detections.DetectionTime(f)+hours(window / 2)]; 
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
    scatter(detection.crest_time(f), detection.crests(f), 'filled', 'r')
    scatter(detection.trough_time(f), detection.troughs(f), 'filled', 'r')
    yline(0);
    % scatter(detection.T1s(f), 0, 'filled', 'b')
    % scatter(detection.T2s(f), 0, 'filled', 'g')
    legend('HDH', 'Crest', 'Trough','y=0')%,'T1', 'T2')
    % title('Ocean-Bottom Pressure');
    % if ~isnat(detection.template(f))
    %     subtitle(['Template: ', datestr(detection.template(f), 'yy-mm-dd HH:MM:SS')]);
    % else
    %     subtitle('Template: NaT')
    % end
    % xlabel('Time');
    % ylabel('Amplitude (psi)');
    text(velocity_window_start(f) + minutes(10), 0.01, ...
    ['period = ' num2str(minutes(detection.period(f)), '%.2f') ' minutes'])

    grid on
    xlim(time_x);
    % ylim([-0.07 0.07])
    f
    detection.Type(f)
    grid on;
end


%% Remove signals that may not be internal waves
bad_indices_original = [41 46 62 63 64 65 67 70 75 78 80 87 97 99 122 134 135 163 164 178 179 183 197 201 208 209 212 236 237 242 246 249 254 258 259 260 270 277 290 291 293 301 309 315 366 401 420 426 443 444 455 460];
bad_indices2 = [41 46 54 62 63 64 65 67 75 78 86 87 97 99 135 164 178 183 197 209 259 260 284 330 401 420 443 444 455] ;
% % filtered_waves = detection_window;
filtered_detections = detection;
filtered_detections(bad_indices2, :) = []; % Removes the row(s)
% remove aug 11 6.57 pm
writetable(filtered_detections, 'all_detections4.csv');
