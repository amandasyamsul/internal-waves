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
%  Checking if OBP sign same as tidal height
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

hdh = load('bb01_dec_1Hz_all_HDH.mat');
m = [1:2]; % month

% Convert the times from datenum to datetime
hdh_times = datetime(hdh.t, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% 
% Filter data to match desired timeframe
hdh_month_idx = (year(hdh_times) == 2020) & ismember(month(hdh_times), m);
hdh_mt = hdh_times(hdh_month_idx);
hdh_md = hdh.d(hdh_month_idx);   
hdh_md = hdh_md/1676128.86; % correction factor to convert pressure to psi

t_otps = (hdh_mt(1):seconds(1):hdh_mt(end));
z_otps = tmd_predict(tide_model, lat, lon, t_otps,'h');

%% Decimate & filter OBP data (retain only 1/10 data points)
close all;
decim = 10;
d = hdh.d / 1676128.86; % correction factor to convert pressure to psi
t = hdh.t;

d2=d(1:end/2);
d2=detrend(d2);
t2=t(1:end/2);


d3=decimate(d2,decim);
t3 = t2(1:decim:end); 
t3_datetime = datetime(t3, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');

figure(10)
plot(t3_datetime,d3)
grid on
title('decimated obp data')

dt=1/(3600/decim);
maxT=35;
minT=11;
[b,a]=butter(2, [2*dt/maxT 2*dt/minT],'bandpass');
dbp=filtfilt(b,a,d3);

figure(11)
plot(t3_datetime,dbp)
grid on
title('decimated & lowpassed obp data')

%% Decimate & filter OTPS data
dotps = z_otps;
totps = t_otps;

d2otps = detrend(dotps);

d3otps = decimate(d2otps,10);
t3otps = totps(1:10:end); 

figure(20)
plot(t3otps,d3otps)
grid on
title('decimated otps model')

dbp_otps=filtfilt(b,a,d3otps);

figure(21)
plot(t3otps,dbp_otps)
grid on
title('decimated & lowpassed otps model')

%%
close all;
time_x = [datetime(2020,1,1), datetime(2020,2,30)];

figure(222);

% Filtered data plot
ax2 = subplot(2, 1, 1);
hold on;
plot(t3_datetime, dbp);
title('HDH decimated & filtered');
yline(0, 'r')
xlabel('Time');
ylabel('Amplitude (psi)');
xlim(time_x);
grid on;

% Filtered data plot

ax3 = subplot(2, 1, 2);
hold on;
plot(t3otps,dbp_otps)
title('OTPS model decimated & filtered');
yline(0, 'r')
xlabel('Time');
ylabel('Tidal height (m)');
xlim(time_x);
% datetick
grid on;

linkaxes([ax2 ax3], 'x')
%%
close all;

figure()
plot(t3otps,dbp_otps,t3_datetime,dbp*1e2)
grid on;
%%
% 
% t_gauge = [datetime(2025,6,1):hours(1):datetime(2025,6,15)];
% z_gauge = tmd_predict(tide_model, 13.7, 109.250, t_gauge,'h');
% 

%%
function ybp=low(x,dt,T)
    % Usage:  bp(x, dt,T)
    % lowpasses signal from T secs
    % Using 1 pass 4th order butterworth filter
    
    [b,a]=butter(4,2*dt/T);
    ybp=filter(b,a,x);
end
 