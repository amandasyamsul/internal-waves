addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

obs_coords = [21.00116 117.40267];
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

% t_interp = (hdh_mt(1):seconds(1):hdh_mt(end));
t_tides = [filtered_waves.detection_time(1):hours(10):filtered_waves.detection_time(end)];
z_2020 = tmd_predict(tide_model, lat, lon, t_tides,'h');



%%
close all;

% Z=fft(z_2020);
% N=length(Z);
% dt=1;
% figure(10)
% f=[0:N-1]/(N*dt);
% T=1./f/3600;
% T=T(1:end/2);
% Z=Z(1:end/2);
% 
% loglog(T, abs(Z));
% xlabel('Period (hours)')
% title('OTPS model')

figure(1)

NFFT=4194304; % 2^log2(3600*24*28)
[pxx_model, f_model] = pwelch(z_2020,[],[],NFFT);
f_model_hz=f_model/(2*pi);
T_model=1./f_model_hz/3600;

loglog(T_model, abs(pxx_model) )
grid on;
xlabel('Period (hours)')
title('OTPS model')

% OBP data
figure(2)

detrended_data = detrend(hdh_md);
[pxx, f] = pwelch(detrended_data,[],[],NFFT);
f_hz=f/(2*pi);
T=1./f_hz/3600;
loglog(T, abs(pxx))
grid on
xlabel('Period (hours)')
title('OBP data')

%% filtered data
figure(30)

fd = hdh_processed;
[pxx2, f2] = pwelch(fd,[],[],NFFT);
f2_hz=f2/(2*pi);
T2=1./f2_hz/3600;
loglog(T2, abs(pxx2))
grid on
xlabel('Period (hours)')
title('OBP filtered data')

%%
figure(4)

fs = 1; 
[S, F, T] = spectrogram(fd);

imagesc(T, F, 10*log10(abs(S).^2));
axis xy;
c = colorbar();
c.Label.String = 'Power density';

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of OBP filtered data');

% figure(5)
% waterplot(S,F,T)
%%
Z_obp=fft(d_interp);
N=length(Z_obp);
figure(20)
f=[0:N-1]/(N*dt);
T=1./f_hz/3600;
T=T(1:end/2);
Z_obp=Z_obp(1:end/2);

loglog(T, abs(Z_obp));
xlabel('Period (hours)')
title('OBP data')


function waterplot(s,f,t)
% Waterfall plot of spectrogram
    waterfall(f,t,abs(s)'.^2)
    set(gca,XDir="reverse",View=[30 50])
    xlabel("Frequency (Hz)")
    ylabel("Time (s)")
end