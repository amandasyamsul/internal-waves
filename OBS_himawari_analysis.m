% Amanda Syamsul
% March 7th, 2025
% Analysis on Internal Wave propagation proximal to OBS
% Relevant outputs: Figs. 2b, 4ab, 5 (without model), 6a/b/c (needs d), S6

close all; clear all; clc

data_dir = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(data_dir);

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/data/');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/seismic');
addpath('/Users/amandasyamsul/Documents/MATLAB/DrosteEffect-BrewerMap-3.2.5.0')
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

N2 = readtable('daily_max_N2_all.csv');
filtered_waves = readtable('filtered_waves.csv');
all_detections = readtable('all_detections5.csv');
luzon = readtable('luzon.csv');

all_data_2020 = readtable('OIW_2020_data.csv');
OBS_latband_data_2020 = readtable('IW_OBS_latband_2020.csv');
avglat = readtable('avg_OBS_latband_data_2020.csv');
OBS_data_2020 = readtable('IW_OBS_2020.csv'); % square around OBS

OBS_latband_data_2020_adjusted = readtable('OBS_latband_2020_adjusted.csv');
fs = 13;

% Choose which datasets to use for this analysis
square_around_OBS = OBS_data_2020;
latband_across_OBS = OBS_latband_data_2020_adjusted;

dataset = latband_across_OBS;
plot_year = 2020;
obs_coords = [21.00116 117.40267];

% Define geographical bounds (longitude and latitude)
lon1 = 116.2; lon2 = 117.7;  
lat1 = 20.0; lat2 = 21.5;   

% Define x-y coordinate system bounds (based on satellite image pixels)
x_min = 0; x_max = 500;      % X-coord in image space
y_min = -400; y_max = 0;     % Y-coord in image space

% Computte the range of coordinates
image_width = x_max - x_min;    % Total width of the image in x-coordinates
image_height = y_max - y_min;   % Total height of the image in y-coordinates

lon_obs = lon1 + (square_around_OBS.x_coord - x_min) / image_width * (lon2 - lon1);
lon_obs_latband = lon1 + (latband_across_OBS.x_coord - x_min) / image_width * (lon2 - lon1);
square_around_OBS.longitude = lon_obs;
latband_across_OBS.longitude = lon_obs_latband;

bathymetry1 = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "elevation");
bathymetry = bathymetry1';
lon = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "lon");
lat = ncread("gebco_2024_n21.5_s20.0_w116.2_e117.7.nc", "lat");

% Find indices within the latitude band
lat_idx = lat >= 21.0012 & lat <= 21.1012;

% Subset bathymetry to only those latitudes
bathymetry_latband = bathymetry(lat_idx, :);

% Now compute the average along latitude dimension (rows)
avg_depth_along_lon = mean(bathymetry, 1, 'omitnan'); 
avg_depth_along_lon_latband = mean(bathymetry_latband, 1, 'omitnan'); 

% alldet_small = all_detections(:, {'crest_time','crests', 'trough_time', 'troughs', 'period'});  
% writetable(alldet_small, 'OBP_measurements.csv');

%% Manuscript Figure 2b: Plot average depth versus longitude ✅
% close all;

figure();
plot(lon, avg_depth_along_lon_latband, 'LineWidth', 2);
hold on
plot(lon, avg_depth_along_lon, 'LineWidth', 2);
xline(obs_coords(2), 'LineWidth',2)
xline(117.0792445, 'k--', 'LineWidth',2)
legend('Latitude 21.0012\circN to 21.1012\circN','Latitude 20.0\circN to 21.5\circN', 'OBP', 'Cutoff longitude')
xlabel('Longitude', 'FontSize', fs)
ylabel('Average depth (m)', 'FontSize', fs)
% title('Average Depth Along Longitude', 'FontSize', fs, 'FontAngle', 'italic');
grid on;

%% Manuscript Figure S6: Correlating propagation speed with amplitude & period
% close all; 
clc;

start_date = datetime(plot_year, 5, 1); 
end_date = datetime(plot_year, 8, 31); 
range = (dataset.datetime >= start_date) & (dataset.datetime <= end_date);
data_lon = dataset.longitude(range);
data_velocity = dataset.velocity(range);

% Match longitudes to corresponding average depths
ccc = knnsearch(lon(:), data_lon(:)); 
depths_at_data_lon = avg_depth_along_lon(ccc); % Get average depths

% Bin velocities by depth
bin_edges = -1000:25:0; % depth bins
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[N, edges, bin_indices] = histcounts(depths_at_data_lon, bin_edges);

% Average velocity per depth bin
avg_velocity = arrayfun(@(b) mean(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));
std_velocity = arrayfun(@(b) std(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));

ttt = dataset(range, :);
ttt_speed = ttt.velocity;
ttt.amplitude = NaN(height(ttt), 1);
ttt.period = NaN(height(ttt), 1);
unique_templates = unique(all_detections.template);

max_time_diff = hours(4); 
matched_amps = NaN(height(ttt), 1);
matched_periods = NaN(height(ttt), 1);

figure()
subplot(1,2,1)
add_reference_lines()
for t = 1:length(unique_templates)
    
    amps_for_this_template = NaN(height(ttt), 1);
    periods_for_this_template = NaN(height(ttt), 1);

    this_template = unique_templates(t);
    filtered_waves_by_temp = all_detections(all_detections.template == this_template, :);
    
    for i = 1:height(ttt)
        time_diffs = abs(filtered_waves_by_temp.DetectionTime - ttt.datetime(i));
        [min_diff, closest_idx] = min(time_diffs);
    
        if min_diff <= max_time_diff
            matched_amps(i) = filtered_waves_by_temp.amps(closest_idx);
            matched_periods(i) = seconds(filtered_waves_by_temp.period(closest_idx));

            amps_for_this_template(i) = filtered_waves_by_temp.amps(closest_idx);
            periods_for_this_template(i) = seconds(filtered_waves_by_temp.period(closest_idx));

            % matched_amps is a cumulative table of all templates, while 
             % amps_for_this_template is reset every iteration
        
        end
    end
    
    hold on

    for b = 1:length(bin_centers)
        idx = bin_indices == b;
    
        velocities = ttt_speed(idx);
        amplitudes = amps_for_this_template(idx);
    
        valid = ~isnan(velocities) & ~isnan(amplitudes);
        velocities = velocities(valid);
        amplitudes = amplitudes(valid);
    
        jitter = (rand(size(velocities)) - 0.5) * 5;
        scatter(bin_centers(b) + jitter, velocities, amplitudes/3, amplitudes, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    xlim([-1000 -500]);
    ylim([0 5]);
    xlabel('elevation (m/s)')
    ylabel('propagation speed (m/s)')
    % cb = colorbar('southoutside');                 
    % cb.Label.String = 'Amplitude (psi)';
    caxis([100 500]);
    colormap(brewermap([],"YlOrRd"))
    grid on;
    set(gca, 'XDir','reverse')
    hold off;
end

ttt.amplitude = matched_amps;
ttt.period = matched_periods;

% PART 3B: Day-by-day moving average
window_size=5;
% Step 1: Get just the date part (no time)
dates_only = dateshift(ttt.datetime, 'start', 'day');

% Step 2: Get unique dates and group indices
[unique_dates, ~, day_group_ids] = unique(dates_only);

avg_velocity    = [];
avg_amplitude   = [];
avg_depth       = [];

% Loop over each day
for d = 1:length(unique_dates)
    idx_today = find(day_group_ids == d);
    
    % Skip if too few points for window
    if length(idx_today) < window_size
        continue
    end
    
    % Extract data for this day
    v_day = ttt_speed(idx_today);
    a_day = ttt.amplitude(idx_today);
    d_day = depths_at_data_lon(idx_today);
    
    valid_len = length(idx_today) - window_size + 1;
    
    avg_v_day = NaN(valid_len,1);
    avg_a_day = NaN(valid_len,1);
    avg_d_day = NaN(valid_len,1);

    for i = 1:valid_len
        v_window = v_day(i:i+window_size-1);
        a_window = a_day(i:i+window_size-1);
        d_window = d_day(i:i+window_size-1);
        
        avg_v_day(i) = nanmean(v_window);
        avg_a_day(i) = nanmean(a_window);
        avg_d_day(i) = nanmean(d_window);
    end

    avg_velocity    = [avg_velocity; avg_v_day];
    avg_amplitude   = [avg_amplitude; avg_a_day];
    avg_depth       = [avg_depth; avg_d_day];
end

% Bin by depth
bin_edges = -1000:25:0;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[~, ~, bin_indices] = histcounts(avg_depth, bin_edges);

subplot(1,2,2)
add_reference_lines()
xlabel('elevation (m/s)')
ylabel('propagation speed (m/s)')
hold on
plot_moving_averages(bin_centers,bin_indices,avg_velocity,avg_amplitude)
set(gca, 'XDir','reverse')

%% Plot propagation speed binned by depth
% close all;

fig = figure();
subplot(1,2,1)
hold on;
for b = 1:length(bin_centers)
    % Get all propagation speed values in the bin
    velocities_in_bin = ttt_speed(bin_indices == b);
    
    jitter = (rand(size(velocities_in_bin)) - 0.5) * 5; % Adjust jitter magnitude as needed
    scatter(bin_centers(b) + jitter, velocities_in_bin, 10, 'filled');
    hold on
end
xline(-900, 'k--', 'LineWidth',2)
xlabel('Depth (m)');
ylabel('propagation speed (m/s)');
set(gca, "XDir", "reverse")
title('Propagation speeds binned by depth')
subtitle('Latitude 21.0012 - 21.1012')
ylim([0 5])
grid on;
hold off;

subplot(1,2,2)
% Bin velocities by depth
bin_edges = -1000:25:0; % depth bins
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
[N, edges, bin_indices] = histcounts(depths_at_data_lon, bin_edges);

% Average velocity per depth bin
avg_velocity = arrayfun(@(b) mean(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));
std_velocity = arrayfun(@(b) std(data_velocity(bin_indices == b), 'omitnan'), 1:length(bin_centers));

hold on;
errorbar(bin_centers, avg_velocity, std_velocity, 'k.', 'CapSize', 0, 'HandleVisibility', 'off');
scatter(bin_centers, avg_velocity, 80, 'filled');
xline(-900, 'k--', 'LineWidth',2)
% title(num2str(years(i)), 'FontSize', fs, 'FontAngle', 'italic');
xlim([-1000, -400])
ylim([0,5])
grid on;
title('Avg. propagation speeds binned by depth')
subtitle('Latitude 21.0012 - 21.1012')
xlabel('Depth (m)');
ylabel('propagation speed (m/s)');

set(gca, "XDir", "reverse")

%% Manuscript Figure 4 a,b: T and A over time 
% close all;
obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);
figure()

ax1 = subplot(2,1,1); hold on

www  = unique(all_detections.Type);     % still a cell array of char
cols = lines(numel(www));               % colors for each type
shapes = {'o','s','d','^','v','>'};     % circle, square, diamond, triangle up, down, right

for k = 1:numel(www)
    idx = strcmp(all_detections.Type, www{k});   % indices for this type
    
    % plot with specific marker and color
    plot(all_detections.DetectionTime(idx), seconds(all_detections.period(idx)), ...
        shapes{mod(k-1, numel(shapes))+1}, ...
        'MarkerFaceColor', cols(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', www{k});
end

grid on
% ylabel('Period (s)')
legend('Location','northeast')

ax2 = subplot(2,1,2); hold on
for k = 1:numel(www)
    idx = strcmp(all_detections.Type, www{k});
    
    plot(all_detections.DetectionTime(idx), all_detections.amps(idx), ...
        shapes{mod(k-1, numel(shapes))+1}, ...
        'MarkerFaceColor', cols(k,:), ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', www{k});
end

grid on
% ylabel('Amplitude (psi)')
legend('Location','northeast')

linkaxes([ax1 ax2], 'x')

%% Manuscript Figure 5: T vs A
% close all;

clc;
all_detections = readtable('all_detections5.csv');

data_to_bin = all_detections.period; % convert to s
ampli = all_detections.amps;

% Define bin edges
startT = min(data_to_bin);
endT   = max(data_to_bin);
edgesT = startT:minutes(1):endT;

[N,~,binIdx] = histcounts(data_to_bin, edgesT);

amps_avg = accumarray(binIdx((binIdx>0)), ampli((binIdx>0)), [], @mean);
amps_std = accumarray(binIdx(binIdx>0), ampli(binIdx>0), [], @std);
amps_count = accumarray(binIdx(binIdx>0), 1, [], @sum);
amps_sem = amps_std ./ sqrt(amps_count);  % standard error

binDates = edgesT(1:end-1)';  % start of each bin

Ts = seconds(edgesT(1:end-1));

figure()
% unbinned version
Tvals = seconds(data_to_bin);
Avals = ampli;
% Fit power law: log A ~ alpha log T + const
% p_unbin = polyfit(log10(Tvals), log10(Avals), 1);
legend()
loglog(Tvals, Avals, '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 18, 'DisplayName','Unbinned T vs A');
hold on
p_unbin = polyfit(log10(Tvals), log10(Avals), 1)

hold on
errorbar(Ts, amps_avg, amps_sem, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 8, ...
         'CapSize', 6, 'LineWidth', 2, 'DisplayName','Binned T vs A')  % error bars
xlabel('Period (s)')
ylabel('Amplitude (Pa)')
grid on
p = polyfit(log10(Ts), log10(amps_avg), 1)
xfit = linspace(min(Ts), max(Ts), 100);
yfit = 10.^(polyval(p, log10(xfit)));

loglog(xfit, yfit, 'b--','LineWidth',3, 'DisplayName','Slope = -1.34')


x_unbin = log10(Tvals(:));
y_unbin = log10(Avals(:));

lm_unbin = fitlm(x_unbin, y_unbin);   % linear model
disp(lm_unbin)

% Add in nonhydrostatic model for A vs T

g=9.8; % gravity - never change
rho=1e3; % density of water - never change
H=619; % total depth 

h1_grid=[80:5:140];
coef = 0.0002;
a_grid = [-180:10:-100];
cmap = cool(length(a_grid));  
% figure()

for aa = 1:length(a_grid)
    a = a_grid(aa);

    for hh = 1:length(h1_grid)
        h1 = h1_grid(hh);
        h2=H-h1; % lower layer
        gp=g*coef; 
    
        s=[-1.5e3:1.5e3]; % x coordinate of the region (meters)
        [c0store,kstore,alpha,beta,I4,I5]=ComputeCoefficients(h1,h2,gp,g); 
        % then pressure
        [pH,pSH,pNH,L,T,c1,c0]=ComputePressure(s,h1,h2,a,gp,g,rho,kstore,I4,I5,c0store,alpha,beta);
        Total= pH+pSH+pNH;

        T_grid(hh) = T; 
        A_grid(hh) = max(Total);
    
    end
   
    loglog(T_grid, A_grid, 'o', ...
        'MarkerSize',10, ...
        'MarkerFaceColor',cmap(aa,:), ...
        'MarkerEdgeColor',cmap(aa,:), ...
        'DisplayName',sprintf('a = %g m', a_grid(aa)))
    hold on
    grid on
end

xlim([250 900])
ylim([0 1e3])
legend()


%% Manuscript Figure 6: Plot fortnightly cycle
% close all;
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/OTPS/chadagreene-Tide-Model-Driver-fd89bd6');

obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

% OTPS (unchanged)
t_otps = (all_detections.DetectionTime(1):hours(1):all_detections.DetectionTime(end));
z_otps = tmd_predict(tide_model, lat, lon, t_otps, 'h');

startDate1 = min(all_detections.DetectionTime) - days(2);
endDate1   = datetime(2020, 2, 27);
edges1 = startDate1:days(13.66):all_detections.DetectionTime(end);

figure()
ax(1) = subplot(4,1,1); hold on
plot(t_otps, z_otps, 'Color',	[0.6350, 0.0780, 0.1840])
xline(edges1, 'k')
ylabel('Tidal height (m)')

ax(2) = subplot(4,1,2); hold on
scatter(all_detections.DetectionTime, seconds(all_detections.period), 12, 'k', 'filled')
xline(edges1, 'k')
ylabel('Period (s)')

data_to_bin = seconds(all_detections.period);

% Bin once using the combined edges (left-inclusive to avoid double-counting)
[~,~,binIdx] = histcounts(all_detections.DetectionTime, edges1);
valid = (binIdx > 0) & ~isnan(data_to_bin);
fortnightly_avg = accumarray(binIdx(valid), data_to_bin(valid), [numel(edges1)-1, 1], @mean, NaN);

binDates = edges1(1:end-1)';  % start of each bin

ax(3) = subplot(4,1,3); hold on
scatter(binDates, fortnightly_avg, 28, 'k', 'filled')
plot(binDates, fortnightly_avg, '-')
grid on
ylabel('Average period (s)')

ax(4) = subplot(4,1,4); hold on
scatter(N2.date, N2.max_N2, 'filled')
grid on
ylabel('Max N_2')

linkaxes(ax, 'x')

%% Functions

function plot_moving_averages(bin_centers,bin_indices,avg_velocity,avg_amplitude)
    hold on;
    for b = 6:length(bin_centers)
        idx = bin_indices == b;
    
        v = avg_velocity(idx);
        a = avg_amplitude(idx);
    
        % Jitter for depth
        jitter = (rand(size(v)) - 0.5) * 5;
        scatter(bin_centers(b) + jitter, v, a/3, a, 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    xlim([-1000 -500]);
    ylim([1 3]);        
    cb.Label.String = 'Amplitude (Pa)'; 
    caxis([100 500]);
    colormap(brewermap([],"YlOrRd"))
    grid on;
end

function add_reference_lines()
% Mean speed on continental slope from Ramp (2010)
    y0 = 2.22;
    dy = 0.18;
    
    ymin = y0 - dy;
    ymax = y0 + dy;
    
    xmin = -2491;
    xmax = -350;
    
    patch([xmin xmax xmax xmin], [ymin ymin ymax ymax], [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
end

function bubblelegend(sizeVals, labels, color, location)
% BUBBLELEGEND Create a size legend for scatter/bubble plots.
%
% bubblelegend(sizeVals, labels, color, location)
%   sizeVals - vector of marker sizes (same units as scatter SizeData)
%   labels   - cell array of text labels (same length as sizeVals)
%   color    - RGB triplet or color name (e.g., 'k', [0 0 0])
%   location - legend location string ('northwest','northeast', etc.)
%
% Example:
%   bubblelegend([20 40 60], {'Small','Medium','Large'}, 'k', 'northeast')

    hold on
    ax = gca;
    xL = xlim(ax);
    yL = ylim(ax);

    xPos = xL(2) + 0.02 * range(xL);
    yPosStart = yL(2) - 0.05 * range(yL);

    for i = 1:length(sizeVals)
        yPos = yPosStart - (i-1) * 0.07 * range(yL);
        scatter(xPos, yPos, sizeVals(i), color, 'filled');
        text(xPos + 0.03 * range(xL), yPos, labels{i}, ...
            'VerticalAlignment','middle');
    end
    ax.XLim = [xL(1) xL(2) + 0.2*range(xL)];
end

%% Model functions
function [c0store,kstore,alpha,beta,I4,I5]=ComputeCoefficients(h1,h2,gp,g)
    % computes c0,alpha, beta for each wavenumber k for free surface case
    H=h1+h2;
    dz=0.1;
    z=[0:dz:H];
    ii=0;
    for ik=-4:0.01:-1, % this creates steps, and we can smooth out using 0.001
        ii=ii+1;
        k=10^ik;
        kstore(ii)=k;
        
        % For each wavenumber k, we compute the linear wave speed c0 that a small-amplitude internal wave would travel at in a two-layer ocean with a free surface
        % hyperbolic tangents for later
        T1=tanh(k*h1);
        T2=tanh(k*h2);
    
        % Ac^4+Bc^2+C=0 (quadratic eqn where x=c^2)
        A=k^2*(1+T1*T2);
        B=-k*(g*(T1+T2)+gp*T2);
        C=g*gp*T1*T2;
        c02=(-B-sqrt(B^2-4*A*C))/(2*A);  % + is surface. - is internal
        c0=sqrt(c02); 
        c0store(ii)=c0;
        beta1=-g/(k*c02); % this is not the same as beta (dispersion)
        
        % phi = vertical structure (how horizontal v changes with depth)
        phi1 = (cosh(k*(H-z))+beta1*sinh(k*(H-z)))./(cosh(k*h1)+beta1*sinh(k*h1)); % top layer
        phi2 = sinh(k*z)./(sinh(k*h2)); % bottom layer
        II2=(z<=h2); % idx for bottom layer
        II1=(z>h2); % idx for top layer
        phi=[phi2(II2) phi1(II1)];
    
        % z at the end means dphi/dz
        phi1z=(-k*sinh(k*(H-z))-beta1*k*cosh(k*(H-z)))./(cosh(k*h1)+beta1*sinh(k*h1)); % analytical derivative
        phi2z=k*cosh(k*z)./(sinh(k*h2));
        phiz=[phi2z(II2) phi1z(II1)];
        
        % integrals
        I1=trapz(phiz.^3)*dz;
        I2=trapz(phiz.^2)*dz;
        I3=trapz(phi.^2)*dz;
        
        % I4 and I5 are for pressure calculations
        I4(ii)=trapz(phi)*dz;
        I5(ii)=trapz(phi.*phiz)*dz;
        
        % alpha = nonlinearity, beta = dispersion
        % we are assuming negligible current velocity (ie. U=0)
        alpha(ii)=(3/2)*c0*I1/I2; % note sign switch from Moum and Nash U-c0 
        beta(ii)=(1/2)*c0*I3/I2;
    end
end

function [pH,pSH,pNH,L,T, c1, c0]=ComputePressure(s,h1,h2,a,gp,g,rho,kstore,I4,I5,c0store,alpha,beta)
    H=h1+h2;
    L_array=sqrt(12*beta./(a.*alpha)); % equation before 8.17 Gerkema 
    [~,ik]=min(abs(1./kstore-L_array));
    L=L_array(ik); % wavelength
    k=1/L;
    
    % Shape function
    sech_sL = 1./cosh(s/L);
    A = a * sech_sL.^2;
    As = (-2*a/L) .* sech_sL.^2 .* tanh(s./L);
    Ass = (4*a/L^2) * sech_sL.^2 - (6*a/L^2)*sech_sL.^4;
    
    % shallow-water dispersion !! (simplification consistent with calculation of L)
    %Bref=-(g*H+gp*h2);
    %Cref=g*gp*h1*h2;
    %c0 = 0.5*(-Bref - sqrt(Bref^2 - 4*Cref));
    % no longer shallow water
    c0=c0store(ik);
    %c1shallow = c0 * (1 + a*(h1-h2)/(2*h1*h2));
    c1=c0+a*alpha(ik)/3; % c1 = nonlinear velocity, a = amplitude
    %[c1 c0 c1shallow]
    
    % Pressure
    I4_k=I4(ik);
    I5_k=I5(ik);
    
    D = cosh(k*h1) + (-g/(c0^2*k)) * sinh(k*h1); % required to incorporate free surface
    pH = rho * A * gp;
    pSH = rho * A * g/D;
    pNH = rho * c1^2 * (Ass * I4_k );%+ (As.^2 - A.*Ass)*I5_k); % negligble last term (advective)
    T=L/c1;
end
