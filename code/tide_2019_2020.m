function [l_time, l_vel_u, l_vel_v, d_time, d_height] = tide_2019_2020()

%% Load Luzon E-W tidal velocity at 20.6000/121.9000
fid = fopen('Luzon_2019_2020.txt'); %UTC time 
A = textscan(fid, '%s %s %s %s %s %s %s', 'Headerlines', 1);
d_tide = str2double(split(A{3}, '.'));
t_tide = str2double(split(A{4},':'));
l_time = datenum(d_tide(:,3), d_tide(:,1), d_tide(:,2), t_tide(:,1), t_tide(:,2), t_tide(:,3));
l_vel_u = str2double(A{5}); %cm/s; E-W
l_vel_v = str2double(A{6}); %cm/s; N-S

%% Load Dongsha tidal height at 21.0012/117.4027
fid = fopen('Dongsha_2019_2020.txt'); %UTC time 
A = textscan(fid, '%s %s %s %s %s %s', 'Headerlines', 1);
d_tide = str2double(split(A{3}, '.'));
t_tide = str2double(split(A{4},':'));
d_time = datenum(d_tide(:,3), d_tide(:,1), d_tide(:,2), t_tide(:,1), t_tide(:,2), t_tide(:,3));
d_height = str2double(A{5}); %m
