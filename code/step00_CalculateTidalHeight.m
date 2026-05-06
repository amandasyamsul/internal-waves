% Amanda Syamsul
% Calculating tmd_predict output
% Relevant outputs: t_otps, z_otps (used in step02_HimawariOBPAnalysis.m)

% NOTE: TPXO9_atlas_v5.nc is needed for this script to run
close all; clear all; clc

project_root = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(project_root);

% Add project folders
addpath(fullfile(project_root, 'code'));
addpath(fullfile(project_root, 'data'));
addpath(fullfile(project_root, 'OTPS'));
addpath(fullfile(project_root, 'OTPS', 'chadagreene-Tide-Model-Driver-fd89bd6'));

all_detections = readtable('all_detections_new.csv');
obs_coords = [21.00116 117.40267]; % 21N 117E 
tide_model = 'TPXO9_atlas_v5.nc';
lat = obs_coords(1);
lon = obs_coords(2);

% OTPS (unchanged)
t_otps = (all_detections.DetectionTime(1):hours(1):all_detections.DetectionTime(end));
z_otps = tmd_predict(tide_model, lat, lon, t_otps, 'h');

%% Save outputs to file
col_names = {'t_otps', 'z_otps'};
tmd_output = table(t_otps', z_otps', 'VariableNames', col_names);

% writetable(tmd_output, 'tmd_output.csv');