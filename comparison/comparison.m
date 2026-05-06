% Amanda Syamsul
% April 11th 2024
% Fig. S5 - Comparing automated measurements to manual measurements

close all; clear all; clc

project_root = '/Users/amandasyamsul/Documents/MATLAB/OIW';
cd(project_root);

% Add project folders
addpath(fullfile(project_root, 'code'));
addpath(fullfile(project_root, 'comparison'));

comparison_data = readtable('comparison_data.csv');
manual04 = readtable('manual0704.csv');
manual06 = readtable('manual0706.csv');
manual11 = readtable('manual0711.csv');
manual13 = readtable('manual0713.csv');

% Stack all manual data
all_manual_velocity = [manual04.v; manual06.v; manual11.v; manual13.v];
%%

close all;

figure()

subplot(4,1,1)
plot(comparison_data.time, all_manual_velocity, 'Marker', 'o', 'MarkerFaceColor','blue', 'DisplayName', 'manual measurements');
hold on
errorbar(comparison_data.time, comparison_data.velocity, comparison_data.std_dev,'o', 'MarkerFaceColor','red','DisplayName', 'automatic measurements');

xlim([datetime(2020, 7, 4, 3, 0, 0), datetime(2020, 7, 4, 9, 0, 0)])
legend show;
ylim([0,5])
grid on;

subplot(4,1,2)
plot(comparison_data.time, all_manual_velocity, 'Marker', 'o', 'MarkerFaceColor','blue')
hold on
errorbar(comparison_data.time, comparison_data.velocity, comparison_data.std_dev,'o', 'MarkerFaceColor','red')
xlim([datetime(2020, 7, 6, 3, 0, 0), datetime(2020, 7, 6, 9, 0, 0)])
ylim([0,5])
grid on;

subplot(4,1,3)
plot(comparison_data.time, all_manual_velocity, 'Marker', 'o', 'MarkerFaceColor','blue')
hold on
errorbar(comparison_data.time, comparison_data.velocity, comparison_data.std_dev,'o', 'MarkerFaceColor','red')
xlim([datetime(2020, 7, 11, 3, 0, 0), datetime(2020, 7, 11, 9, 0, 0)])
ylim([0,5])
grid on;

subplot(4,1,4)
plot(comparison_data.time, all_manual_velocity, 'Marker', 'o', 'MarkerFaceColor','blue')
hold on
errorbar(comparison_data.time, comparison_data.velocity, comparison_data.std_dev,'o', 'MarkerFaceColor','red')
xlim([datetime(2020, 7, 13, 3, 0, 0), datetime(2020, 7, 13, 9, 0, 0)])
ylim([0,5])
grid on;

hold off
