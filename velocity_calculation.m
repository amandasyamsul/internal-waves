% Amanda Syamsul
% April 11th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll

close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');

% setup filename file
% !ls -1 0*svg > Filenames

%% Read in CSV files

luzon = readtable('luzon.csv');
% 
% april = readtable('aprildata.csv');
% 
may = readtable('maydata.csv');
may_away = readtable('may_away.csv');
may_near = readtable('may_near.csv');
% may_10km_60 = readtable('may_10km_60.csv');
% may_5km_60 = readtable('may_5km_60.csv');
% may_10km_60c = readtable('may_10km_60c.csv');
% may_5km_60c = readtable('may_5km_60c.csv');
% 
% june = readtable('junedata.csv');
% 
% july = readtable('julydata.csv');
% 
% august = readtable('augdata.csv');
% 
% sept = readtable('septdata.csv');

% [l_time, l_vel, d_time, d_height] = tide_2019_2020();
% 
% l_dates = datetime(l_time, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% d_dates = datetime(d_time, 'ConvertFrom', 'datenum', 'Format', 'd-MMM-y HH:mm:ss');
% %%
% l_vel = l_vel ./ 100; % Convert cm/s to m/s
% luzon = table(l_dates, l_vel);
% dongsha = table(d_dates, d_height);
%
% writetable(luzon,'luzon.csv')
% writetable(dongsha, 'dongsha.csv')

%%

[velocity,velocity_std,back_azimuth,back_azimuth_std,t_curr_array]=calculate_velocity('May2020', false, true, 10000, 0.6);

% Saving data to table & SVG 

% Define column names
col_names = {'time', 'velocity', 'velocity_std', 'back_azimuth', 'back_azimuth_std'};
T = table(t_curr_array', velocity', velocity_std', back_azimuth', back_azimuth_std', 'VariableNames', col_names);

writetable(T,'may_near.csv')

%% Plots!
close all; 
time_x = [datetime(2020, 5, 1, 1, 0, 0), datetime(2020, 5, 31, 7, 0, 0)];

%%

figure(1)
% Plot 2: Tidal velocity
subplot(2, 1, 1); % Second subplot in a 2x2 grid
plot(luzon.l_dates, luzon.l_vel)
xlim(time_x - hours(44.5));
title('Tidal velocities in the Luzon Strait in April 2020 (44.5 hour lag)');
xlabel('time');
ylabel('velocity (m/s)');

% Plot 1: Internal wave velocity
subplot(2, 1, 2); % First subplot in a 2x2 grid
% plot_avg_velocity(may_obs_10km, '10 km E-W of OBS — cropped waves out of bounds')
plot_avg_velocity(april, 'velocity')
% plot_avg_velocity(may_10km_60c, '10 km E-W of OBS - waves 60% in bounds above Dongsha')
xlim(time_x);
title('Internal wave velocities in April 2020');

%%

close all;

% Plot velocities with different constraints
figure()
hold on;
plot_velocity(may, 'total distance')
plot_velocity(may_away, 'waves > 42 km E away from atoll')
plot_velocity(may_near, 'waves < 42 km E from atoll')
xlim(time_x);
title('Internal wave velocities in May 2020');

%%
close all; 

figure(2);
hold on
plot_azimuth(sept, 'total distance')

figure(3);
plot_azimuth2(sept)

%% Daily velocity plots

% maydays = [5, 6, 7, 9, 10, 12, 26, 27, 29];
% junedays = [4 11 19 20 21 22 24 25 27 28];
% julydays = [6 11 12 13 27];

julydays = [13 27];

figure(4)

for p = 1:length(julydays)
    time_x = [datetime(2020, 7, julydays(p), 1, 0, 0), datetime(2020, 7, julydays(p), 11, 0, 0)];
    
    % Plot 1: Tidal velocity
    subplot(length(julydays)*2, 1, 2*p-1); % Odd index for tidal velocity
    plot(luzon.l_dates, luzon.l_vel)
    xlim(time_x - hours(44.5));
    % title(['Tidal velocities in the Luzon Strait in May 2020 (44.5 hour lag) - Day ', num2str(maydays(p))]);
    % xlabel('time');
    % ylabel('velocity (m/s)');
    
    % Plot 2: Internal wave velocity
    subplot(length(julydays)*2, 1, 2*p); % Even index for internal wave velocity
    plot_velocity(july_obs_5km)
    xlim(time_x);
    % title(['Internal wave velocities in May 2020 - Day ', num2str(maydays(p))]);
end

%% Sensitivity tests

data = may_5km_60;

data.DateOnly = dateshift(data.time, 'start', 'day');

% Group by date and calculate mean velocity and standard deviation
dailyStats = groupsummary(data, 'DateOnly', {'mean', 'std'}, 'velocity');

% Extract necessary columns
avgTime3 = dailyStats.DateOnly;
avgVelocity3 = dailyStats.mean_velocity;

%%
clc
close all

data1 = avgVelocity1;
data2 = avgVelocity3;

% Find common datetime values
[common_times, idx1, idx2] = intersect(avgTime1, avgTime3);

% Extract matched data based on common times
matched_data1 = data1(idx1);
matched_data2 = data2(idx2);

% F-test for variance
[h_var, p_var] = vartest2(matched_data1, matched_data2);

% t-test for means
[h_mean, p_mean] = ttest2(matched_data1, matched_data2);

% Bland-Altman Plot
mean_data = (matched_data1 + matched_data2) / 2;
diff_data = matched_data1 - matched_data2;
mean_diff = mean(diff_data);
std_diff = std(diff_data);

figure;
scatter(mean_data, diff_data, 'o', 'Filled');
xlabel('Mean of Two Methods');
ylabel('Difference of Two Methods');
title('Bland-Altman Plot');
hold on;
yline(mean_diff, 'r-', 'LineWidth', 2, 'DisplayName', 'Mean Difference');
yline(mean_diff + 1.96 * std_diff, 'r--', 'LineWidth', 2, 'DisplayName', 'Mean + 1.96 SD');
yline(mean_diff - 1.96 * std_diff, 'r--', 'LineWidth', 2, 'DisplayName', 'Mean - 1.96 SD');
legend('Data Points', 'Mean Difference', 'Mean + 1.96 SD', 'Mean - 1.96 SD');
grid on;
hold off;

% Display results
fprintf('P-value from F-test comparing variances: %f\n', p_var);
fprintf('P-value from t-test comparing means: %f\n', p_mean);

%%
function plot_azimuth2(data)
    % Convert back azimuth to radians for polar plot
    azimuths_rad = deg2rad(data.back_azimuth);
    
    % Scale normalized time to desired radius range, e.g., [0, 2]
    radius = data.velocity;
    
    p = polaraxes; % Create a polar axes context
    hold on;

    % Adjust the polar plot to set 0 degrees to North and clockwise direction
    p.ThetaDir = 'clockwise';
    p.ThetaZeroLocation = 'top';

    % Create polar scatter plot with velocity as radius using scatter
    scatter(p, azimuths_rad, radius, 30, 'k', 'filled');
    
    % Add a custom legend to indicate that the radius represents velocity
    h = legend('Back azimuth',  'Location', 'west', 'Orientation', 'horizontal');
    title(h, 'Radius = Velocity (m/s)');

    % Enhance the plot
    title('Direction of Wave Source');
    legend show;
    hold off;
end


function plot_azimuth(data, label)
    % Create a scatter plot with azimuths
    
    % Convert datetime to date only (removing time part)
    data.DateOnly = dateshift(data.time, 'start', 'day');

    % Group by date and calculate mean velocity and standard deviation
    dailyStats = groupsummary(data, 'DateOnly', {'mean', 'std'}, 'back_azimuth');
    
    % Extract necessary columns
    avgTime = dailyStats.DateOnly;
    avg_back_azimuth = dailyStats.mean_back_azimuth;
    std_back_azimuth = dailyStats.std_back_azimuth;
   
    hold on
    % Plot mean back azimuth with error bars
    errorbar(avgTime, avg_back_azimuth, std_back_azimuth, 'r-', 'LineWidth', 1, 'DisplayName', label);
    
    % Plot individual back azimuth values
    scatter(data.time, data.back_azimuth, 20, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
    
    grid on;
    
    % Label the axes
    title('Direction of Wave Source');
    xlabel('Time');
    % ylabel('Back Azimuth (degrees)');
    
    % Show legend
    % legend show;
    
    hold off;
end

function plot_avg_velocity(data, caption)

    % Convert datetime to date only (removing time part)
    data.DateOnly = dateshift(data.time, 'start', 'day');

    % Group by date and calculate mean velocity and standard deviation
    dailyStats = groupsummary(data, 'DateOnly', {'mean', 'std'}, 'velocity');
    
    % Extract necessary columns
    avgTime = dailyStats.DateOnly;
    avgVelocity = dailyStats.mean_velocity;
    stdDevVelocity = dailyStats.std_velocity;

    % Plotting
    errorbar(avgTime, avgVelocity, stdDevVelocity, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', caption);
    hold on;
    % plot(avgTime, avgVelocity, 'o-', 'LineWidth', 1, 'DisplayName', caption, 'MarkerFaceColor', 'blue');
    % scatter(avgTime, avgVelocity, 'ko', 'Filled','DisplayName', 'average velocity');
    xlabel('Time');
    ylabel('Velocity (m/s)');
    legend show;
    grid on;
    set(gca, 'FontSize', 12);
end

function plot_velocity(data, caption)
    % Plot velocity vs time with error bars
    errorbar(data.time, data.velocity, data.velocity_std, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'error bars');
    hold on
    plot(data.time, data.velocity,'-','LineWidth', 1,'DisplayName', caption);
    xlabel('Time');
    ylabel('Velocity (m/s)');
    legend show;
    grid on;
    set(gca, 'FontSize', 12);
end

function [velocity,velocity_std,back_azimuth,back_azimuth_std,t_curr_array]=calculate_velocity(data_file, verbose, only_obs, constraint, prop)

    data = readtable(data_file, 'ReadVariableNames', false); 
    
    data_array = table2array(data);

    for f=1:length(data_array)-1
        file1 = data_array{f};
        file2 = data_array{f+1};
    
        prev = loadsvg(file1,1,1);
        curr = loadsvg(file2,1,1);
        obs = loadsvg('OBS.svg',1,1);
    
        % Failsafe in case there is a smaller vector in the svg that should
        % be ignored
        if (length(prev)>1)
            for icheck=1:length(prev)
                Ncheck(icheck)=length(prev{icheck});
            end
            [~,Ilong]=max(Ncheck);
            prev=prev(Ilong); % pick the longest vector
        end
    
        if (length(curr)>1)
            for icheck=1:length(curr)
                Ncheck(icheck)=length(curr{icheck});
            end
            [~,Ilong]=max(Ncheck);
            curr=curr(Ilong); % pick the longest vector
        end
    
        close all; 
        
        dt_steps = readmatrix('dt-picking.csv');
        
        % Extract last three digits from file names
        prev_num = str2double(file1(end-14:end-11));
        curr_num = str2double(file2(end-14:end-11));

        % Find the corresponding dt values from the CSV data
        prev_dt = find(dt_steps(:,1) == prev_num);
        curr_dt = find(dt_steps(:,1) == curr_num);
        dt = (curr_dt - prev_dt) * 600;

        % Check if dt is less than 0 (between two diff. days), and if so, skip to the next iteration
        if dt < 0
            continue;
        end
        
        if verbose
          disp(['Time difference (dt) between SVG files: ' num2str(dt) ' seconds']);
        end
      
        % Set datetime format for previous wave
        % Extract month, day, and time parts from the filename
        year_prev = str2double(file1(1:4));
        month_prev = str2double(file1(5:6));
        day_prev = str2double(file1(7:8));
        time_prev = file1(9:12);
        
        % Determine hour and minute from the time part
        if length(time_prev) == 4
            hour_prev = str2double(time_prev(1:2));  % First two digits are the hour
            min_prev = str2double(time_prev(3:4)); % Last two digits are the minutes
        else
            hour_prev = str2double(time_prev(1:2)); % First two digits are the hour
            min_prev = str2double(time_prev(3:4)); % Last two digits are the minutes
        end
        
        % Construct the datetime object assuming time is in AM UTC
        t_prev = datetime(year_prev, month_prev, day_prev, hour_prev, min_prev, 0, 'TimeZone', 'UTC');

        % Set datetime format for previous wave
        % Extract month, day, and time parts from the filename
        year_curr = str2double(file2(1:4));
        month_curr = str2double(file2(5:6));
        day_curr = str2double(file2(7:8));
        time_curr = file2(9:12);
        
        % Determine hour and minute from the time part
        if length(time_curr) == 4
            hour_curr = str2double(time_curr(1:2));  % First two digits are the hour
            min_curr = str2double(time_curr(3:4)); % Last two digits are the minutes
        else
            hour_curr = str2double(time_curr(1:2)); % First two digits are the hour
            min_curr = str2double(time_curr(3:4)); % Last two digits are the minutes
        end
        
        % Construct the datetime object assuming time is in AM UTC
        t_curr = datetime(year_curr, month_curr, day_curr, hour_curr, min_curr, 0, 'TimeZone', 'UTC');        
        
        % Define latitude range covered by the image in meters
        % Using scale bar in satellite images, 25 km = 70.03 pixels
        meters_per_pix = 25000/70.03;
        
        % close all; clc
        
        original_x_prev=prev{1,1}(:,1);
        original_y_prev=prev{1,1}(:,2);
        original_x_curr=curr{1,1}(:,1);
        original_y_curr=curr{1,1}(:,2);


        %% Interpolate data
    
        % Previous wave:
        % Combine x and y into a matrix and sort rows by the first column (x)
        data = [original_x_prev, original_y_prev];
        data = sortrows(data, 1);
    
        % Find unique x values and average corresponding y values
        [y_unique, ia, ic] = unique(data(:,2));
        x_avg = accumarray(ic, data(:,1), [], @mean);
    
        % Interpolation query points
        yq = linspace(min(y_unique), max(y_unique), numel(y_unique)*10);  % A more reliable query range
    
        % Interpolate
        vq1 = interp1(y_unique, x_avg, yq);
    
        % Combine and sort x values
        x_combined = [data(:,1); vq1'];
        y_combined = [data(:,2); yq'];
    
        % Sort the combined arrays by x values
        [y_sorted, sort_index] = sort(y_combined);
        x_sorted = x_combined(sort_index);
    
        x_prev = x_sorted;
        y_prev = y_sorted;
    
        % Curent wave:
    
        % Combine x and y into a matrix and sort rows by the first column (x)
        data = [original_x_curr, original_y_curr];
        data = sortrows(data, 1);
    
        % Find unique x values and average corresponding y values
        [y_unique, ia, ic] = unique(data(:,2));
        x_avg = accumarray(ic, data(:,1), [], @mean);
    
        % Interpolation query points
        yq = linspace(min(y_unique), max(y_unique), numel(y_unique)*10);  % A more reliable query range
    
        % Interpolate
        vq1 = interp1(y_unique, x_avg, yq);
    
        % Combine and sort x values
        x_combined = [data(:,1); vq1'];
        y_combined = [data(:,2); yq'];
    
        % Sort the combined arrays by x values
        [y_sorted, sort_index] = sort(y_combined);
        x_sorted = x_combined(sort_index);
    
        x_curr = x_sorted;
        y_curr = y_sorted;

        % Cropping the longer vector to match wave lengths
        
        % Calculate y-lengths (ranges)
        yLengthPrev = max(y_prev) - min(y_prev);
        yLengthCurr = max(y_curr) - min(y_curr);
        
        % Compare y-lengths and choose the vector with the shorter y
        if yLengthCurr > yLengthPrev
            short_x = x_prev;
            short_y = y_prev;
            chosen_name = 'Previous';
            long_x = x_curr;
            long_y = y_curr;
        else
            short_x = x_curr;
            short_y = y_curr;
            chosen_name = 'Previous';
            long_x = x_prev;
            long_y = y_prev;
        end
        
        % Find the minimum and maximum y-values of the shorter vector
        ymin = min(short_y);
        ymax = max(short_y);
        
        % Find the indices in the longer vector where the y-values are within the range of the shorter vector
        indices = long_y >= ymin & long_y <= ymax;
        
        % Crop the longer vector to match the range and y-max/y-min of the shorter vector
        cropped_long_y = long_y(indices);
        cropped_long_x = long_x(indices);

        % Compare y-lengths and replace longer wave vector with cropped version
        if yLengthCurr > yLengthPrev
            y_curr = cropped_long_y;
            x_curr = cropped_long_x;
        else
            y_prev = cropped_long_y;
            x_prev = cropped_long_x;
        end
                
        if verbose
            % Visualization for verification
            figure;
            plot(x_prev, y_prev, 'b-', 'DisplayName', 'Previous Wave');
            hold on;
            plot(x_curr, y_curr, 'g-', 'DisplayName', 'Current Wave');
            plot(cropped_long_x, cropped_long_y, 'r--', 'DisplayName', sprintf('Adjusted %s Wave', chosen_name));
            plot(392,159, 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);
            legend show;
            xlabel('X Axis');
            ylabel('Y Axis');
            title(sprintf('Comparison of Adjusted and Original Vectors'));
        end


        %% Use to only calculate near OBS

        if only_obs
            %% Define the boundary based on constraint and meters_per_pix
            bound = constraint / meters_per_pix;

            %% Check the proportion of x_prev and x_curr within the bounds of OBS
            % in_bounds_prev = (x_prev >= (392-bound)) & (x_prev <= (392+bound));
            % in_bounds_curr = (x_curr >= (392-bound)) & (x_curr <= (392+bound));

            % Away from atoll
            in_bounds_prev = (x_prev >= 258) & (x_prev <= 377);
            in_bounds_curr = (x_curr >= 258) & (x_curr <= 377);

            % Calculate the proportion of points within the bounds
            proportion_in_bounds_prev = sum(in_bounds_prev) / length(x_prev);
            proportion_in_bounds_curr = sum(in_bounds_curr) / length(x_curr);

            % Check if at least prop% of the points are within the bounds for both waves
            if proportion_in_bounds_prev < prop || proportion_in_bounds_curr < prop
                continue; % Skip iteration if less than prop% of the points are within the bounds
            end

            %% Comment out code below if using entire y-axis
            % % Filter the arrays to keep only values above Dongsha (y=251)
            % obs_filter_prev = y_prev <= 251; % less than because points are flipped
            % obs_filter_curr = y_curr <= 251;
            % 
            % % Apply the filters
            % x_prev = x_prev(obs_filter_prev);
            % y_prev = y_prev(obs_filter_prev);
            % x_curr = x_curr(obs_filter_curr);
            % y_curr = y_curr(obs_filter_curr);

            % % Skip iteration if any of the arrays are empty after filtering
            % if isempty(x_prev) || isempty(y_prev) || isempty(x_curr) || isempty(y_curr)
            %     continue;
            % end
            %%

            if verbose
                figure;
                plot(x_prev, y_prev, '-*', 'DisplayName', 'Previous Wave');
                hold on;
                plot(x_curr, y_curr, '-*', 'DisplayName', 'Current Wave');
                plot(392, 159, 'p', 'DisplayName', 'OBS Station', 'MarkerFaceColor', 'red', 'MarkerSize', 15);
                plot(199, 251, 'p', 'DisplayName', 'Dongsha', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);
                legend show;
            end
        end

        % Failsafe for OBS -- skip if one wave is too short
        if length(x_prev) < 0.3*length(x_curr)
            continue
        end

        if length(x_curr) < 0.3*length(x_prev)
            continue
        end

        %% Finding measurement points:
    
        % Define the fractions of the line where you need the indices
        fractions = [1/10, 2/10, 3/10, 4/10, 1/2, 6/10, 7/10, 8/10, 9/10];
    
        % Preallocate the matrix to store results
        calc_points = zeros(length(fractions), 2);

        % Loop over each fraction
        for i = 1:length(fractions)
            % Calculate the indices for the current and previous lines
            index_prev = round(length(y_prev) * fractions(i));
            index_curr = round(length(y_curr) * fractions(i));

            
            % Test code to see if this solves issue of rounding to 0
            if index_prev==0
                continue
            end

            if index_curr==0
                continue
            end
    
            % Store the results in the matrix
            calc_points(i, :) = [index_prev index_curr];
        end
            
        for i=1:length(fractions)
            % Calculate the slope at varying points on the previous wave
            gap=10;
    
            % edit gap for waves that have too few points!
            if (calc_points(i,1) - gap)<1
                gap=round(calc_points(i,1)/2);
            end
    
            if (calc_points(i,1) - gap)==0
                gap=0;
            end
            
            % Failsafe if vector is still too short -- skip to next iteration
            if (calc_points(i,1) - gap)==0
                continue
            end

            x1a = x_prev(calc_points(i,1) - gap);
            y1a = y_prev(calc_points(i,1) - gap);
            
            gap = 10;
            if (calc_points(i,1) + gap) > length(x_prev)
                gap = 1; % this has to be fixed!!!!
                if (calc_points(i,1) + gap) > length(x_prev)
                    continue;
                end
            end

            x2a = x_prev(calc_points(i,1) + gap);
            y2a = y_prev(calc_points(i,1) + gap);
            slope_prev(i) = (y2a - y1a) / (x2a - x1a);
            
            % Calculate the slope at varying points on the current wave
            gap=10;
            if (calc_points(i,2) - gap)<1
                gap=round(calc_points(i,2)/2);
            end
    
            if (calc_points(i,2) - gap)==0
                gap=0;
            end
    
            x1a= x_curr(calc_points(i,2) - gap);
            y1a= y_curr(calc_points(i,2) - gap);
            
            gap=10;

            if (calc_points(i,2) + gap)>length(x_curr)
                gap=1; % this has to be fixed!!!!
            end
            
            % Failsafe for OBS
            if (calc_points(i,2) + gap)>length(x_curr)
                continue
            end

            x2a = x_curr(calc_points(i,2) + gap);
            y2a = y_curr(calc_points(i,2) + gap);
            slope_curr(i) = (y2a - y1a) / (x2a - x1a);
        
            % Measure perp. dist from slope mid point to closest point on wave of next image
        
            % Define the coordinates of the points
            
            m1(i) = slope_prev(i);
            m2(i) = slope_curr(i);
        
            x1(i) = x_prev(calc_points(i,1));
            y1(i) = y_prev(calc_points(i,1));
        
            x2(i) = x_curr(calc_points(i,2));
            y2(i) = y_curr(calc_points(i,2));
    
            xs_num(i) = y2(i) - (m2(i)*x2(i)) -(x1(i)/m1(i)) -y1(i);
            xs_den(i) = -(1/m1(i)) -m2(i);
            xs(i) = xs_num(i)/xs_den(i);
            ys(i) = m2(i) * (xs(i) - x2(i)) + y2(i);
            ms(i) = (ys(i) - y1(i)) / (xs(i) - x1(i));
            
            xvals = 200:500;
            y_perp(:,i) = -(1/m1(i)) * (xvals - x1(i)) + y1(i);
    
            point1(i,:) = [x1(i) y1(i)];
            point2(i,:) = [xs(i) ys(i)];
        
            % Calculate velocity using all points
            if any(isnan(point1(i,:))) || any(isnan(point2(i,:)))
                continue
            end
                
            distance_pix(i) = norm(point1(i,:) - point2(i,:));
            distance_m(i) = distance_pix(i) * meters_per_pix;
            v(i) = distance_m(i)/dt;
        end
        
        try
            % m* is the slope of the line connecting both waves

            % Calculate the slope angle in degrees

            % Calculate the slope angle in degrees
            mean_ms = mean(ms);

            % Assuming mean_ms is a slope value, convert it to an angle in degrees
            ms_degree = atand(mean_ms);
            
            % Convert the slope angle to azimuth
            % Assume north is 0 degrees, east is 90 degrees, south is 180 degrees, west is 270 degrees
            if ms_degree >= 0
                back_az = 90 + ms_degree;
            else
                back_az = 90 - abs(ms_degree);
            end

            back_azimuth(f+1) = mean(back_az);
            back_azimuth_std(f+1) = std(back_az);
            velocity(f+1) = mean(v);
            velocity_std(f+1)=std(v);
        catch
            continue
        end

        % Store the datetime object
        t_curr_array(f+1) = t_curr;

        if verbose
            disp(t_curr);
        end
    
        %% logic testing
        if verbose
            figure(10)
            clf
            hold on
            
            plot(x_prev,y_prev,'b'); leg7 = "previous wave";
            plot(x_curr,y_curr,'k'); leg8 = "current wave";
            
            for p = 1:length(fractions)

                if calc_points(p, 1) == 0
                    continue
                end
                
                if length(fractions) > length(p)
                    continue
                end

                plot(x_prev(calc_points(p, 1)), y_prev(calc_points(p,1)), 'bo')
                plot(x_curr(calc_points(p, 2)), y_curr(calc_points(p,2)), 'ko')
                plot(xs(p),ys(p),'m*'); leg4 = "perp. intersect";
                plot(x1(p),y1(p),'b*'); leg4 = "";
                plot(x2(p),y2(p),'k*'); leg4 = "perp. intersect";
                y_perp(:,p) = -(1/m1(p)) * (xvals - x1(p)) + y1(p);
                plot(xvals,y_perp(:,p), 'r--'); leg3 = "perp. line";
                plot(392,159, 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);
                plot(199, 251, 'p', 'DisplayName', 'Dongsha', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);
            
            end
        end

        if verbose
            disp(['Average velocity between SVG files: ' num2str(velocity) 'm/s'])
            disp(['next!'])
        end
    end
end

function compare(manual_velocity,auto_velocity,auto_velocity_error)
    manual = readtable(manual_velocity);
    lim = height(manual);
    avg_manual = mean(manual.v(2:end));
    avg_auto = mean(auto_velocity(2:end));
    diff_avg = abs(avg_manual-avg_auto);

    figure(20)
    hold on
    plot(manual.seconds(1:lim),manual.v(1:lim), 'DisplayName', 'manual measurements');
    errorbar(manual.seconds(1:lim), auto_velocity(1:lim), auto_velocity_error(1:lim),'DisplayName', 'automatic measurements');
    xlabel('time (seconds)');
    ylabel('manually measured velocity (m/s)');
    title('comparison between manually and automatically measured velocities');
    disp(['diff. in avg. velocity = ', num2str(diff_avg), ' m/s'])
    legend show;
    hold off
end

function filtered_array = filter_continuous_ones(array)
    
    % Find the start and end indices of each contiguous segment of 1s
    d = diff([0; array; 0]);
    start_indices = find(d == 1);
    end_indices = find(d == -1) - 1;
    
    % Calculate the length of each segment
    segment_lengths = end_indices - start_indices + 1;
    
    % Find the index of the longest segment
    [~, longest_segment_index] = max(segment_lengths);
    
    % Get the start and end indices of the longest segment
    longest_segment_start = start_indices(longest_segment_index);
    longest_segment_end = end_indices(longest_segment_index);
    
    % Initialize the output array with zeros
    filtered_array = zeros(size(array));
    
    % Set the elements of the longest segment to 1
    filtered_array(longest_segment_start:longest_segment_end) = 1;

    filtered_array = logical(filtered_array);
end
