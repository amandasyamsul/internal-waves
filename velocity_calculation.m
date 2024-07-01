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
may = readtable('maydata.csv');
may_obs = readtable('may_obs.csv')
% june = readtable('junedata.csv');
% june_obs = readtable('june_obs.csv');
% july = readtable('julydata.csv');
% july_obs = readtable('july_obs.csv');
% august = readtable('augdata.csv');

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

% [velocity,velocity_error,mean_ms,t_curr_array] = calculate_velocity('May2020', true, true);

%%

% Calculate the slope angle in degrees
ms_degree = atand(mean_ms);
% Convert the slope angle to azimuth
% Assume north is 0 degrees, east is 90 degrees, etc.
if ms_degree >= 0
    back_azimuth = 90 - ms_degree;
    azimuth = back_azimuth + 180;
else
    back_azimuth = 90 + abs(ms_degree);
    azimuth = back_azimuth + 180;
end

% % Ensure azimuth is within 0 to 360 degrees
% if azimuth < 0
%     azimuth = azimuth + 360;
% elseif azimuth >= 360
%     azimuth = azimuth - 360;
% end

% Create a scatter plot with azimuths
figure(1);
scatter(t_curr_array, back_azimuth, 20, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
grid on;

% Label the axes
title('Direction of wave source')
xlabel('Time');
ylabel('Back azimuth (degrees)');

%% Saving data to table & SVG

% Define column names
col_names = {'time', 'velocity', 'std_dev', 'azimuth', 'back_azimuth'};
T = table(t_curr_array', velocity', velocity_error', azimuth', back_azimuth', 'VariableNames', col_names);

writetable(T,'may_obs.csv')

%% Plots!

figure(2)

time_x = [datetime(2020, 5, 7, 3, 0, 0), datetime(2020, 5, 7, 7, 0, 0)];

% Plot 1: Internal wave velocity
subplot(2, 1, 1); % First subplot in a 2x2 grid
plot_velocity(may_obs)
xlim(time_x);
title('Internal wave velocities in May 2020');

% Plot 2: Tidal velocity
% Create a logical index for timeframe of choice
luzon.l_dates = luzon.l_dates + hours(50); % timelag for waves to propagate from luzon strait
timeframe = month(luzon.l_dates) == 5;

subplot(2, 1, 2); % Second subplot in a 2x2 grid
plot(luzon.l_dates(timeframe), luzon.l_vel(timeframe))
xlim(time_x);
title('Tidal velocities in the Luzon Strait in May 2020 (50 hour lag)');
xlabel('time');
ylabel('velocity (m/s)');

%% Comparing to manually measured values

function plot_avg_velocity(data)

    % Convert datetime to date only (removing time part)
    data.DateOnly = dateshift(data.time, 'start', 'day');

    % Group by date and calculate mean velocity and standard deviation
    dailyStats = groupsummary(data, 'DateOnly', {'mean', 'std'}, 'velocity');
    
    % Extract necessary columns
    avgTime = dailyStats.DateOnly;
    avgVelocity = dailyStats.mean_velocity;
    stdDevVelocity = dailyStats.std_velocity;

    % Plotting
    errorbar(avgTime, avgVelocity, stdDevVelocity, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'error bars');
    hold on;
    plot(avgTime, avgVelocity, 'r-', 'LineWidth', 1, 'DisplayName', 'trend');
    scatter(avgTime, avgVelocity, 'ko', 'Filled','DisplayName', 'average velocity');
    xlabel('Time');
    ylabel('Velocity (m/s)');
    legend show;
    grid on;
    set(gca, 'FontSize', 12);
end

function plot_velocity(data_file)
    % Plot velocity vs time with error bars
    errorbar(data_file.time, data_file.velocity, data_file.std_dev, 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'error bars');
    hold on
    plot(data_file.time, data_file.velocity,'r-','LineWidth', 1,'DisplayName', 'velocity');
    xlabel('Time');
    ylabel('Velocity (m/s)');
    legend show;
    grid on;
    set(gca, 'FontSize', 12);
end

function [velocity,velocity_error,mean_ms,t_curr_array]=calculate_velocity(data_file, verbose, only_obs)

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

        if verbose
            figure;
            plot(x_prev, y_prev, '-*', 'DisplayName', 'Previous Wave');
            hold on;
            plot(x_curr, y_curr, '-*', 'DisplayName', 'Current Wave');
            plot(392,159, 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);            
            legend show;
        end

        if only_obs
            bound = 10000/meters_per_pix;
            yPrev_filter = y_prev >= (159-bound) & y_prev <= (159+bound);
            yCurr_filter = y_curr >= (159-bound) & y_curr <= (159+bound);

            x_prev = x_prev(yPrev_filter);
            y_prev = y_prev(yPrev_filter);
            x_curr = x_curr(yCurr_filter);
            y_curr = y_curr(yCurr_filter);

            % Skip iteration if any of the arrays are empty after filtering
            if isempty(y_prev) || isempty(y_curr)
                continue;
            end

            if verbose
                figure;
                plot(x_prev, y_prev, '-*', 'DisplayName', 'Previous Wave');
                hold on;
                plot(x_curr, y_curr, '-*', 'DisplayName', 'Current Wave');
                plot(392,159, 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);
                legend show;
            end

        end

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
            %% Filter the arrays to keep only values within 10 km radius from OBS
            
            bound = 10000/meters_per_pix;

            obs_filter_prev = x_prev >= (392-bound) & x_prev <= (392+bound);
            obs_filter_curr = x_curr >= (392-bound) & x_curr <= (392+bound);
    
            % Apply the filters
            x_prev = x_prev(obs_filter_prev);
            y_prev = y_prev(obs_filter_prev);
            x_curr = x_curr(obs_filter_curr);
            y_curr = y_curr(obs_filter_curr);
    
            % Skip iteration if any of the arrays are empty after filtering
            if isempty(x_prev) || isempty(y_prev) || isempty(x_curr) || isempty(y_curr)
                continue;
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
            mean_ms(f+1) = mean(ms);
            time(f+1) = curr_num;
            avg_v = mean(v);
            velocity(f+1) = avg_v;
            velocity_error(f+1)=std(v);
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
            
            end
        end

        if verbose
            disp(['Average velocity between SVG files: ' num2str(avg_v) 'm/s'])
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

