% Amanda Syamsul
% April 11th 2024
% Comparing automated measurements to manual measurements

close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');

% setup filename file
% !ls -1 0*svg > Filenames
%%

% Read in CSV files
comparison_data = readtable('comparison_data.csv');
manual04 = readtable('manual0704.csv');
manual06 = readtable('manual0706.csv');
manual11 = readtable('manual0711.csv');
manual13 = readtable('manual0713.csv');

% Stack all manual data
all_manual_velocity = [manual04.v; manual06.v; manual11.v; manual13.v];

%%

figure()
hold on
plot(comparison_data.time, all_manual_velocity, 'Marker', '*', 'DisplayName', 'manual measurements');
errorbar(comparison_data.time, comparison_data.velocity, comparison_data.std_dev,'Marker', '*','DisplayName', 'automatic measurements');
xlabel('time (seconds)');
ylabel('velocity (m/s)');
title('Comparison between Manually & Automatically Measured Internal Wave Velocities');
disp(['diff. in avg. velocity = ', num2str(diff_avg), ' m/s'])
legend show;
grid on;
hold off

%% Save all comparison data to CSV file

% % Concatenate the arrays vertically
% all_t_curr_array = [t_curr_array4'; t_curr_array6'; t_curr_array11'; t_curr_array13'];
% all_velocity = [auto_velocity4'; auto_velocity6'; auto_velocity11'; auto_velocity13'];
% all_velocity_error = [velocity_error4'; velocity_error6'; velocity_error11'; velocity_error13'];
% 
% %%
% % Define column names
% col_names = {'time', 'velocity', 'std_dev'};
% 
% % Create the table with variable names directly
% T = table(all_t_curr_array, all_velocity, all_velocity_error, 'VariableNames', col_names);
% 
% % Write the table to a CSV file
% writetable(T, 'comparison_data.csv');

[auto_velocity4,velocity_error4,mean_ms,t_curr_array4] = calculate_velocity('0704_chosen', false, false);
auto_velocity(1) = NaN
compare('manual0704',auto_velocity4,t_curr_array4,velocity_error4)

%%
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
            plot(obs_1, '-')
            legend show;
        end

        % Cropping the longer vector to match wave lengths-==========-8
        
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
        
        if verbose
            % Visualization for verification
            figure;
            plot(x_prev, y_prev, 'b-', 'DisplayName', 'Previous Wave');
            hold on;
            plot(x_curr, y_curr, 'g-', 'DisplayName', 'Current Wave');
            plot(cropped_long_x, cropped_long_y, 'r--', 'DisplayName', sprintf('Adjusted %s Wave', chosen_name));
            legend show;
            xlabel('X Axis');
            ylabel('Y Axis');
            title(sprintf('Comparison of Adjusted and Original Vectors'));
        end

        % Compare y-lengths and replace longer wave vector with cropped version
        if yLengthCurr > yLengthPrev
            y_curr = cropped_long_y;
            x_curr = cropped_long_x;
        else
            y_prev = cropped_long_y;
            x_prev = cropped_long_x;
        end
        

        %% Use to only calculate near OBS
        
        if only_obs
            %% Filter the arrays to keep only values between 350 and 450
           
            obs_filter_prev = x_prev >= 350 & x_prev <= 450; 
            obs_filter_curr = x_curr >= 350 & x_curr <= 450;

    
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
        if length(x_prev) < 0.2*length(x_curr)
            continue
        end

        if length(x_curr) < 0.2*length(x_prev)
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
            if (calc_points(i,1) + gap)>length(x_prev)
                gap=1; % this has to be fixed!!!!
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
        
        % m* is the slope of the line connecting both waves
        mean_ms(f+1) = mean(ms);
        time(f+1) = curr_num;
        avg_v = mean(v);
        velocity(f+1) = avg_v;
        velocity_error(f+1)=std(v);

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

                plot(x_prev(calc_points(p, 1)), y_prev(calc_points(p,1)), 'bo')
                plot(x_curr(calc_points(p, 2)), y_curr(calc_points(p,2)), 'ko')
                plot(xs(p),ys(p),'m*'); leg4 = "perp. intersect";
                plot(x1(p),y1(p),'b*'); leg4 = "";
                plot(x2(p),y2(p),'k*'); leg4 = "perp. intersect";
                y_perp(:,p) = -(1/m1(p)) * (xvals - x1(p)) + y1(p);
                plot(xvals,y_perp(:,p), 'r--'); leg3 = "perp. line";
            
            end
        end

        if verbose
            disp(['Average velocity between SVG files: ' num2str(avg_v) 'm/s'])
            disp(['next!'])
        end
    end
end

function compare(manual_velocity,auto_velocity,time,auto_velocity_error)
    manual = readtable(manual_velocity);
    lim = height(manual);
    avg_manual = mean(manual.v(2:end));
    avg_auto = mean(auto_velocity(2:end));
    diff_avg = abs(avg_manual-avg_auto);

    figure(20)
    hold on
    plot(time',manual.v, 'Marker', '*', 'DisplayName', 'manual measurements');
    errorbar(time, auto_velocity, auto_velocity_error,'Marker', '*','DisplayName', 'automatic measurements');
    xlabel('time (seconds)');
    ylabel('velocity (m/s)');
    title('Comparison between Manually & Automatically Measured Internal Wave Velocities');
    disp(['diff. in avg. velocity = ', num2str(diff_avg), ' m/s'])
    legend show;
    grid on;
    hold off
end