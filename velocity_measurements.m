% Amanda Syamsul
% April 11th 2024
% Analysis on OIW Position & Velocity Near Dongsha Atoll in July-Dec 2020
close all; clear all; clc

addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/loadsvg_2');
addpath('/Users/amandasyamsul/Documents/MATLAB/OIW/annotated');

% setup filename file
% !ls -1 0*svg > Filenames


%%
% [velocity_error, velocity] = calculate_velocity('0704_chosen', true)

% Initialize arrays to store velocities and velocity errors
velocity_error_array = [];
velocity_array = cell(1, 13);  % Use a cell array to store velocities of varying sizes

% Loop through the days
for day = 1:13
    % Construct the file name based on the day
    if day < 10
        file_name = sprintf('070%d_chosen', day);
    else
        file_name = sprintf('07%d_chosen', day);
    end
    
    % Construct the full path to the file in the 'annotated' directory
    full_file_path = fullfile('annotated', file_name);
    
    % Check if the file exists
    if exist(full_file_path, 'file')
        % Calculate velocity and store the results
        [velocity_error, velocity] = calculate_velocity(file_name, true);
        
        % Store the results in the arrays
        % velocity_error_array = [velocity_error_array; velocity_error];
        velocity_array{day} = velocity';  % Store each velocity as a cell
    else
        % Skip if the file does not exist
        fprintf('File %s does not exist. Skipping...\n', full_file_path);
    end
end

% Display or process the results as needed
disp('Velocity errors:');
disp(velocity_error_array);

disp('Velocities:');
disp(velocity_array);

%%

velocity(1) = NaN;

x_nums = 1:length(velocity);
figure()
hold on
errorbar(x_nums, velocity,velocity_error,'DisplayName', 'automatic measurements');
xlabel('x axis');
ylabel('velocity (m/s)');
title('automatically measured velocities');
legend show;
hold off

%% Comparing to manually measured values

% compare('dist0711.csv',velocity,velocity_error)

%% Functions

function [velocity_error,velocity]=calculate_velocity(data_file,print_or_not)

    verbose = print_or_not; % default is to not print any statements

    data = readtable(data_file, 'ReadVariableNames', false); % using if
    
    data_array = table2array(data);
    
    for f=1:length(data_array)-1
        file1 = data_array{f};
        file2 = data_array{f+1};
    
        prev = loadsvg(file1,1,1);
        curr = loadsvg(file2,1,1);
    
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
        prev_num = str2double(file1(end-6:end-4));
        curr_num = str2double(file2(end-6:end-4));
        
        % Find the corresponding dt values from the CSV data
        prev_dt = find(dt_steps(:,1) == prev_num);
        curr_dt = find(dt_steps(:,1) == curr_num);
        
        dt = (curr_dt - prev_dt) * 600;
        
        if verbose
          disp(['Time difference (dt) between SVG files: ' num2str(dt) ' seconds']);
        end

        % Define latitude range covered by the image in meters
        lat_degrees = 1.5;
        r_earth = 6378100; % meters
        lat_meters = (2*pi*r_earth*lat_degrees)/360;
        meters_per_pix=366.6227;
        
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
        %%
    
        figure;
        plot(x_prev, y_prev, '-*', 'DisplayName', 'Previous Wave');
        hold on;
        plot(x_curr, y_curr, '-*', 'DisplayName', 'Current Wave');
    
        %%
        
        % Cropping the longer vector
        
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
        
        % Compare y-lengths and replace longer wave vector with cropped version
        if yLengthCurr > yLengthPrev
            y_curr = cropped_long_y;
            x_curr = cropped_long_x;
        else
            y_prev = cropped_long_y;
            x_prev = cropped_long_x;
        end
        
        %% Finding measurement points:
    
        % Define the fractions of the line where you need the indices
        fractions = [1/10, 2/10, 3/10, 4/10, 1/2, 6/10, 7/10, 8/10, 9/10];
    
        % Preallocate the matrix to store results
        calc_points = zeros(length(fractions), 2);
    
        % Loop over each fraction
        for i = 1:length(fractions)
            % Calculate the indices for the current and previous lines
            index_curr = round(length(y_curr) * fractions(i));
            index_prev = round(length(y_prev) * fractions(i));
    
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
        
        avg_v = mean(v);
        velocity(f+1) = avg_v;
        velocity_error(f+1)=std(v);
    
        %% logic testing
    
        figure(10)
        clf
        hold on
        
        plot(x_prev,y_prev,'b'); leg7 = "previous wave";
        plot(x_curr,y_curr,'k'); leg8 = "current wave";
        
        for p = 1:length(fractions)
            plot(x_prev(calc_points(p, 1)), y_prev(calc_points(p,1)), 'bo')
            plot(x_curr(calc_points(p, 2)), y_curr(calc_points(p,2)), 'ko')
            plot(xs(p),ys(p),'m*'); leg4 = "perp. intersect";
            plot(x1(p),y1(p),'b*'); leg4 = "";
            plot(x2(p),y2(p),'k*'); leg4 = "perp. intersect";
            y_perp(:,p) = -(1/m1(p)) * (xvals - x1(p)) + y1(p);
            plot(xvals,y_perp(:,p), 'r--'); leg3 = "perp. line";
        
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