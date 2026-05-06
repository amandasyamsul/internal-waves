function [velocity,velocity_std_dev,velocity_std_err,back_azimuth,t_curr_array, mid_x]=calculate_velocity(data_file, verbose, constrained, prop)

    data = readtable(data_file, 'ReadVariableNames', false); 
    
    data_array = table2array(data);

    for f=1:length(data_array)-1
        
        file1 = data_array{f};
        file2 = data_array{f+1};
    
        try
            prev = loadsvg(file1,1,0);
        catch
            continue; % Skip to the next iteration if file1 causes an error
        end
        
        try
            curr = loadsvg(file2,1,0);
        catch   
            continue; % Skip to the next iteration if file2 causes an error
        end

        obs = loadsvg('OBS.svg',1,0);
    
        % Failsafe in case there is a smaller vector in the svg that should be ignored
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
        
        obs_x1 = 392.5567;
        obs_x2 = 425.8900;
        obs_y1 = 133.0240;
        obs_y2 = 106.3573;

        % Square boundary for OBS velocity measurements
        obs_x = [obs_x1, obs_x2, obs_x2, obs_x1, obs_x1];
        obs_y = [obs_y1, obs_y1, obs_y2, obs_y2, obs_y1];

        obs_coords_im = [400.8900 133.0240];

        if verbose
            % Visualization for verification
            figure;
            plot(x_prev, y_prev, 'b-', 'DisplayName', 'Previous Wave');
            hold on;
            plot(x_curr, y_curr, 'k-', 'DisplayName', 'Current Wave');
            plot(cropped_long_x, cropped_long_y, 'r--', 'DisplayName', sprintf('Adjusted %s Wave', chosen_name));
            plot(obs_coords_im(1),obs_coords_im(2), 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','blue','Marker', '^','MarkerSize',10);
            plot(obs_x, obs_y, 'k-', 'LineWidth', 1, 'DisplayName', 'OBS bounds');
            xlim([300 450])
            ylim([0 250])
            legend show;
            axis ij
            grid on
            title(sprintf('Comparison of Adjusted and Original Vectors'));
        end


        %% Use if statement to calculate with spatial constraints

        if constrained
            %% Define the boundary based on constraint and meters_per_pix

            % Cut off waves before refracting around atoll
            in_bounds_prev = (x_prev >= 258);
            in_bounds_curr = (x_curr >= 258);

            % Calculate the proportion of points within the bounds
            proportion_in_bounds_prev = sum(in_bounds_prev) / length(x_prev);
            proportion_in_bounds_curr = sum(in_bounds_curr) / length(x_curr);

            % Check if at least prop of the points are within the bounds for both waves
            if proportion_in_bounds_prev < prop || proportion_in_bounds_curr < prop
                continue; % Skip iteration if less than prop of the points are within the bounds
            end

            % all waves:                            COMMENT OUT code block below 
            % OBS-constrained (square around OBS):  keep all code
            % latband across OBS:                   don't include x constraints
            obs_filter_prev = (y_prev <= obs_y1) & (y_prev >= obs_y2); %& (x_prev >= obs_x1) & (x_prev <= obs_x2);
            obs_filter_curr = (y_curr <= obs_y1) & (y_curr >= obs_y2); %& (x_curr >= obs_x1) & (x_curr <= obs_x2);
            % Apply the filters
            x_prev = x_prev(obs_filter_prev);
            y_prev = y_prev(obs_filter_prev);
            x_curr = x_curr(obs_filter_curr);
            y_curr = y_curr(obs_filter_curr);
            

            % Skip iteration if any of the arrays are empty after filtering
            if isempty(x_prev) || isempty(y_prev) || isempty(x_curr) || isempty(y_curr)
                continue;
            end

            if verbose
                figure;
                plot(x_prev, y_prev, '-*', 'DisplayName', 'Previous Wave');
                hold on;
                plot(x_curr, y_curr, '-*', 'DisplayName', 'Current Wave');
                plot(obs_coords_im(1),obs_coords_im(2), 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','red','MarkerSize',15);
                plot(obs_x, obs_y, 'k-', 'LineWidth', 1, 'DisplayName', 'OBS bounds');
                xlim([300 500])
                ylim([0 400])
                grid on
                legend show;
                axis ij
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
        num_of_points = 10;
        fractions = [1:num_of_points]./num_of_points;
    
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
            
            % % Failsafe if vector is still too short -- skip to next iteration
            % if (calc_points(i,1) - gap)==0
            %     continue
            % end

            x1a = x_prev(calc_points(i,1) - gap);
            y1a = y_prev(calc_points(i,1) - gap);
            
            gap = 10;
            if (calc_points(i,1) + gap) > length(x_prev)
                gap = 1; 
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
            % m* (ms) is the slope of the line connecting both waves
    
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
            velocity(f+1) = mean(v);
            velocity_std_dev(f+1)=std(v);
    
            % standard error
            n = length(v);
            velocity_std_err(f+1) = velocity_std_dev(f+1) / sqrt(n);
    
            % retain the wave's center x-coords (to be converted to longitude)
            mid_ind = round(length(x_curr)/2);
            mid_x(f+1) = x_curr(mid_ind);

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
            hold on
            
            plot(x_prev,y_prev,'r'); leg7 = "previous wave";
            plot(x_curr,y_curr,'k'); leg8 = "current wave";
            
            for p = 1:length(xs)

                if calc_points(p, 1) == 0
                    continue
                end
                
                % if length(fractions) > length(p)
                %     continue
                % end

                % plot(xs(p),ys(p),'m*'); leg4 = "perp. intersect";
                % plot(x1(p),y1(p),'b*'); leg4 = "";
                % plot(x2(p),y2(p),'k*'); leg4 = "perp. intersect";
                y_perp(:,p) = -(1/m1(p)) * (xvals - x1(p)) + y1(p);
                plot(xvals,y_perp(:,p), 'b--'); leg3 = "perp. line";
                scatter(x_prev(calc_points(p, 1)), y_prev(calc_points(p,1)), 'r', 'filled')
                scatter(x_curr(calc_points(p, 2)), y_curr(calc_points(p,2)), 'k', 'filled')
                
                p
            end
            plot(obs_coords_im(1),obs_coords_im(2), 'p','DisplayName', 'OBS Station', 'MarkerFaceColor','blue','Marker', '^','MarkerSize',10);
            plot(199, 251, 'p', 'DisplayName', 'Dongsha', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);
        end

        if verbose
            disp(['Average velocity between SVG files: ' num2str(velocity) 'm/s'])
            disp(['next!'])
        end
    end
end
