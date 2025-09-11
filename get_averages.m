function averages = get_averages(T)
% Compute daily summary statistics for velocity and back azimuth
%
% INPUT:
%   T - a table with at least the following columns:
%       - datetime: datetime array
%       - velocity: numeric array
%       - back_azimuth: numeric array
%
% OUTPUT:
%   averages - table of daily averages, medians, std, max, min

    % Convert datetime to date (remove time)
    T.date = dateshift(T.datetime, 'start', 'day');

    % Group by date
    [G, dateGroups] = findgroups(T.date);

    % Compute daily statistics
    daily_avg_velocity      = splitapply(@nanmean,  T.velocity,      G);
    daily_median_velocity   = splitapply(@nanmedian,T.velocity,      G);
    daily_stdev_velocity    = splitapply(@nanstd,   T.velocity,      G);

    daily_avg_backazimuth   = splitapply(@nanmean,  T.back_azimuth,  G);
    daily_median_backazimuth= splitapply(@nanmedian,T.back_azimuth,  G);
    daily_stdev_backazimuth = splitapply(@nanstd,   T.back_azimuth,  G);
    daily_max_backazimuth   = splitapply(@nanmax,   T.back_azimuth,  G);
    daily_min_backazimuth   = splitapply(@nanmin,   T.back_azimuth,  G);

    % Package into output table
    averages = table(dateGroups, ...
                     daily_avg_velocity, daily_median_velocity, daily_stdev_velocity, ...
                     daily_avg_backazimuth, daily_median_backazimuth, ...
                     daily_stdev_backazimuth, daily_max_backazimuth, daily_min_backazimuth, ...
                     'VariableNames', {'datetime', ...
                                       'average_velocity','median_velocity','velocity_std_dev', ...
                                       'average_back_azimuth','median_back_azimuth','back_azimuth_std_dev', ...
                                       'max_ba','min_ba'});

    % Replace datetime group index with actual date
    unique_dates = unique(T.date);
    averages.datetime = unique_dates;

end
