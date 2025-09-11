fs = 20;

bathymetry1 = ncread("gebco_2024_n21.5_s20.0_w120.0_e123.0.nc", "elevation");
bathymetry = bathymetry1';
lon = ncread("gebco_2024_n21.5_s20.0_w120.0_e123.0.nc", "lon");
lat = ncread("gebco_2024_n21.5_s20.0_w120.0_e123.0.nc", "lat");

close all;

figure()

surf(lon, lat, bathymetry, 'EdgeColor', 'none');
colormap('turbo');
colorbar;        
shading interp;  
% view(2);         
xlabel('Longitude', 'FontSize', fs)
ylabel('Latitude', 'FontSize', fs)
title('Bathymetry of the South China Sea', 'FontSize', fs, 'FontAngle', 'italic');

figure()
contourf(lon, lat, bathymetry, 50, 'LineStyle', 'none');
hold on
colormap('turbo');
colorbar;
grid on;
legend('elevation', 'OBS')
legend show;
title('Bathymetric Contour Map','FontSize', fs, 'FontAngle', 'italic');
xlabel('Longitude', 'FontSize', fs)
ylabel('Latitude', 'FontSize', fs)
