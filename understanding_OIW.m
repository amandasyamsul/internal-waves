clear all; close all; clc

depth = linspace(0, 700, 100); 
density = 1029 + -(exp(-depth / 50)); 

figure;
plot(density, depth, 'k');

xlabel('Density, \rho (kg/m^3)');
ylabel('Depth, z (m)');

title('(a)');

xlim([1027.5 1029.5]);
ylim([-50 700]);

% Show the plot
grid on;
axis ij


%%

g = 9.81;  
rho_0 = 1029; 

d_rho_dz = gradient(density, depth);

N2 = (g / rho_0) * d_rho_dz;

% Calculate N (Brunt-Väisälä frequency)
N = sqrt(N2);

figure;
plot(N, depth, 'b'); 

xlabel('Brunt-Väisälä Frequency, N (s^{-1})');
ylabel('Depth, z (m)');

title('Brunt-Väisälä Frequency Profile');
xlim([min(N) max(N)]);
ylim([0 700]);

% Invert the y-axis for depth
axis ij
grid on;


%%

g = 9.81;  
rho_0 = 1029; 

h1 = 100;  % Upper layer depth (m)
h2 = 700 - h1;  % Lower layer depth (m)

delta_rho = max(density) - min(density); 

% Long-wave speed c0 (Eqn 8)
c0 = sqrt((g * delta_rho / rho_0) * (h1 * h2) / (h1 + h2));

% Nonlinear coefficient alpha (Eqn 8)
alpha = (3/2) * ( (h1 - h2) / (h1 * h2) ) ;

% Estimate eta_0 (the vertical displacement amplitude, an approximation)
eta_0 = 20;

% Calculate the nonlinear phase speed V (Eqn 14)
V = c0 * (1 + ( (2/3) * alpha * eta_0));

% Display the results
fprintf('Long-wave phase speed (c0): %.2f m/s\n', c0);
fprintf('Nonlinear phase speed (V): %.2f m/s\n', V);
