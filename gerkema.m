%% eqn 8.7
% 
d = 30;
H = [0:10:500];
g = 9.8;

rho = 0.001; 
x = 0.001 % ;percent change
g_prime = g*x;
c2 = g_prime .* ( d.*(H - d) ./ H );

figure()
plot(-H, c2)

%% eqn 8.17

% eta = a * sech(x - Ct / l)^2
% C = c_0 * ( 1 +  ( a( (h1-h2) / (h1*h2) ) / 2 )
% l = 2*h1*h2 / (3a(h1-h2) ) ^ (1/2)
% c_0 ^2 = g'(h1h2 / (h1+h2) )

%% Fig 1: Approach 1
% clear; close all; clc

a = -100;          % amplitude (m), depression for h1<h2
H = 619;          % total depth (m)
h1 = linspace(1, 100, 200);   % upper layer depth
h2 = H - h1;
g = 9.81;

rho1 = 1020;                            % upper-layer density
delta_rho = [0.2 0.5 1.0 2.0 3.0];      % contrasts to test
rho2 = rho1 + delta_rho;

figure();
for i = 1:length(delta_rho)

    gp = g * (delta_rho(i) / rho2(i));   % reduced gravity?

    % Linear wave speed
    c0 = sqrt( gp .* ( (h1.*h2) ./ (h1+h2) ) );
    % Nonlinear phase speed
    C = c0 .* ( 1 + a .* ( (h1-h2) ./ (h1.*h2) )/2 );
    
    % Soliton width
    bad = a.*(h1 - h2) <= 0;
    ell = 2*h1.*h2 ./ sqrt(3*a.*(h1 - h2));  % meters
    ell(bad) = NaN;

    tau = ell ./ C;
    tau(bad) = NaN;
    semilogy(h1,tau/60, 'DisplayName', ['gp = ', num2str(gp), ', \Delta\rho = ', num2str(delta_rho(i))]' )
    hold on
end
legend()
xlabel('Shallow layer depth (m)')
ylabel('Temporal half-width (mins)')

xlabel('Upper-layer thickness h_1 (m)')
title(sprintf('Half-width vs stratification (H=%.0f m, a=%.1f m', H, a))
grid on


%% Fig 2: Approach 2

close all;clc;

h1vect=[30:10:50];          % upper layer depth
for i=1:length(h1vect)

    h1 = h1vect(i);
    a = [-100:-5];          % amplitude (m), depression for h1<h2
    H = 619;          % total depth (m)
    h2 = H - h1;
    g=9.81;
    fracrho=1e-4;
    gp = g*fracrho;        % reduced gravity?
    
    % Linear wave speed
    c0 = sqrt( gp .* ( (h1.*h2) ./ (h1+h2) ) );
    % Nonlinear phase speed
    C = c0 .* ( 1 + a .* ( (h1-h2) ./ (h1.*h2) )/2 );
    
    % Soliton width
    bad = a.*(h1 - h2) <= 0;
    ell = 2*h1.*h2 ./ sqrt(3*a.*(h1 - h2));  % meters
    ell(bad) = NaN;
    
    % Temporal half-width (seconds)
    tau = ell ./ C;
    
    plot(abs(a),tau/60, 'LineWidth',1.5,'DisplayName',['Upper layer thickness = ', num2str(h1vect(i)), 'm'])
    hold on
    xlabel('Amplitude (m)')
    ylabel('Temporal half-width (mins)')
end
legend()
grid on

%% Fig 3: Color map of temporal half-width vs h1 and Δρ

[h1g, drhog] = meshgrid(h1, delta_rho);
h2g = H - h1g;

% Reduced gravity for each Δρ
gp = g .* (drhog ./ (rho1 + drhog));

% Linear wave speed
c0 = sqrt(gp .* ((h1g .* h2g) ./ (h1g + h2g)));

% Nonlinear phase speed
C = c0 .* (1 + a .* ((h1g - h2g) ./ (2 .* h1g .* h2g)));

% Soliton width
bad = a .* (h1g - h2g) <= 0;
ell = 2 .* h1g .* h2g ./ sqrt(3 .* a .* (h1g - h2g));
ell(bad) = NaN;

% Temporal half-width (minutes)
tau = (ell ./ C) ./ 60;
tau(bad) = NaN;

% Plot
figure;
contourf(h1, delta_rho, tau, 30, 'LineColor', 'none');
xlim([0 300])
caxis([0 100])
colormap(parula)
colorbar;
xlabel('Upper-layer thickness h_1 (m)');
ylabel('\Delta\rho (kg/m^3)');
title(sprintf('Temporal half-width (min) vs h_1 and \\Delta\\rho (H=%.0f m, a=%.1f m)', H, a));
grid on;
