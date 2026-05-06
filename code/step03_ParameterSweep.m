% Amanda Syamsul
% March 2, 2026 -- Parameter sweep 
close all; clear all; clc;

% Relevant output: Figure S10 and S11

%% Part 1. Matching A and T by varying g' coefficient, h1, and a
close all; clear all; clc

g=9.8; % gravity - never change
rho=1e3; % density of water - never change
H=619; % total depth 

% h1_grid=[70:10:160];
h1_grid=120;
coef=0.0003;
% a_grid = [-180:10:-50];
a_grid=-130;
cmap=cool(length(a_grid));  
% figure()

for aa = 1:length(a_grid)
    a = a_grid(aa);

    for hh = 1:length(h1_grid)
        h1 = h1_grid(hh);
        h2=H-h1; % lower layer
        gp=g*coef; 
    
        s=[-1.5e3:1.5e3]; % x coordinate of the region (meters)
        [c0store,kstore,alpha,beta,I4,I5]=ComputeCoefficients(h1,h2,gp,g); 
        % then pressure
        [pH,pSH,pNH,L,T,c1,c0]=ComputePressure(s,h1,h2,a,gp,g,rho,kstore,I4,I5,c0store,alpha,beta);
        Total= pH+pSH+pNH;
    
        nn = 320;
        con_window=sin(pi*(1:nn)./nn);
        con_window=con_window./sum(con_window);
        
        T_smooth = conv(Total,con_window,'same');

        crest_idx = find(s==0);
        crest = T_smooth(crest_idx);
        before_peak = T_smooth(1:crest_idx);
        [trough, trough_idx] = min(before_peak);
        s_trough = s(trough_idx);
        s_crest = s(crest_idx);
        L_new = s_crest-s_trough;
        T_grid(hh) = L_new/c1; 
        A_grid(hh) = max(T_smooth)-min(T_smooth); % crest minus trough
    
    end
   
end

% xlim([250 900])
% ylim([0 1e3])
% legend()


%%
figure()

plot(s,pH,s,pSH,s,pNH,'LineWidth',1)
hold on
plot(s,Total,'--', 'LineWidth',2)
plot(s,T_smooth, 'LineWidth',2)
hold off
legend('p_{hydrostatic}','p_{surface hydrostatic}','p_{non-hydrostatic}','p_{Total}', 'p_{Total smoothed}')
ylabel('Pressure (Pa)')
grid on
%pause
% [max(Total) T]

%% Part 2. Matching c(z) by varying H, h1, a
close all; clear all; clc

rho=1e3; % density of water - never change
g=9.8; % gravity - never change
depth = linspace(475, 1200, 10);

a_grid=-130;
h1_grid=120;
coef=0.0003;
gp=g*coef;

for aa = 1:length(a_grid)
    a = a_grid(aa)
    for hh = 1:length(h1_grid)
        h1 = h1_grid(hh);
    
        for zz = 1:length(depth)
            
            H=depth(zz); % total depth (change this for fig 3 only)
            h2=H-h1; % lower layer

            s=[-1.5e3:1.5e3]; % x coordinate of the region (meters)
            [c0store,kstore,alpha,beta,I4,I5]=ComputeCoefficients(h1,h2,gp,g); 
            [pH,pSH,pNH,L,T,c1,c0]=ComputePressure(s,h1,h2,a,gp,g,rho,kstore,I4,I5,c0store,alpha,beta);
            % c1 = nonlinear velocity
            c1_grid(zz) = c1; 
        end
    end
    % params = table(h1_grid', c1_grid(:,1), c1_grid(:,2), c1_grid(:,3), ...
    % 'VariableNames', {'h1', 'c475','c675','c875'});
end

% figure()
hold on
scatter(-depth, c1_grid, 'filled', 'DisplayName','a=-180, h1=140')
xlabel('depth (m)')
ylabel('nonlinear wave speed, c_1')

disp(['min c(z) = ', num2str(min(c1_grid))]);
disp(['max c(z) = ', num2str(max(c1_grid))]);

%% Figure S10 - pH/pSH ratio dependence on h1 and a

g=9.8; % gravity - never change
rho=1e3; % density of water - never change
H=619; % total depth 

h1_grid=[30:10:150];
coef = 0.0002;
a_grid = [-180:10:-10];

r_grid = zeros(length(h1_grid), length(a_grid));

for aa = 1:length(a_grid)
    a = a_grid(aa);

    for hh = 1:length(h1_grid)
        h1 = h1_grid(hh);
        h2 = H - h1;
        gp = g * coef; 
    
        s = [-1.5e3:1.5e3];
        [c0store,kstore,alpha,beta,I4,I5] = ComputeCoefficients(h1,h2,gp,g); 
        
        [pH,pSH,pNH,L,T,c1,c0] = ComputePressure( ...
            s,h1,h2,a,gp,g,rho,kstore,I4,I5,c0store,alpha,beta);

        ratio = min(pH) / max(pSH);

        % store per (h1, a)
        r_grid(hh, aa) = ratio;
    end
end


[A, H1] = meshgrid(a_grid, h1_grid); 

figure()
scatter(A(:), H1(:), 100, r_grid(:), 'filled')
xlim([min(A(:)) - 10, max(A(:)) + 10])
ylim([min(H1(:)) - 10, max(H1(:)) + 10])
xlabel('a (m)')
ylabel('h_1 (m)')
title('p_h/p_{sh}')

%% Figure S11 -  Showing how T varies with h1 (Eqn. 52 of SI)

H = 619;
h1 = 70:10:160;
h2 = H - h1;
a = -130;
g = 9.8;
coef = 0.0003;
gp = g * coef;

tau = sqrt( (4 .* h1 .* h2 .* H) ./ (3 .* a .* (h1 - h2) .* gp) ) ...
    .* ( (2 .* h1 .* h2) ./ ( (2 .* h1 .* h2) + a .* (h1 - h2) ) );

figure()
plot(h1, tau, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'auto')
xlabel('Mixed layer depth, h_1 (m)')
ylabel('Period, \tau (s)')
grid on

%%
figure()
scatter(H_grid, c1_grid, 'filled')%, 'DisplayName',num2str(kstore(ks)))
grid on
xlabel('Depth (m)')
ylabel('c_0')

%% Visualization figs
% figure(1)
subplot(2,1,1)
semilogx(kstore,alpha)
ylabel('\alpha')
xlabel('k')
subplot(2,1,2)
semilogx(kstore,beta)
xlabel('k')
ylabel('\beta')

%%
figure(2)
semilogx(kstore,c0store,'o')
xlabel('k')
ylabel('c')
title(['density contrast = ' num2str(gp/g)])


%%

function [c0store,kstore,alpha,beta,I4,I5]=ComputeCoefficients(h1,h2,gp,g)
    % computes c0,alpha, beta for each wavenumber k for free surface case
    H=h1+h2;
    dz=0.1;
    z=[0:dz:H];
    ii=0;
    for ik=-4:0.01:-1, % this creates steps, and we can smooth out using 0.001
        ii=ii+1;
        k=10^ik;
        kstore(ii)=k;
        
        % For each wavenumber k, we compute the linear wave speed c0 that a small-amplitude internal wave would travel at in a two-layer ocean with a free surface
        % hyperbolic tangents for later
        T1=tanh(k*h1);
        T2=tanh(k*h2);
    
        % Ac^4+Bc^2+C=0 (quadratic eqn where x=c^2)
        A=k^2*(1+T1*T2);
        B=-k*(g*(T1+T2)+gp*T2);
        C=g*gp*T1*T2;
        c02=(-B-sqrt(B^2-4*A*C))/(2*A);  % + is surface. - is internal
        c0=sqrt(c02); 
        c0store(ii)=c0;
        beta1=-g/(k*c02); % this is not the same as beta (dispersion)
        
        % phi = vertical structure (how horizontal v changes with depth)
        phi1 = (cosh(k*(H-z))+beta1*sinh(k*(H-z)))./(cosh(k*h1)+beta1*sinh(k*h1)); % top layer
        phi2 = sinh(k*z)./(sinh(k*h2)); % bottom layer
        II2=(z<=h2); % idx for bottom layer
        II1=(z>h2); % idx for top layer
        phi=[phi2(II2) phi1(II1)];
    
        % z at the end means dphi/dz
        phi1z=(-k*sinh(k*(H-z))-beta1*k*cosh(k*(H-z)))./(cosh(k*h1)+beta1*sinh(k*h1)); % analytical derivative
        phi2z=k*cosh(k*z)./(sinh(k*h2));
        phiz=[phi2z(II2) phi1z(II1)];
        
        % integrals
        I1=trapz(phiz.^3)*dz;
        I2=trapz(phiz.^2)*dz;
        I3=trapz(phi.^2)*dz;
        
        % I4 and I5 are for pressure calculations
        I4(ii)=trapz(phi)*dz;
        I5(ii)=trapz(phi.*phiz)*dz;
        
        % alpha = nonlinearity, beta = dispersion
        % we are assuming negligible current velocity (ie. U=0)
        alpha(ii)=(3/2)*c0*I1/I2; % note sign switch from Moum and Nash U-c0 
        beta(ii)=(1/2)*c0*I3/I2;
    end
end

%%
function [pH,pSH,pNH,L,T, c1, c0]=ComputePressure(s,h1,h2,a,gp,g,rho,kstore,I4,I5,c0store,alpha,beta)
    H=h1+h2;
    L_array=sqrt(12*beta./(a.*alpha)); % equation before 8.17 Gerkema 
    [~,ik]=min(abs(1./kstore-L_array));
    L=L_array(ik); % wavelength
    k=1/L;
    
    % Shape function
    sech_sL = 1./cosh(s/L);
    A = a * sech_sL.^2;
    As = (-2*a/L) .* sech_sL.^2 .* tanh(s./L);
    Ass = (4*a/L^2) * sech_sL.^2 - (6*a/L^2)*sech_sL.^4;
    
    % shallow-water dispersion !! (simplification consistent with calculation of L)
    %Bref=-(g*H+gp*h2);
    %Cref=g*gp*h1*h2;
    %c0 = 0.5*(-Bref - sqrt(Bref^2 - 4*Cref));
    % no longer shallow water
    c0=c0store(ik);
    %c1shallow = c0 * (1 + a*(h1-h2)/(2*h1*h2));
    c1=c0+a*alpha(ik)/3; % c1 = nonlinear velocity, a = amplitude
    %[c1 c0 c1shallow]
    
    % Pressure
    I4_k=I4(ik);
    I5_k=I5(ik);
    
    D = cosh(k*h1) + (-g/(c0^2*k)) * sinh(k*h1); % required to incorporate free surface
    pH = rho * A * gp;
    pSH = rho * A * g/D;
    pNH = rho * c1^2 * (Ass * I4_k );%+ (As.^2 - A.*Ass)*I5_k); % negligble last term (advective)
    T=L/c1;
end