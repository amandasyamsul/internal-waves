d = 30;
H = [0:10:500];
g = 9.8;

rho = 0.001; 
x = 0.001 % ;percent change
g_prime = g*x;
c2 = g_prime .* ( d.*(H - d) ./ H );

figure()
plot(-H, c2)