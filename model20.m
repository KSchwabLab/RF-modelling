close all;
%% Modelling the Q factor
Z0 = 50;
R0= 50;
c=3e8*0.7;  
eps0 = 8.85e-12;
freq = 3.98e9:1e3:4.0e9;
w = freq*2*pi;
w0 = 4e9*2*pi;
beta = w/c;
alpha = beta*2e-5; 
gamma = alpha+beta*1i;

lambda0 = c/4e9;
l=0.01672;

%% Capacitor
A = 0.005^2/4*pi;
d = 1e-4;
C_detector = eps0*A/d;
Zl = 1./(1i.*w*C_detector);
Zin = Z0.*(Zl.*cosh(gamma.*l)+Z0.*sinh(gamma.*l))./(Z0.*cosh(gamma.*l)+Zl.*sinh(gamma.*l));
realZin = real(Zin);
figure
plot(freq/1e6,realZin)

[Q_resonator, wres, R, L, C, Z_lump] = lumped(w, realZin);


C_coupler = 1/sqrt(R*Z0*wres^2);
S11 = (Zin+1./(1i*w*C_coupler)-50)./(Zin+1./(1i*w*C_coupler)+50);
figure;
plot(freq/1e6,abs(S11))