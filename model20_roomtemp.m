clear;
close all;
%% Modelling the Q factor
Z0 = 50;
R0= 50;
c=3e8*0.7;  
eps0 = 8.85e-12;
freq = 1.2e9:1e5:1.5e9;
w = freq*2*pi;
freq0 = 1.47e9;
w0 = freq0*2*pi;
beta0 = w0/c;
alpha0 = beta0*2e-4; % is loss angle of 2e-5 correct?
gamma = alpha0+beta0*1i;

lambda0 = c/1.47e9;
l=lambda0/2;

R = R0./(alpha0*l);
L = 2*R0/pi/w0;
C = pi/(2*Z0*w0);

Zin = Z0./alpha0/l./(1+1i*pi.*(w-w0)./alpha0/l/w0);
C_coupler = eps0*0.0005^2*pi/0.001;
C_coupler = 2e-13;
C_coupler_wanted = 1/sqrt(R*Z0*w0^2);

S11 = (Zin+1./(1i*w*C_coupler)-50)./(Zin+1./(1i*w*C_coupler)+50);

plot (abs(S11))

[Q_S11, ~, ~, ~, ~, ~] = lumped(w,1-abs(S11))
