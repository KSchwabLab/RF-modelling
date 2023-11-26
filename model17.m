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
alpha = beta*2e-5; % is loss angle of 2e-5 correct?
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

figure;
plot(freq/1e6,realZin,'LineWidth',1.5)
hold on;

[Q_resonator, wres, R, L, C, Z_lump] = lumped(w, realZin);


A2 = 0.005^2/4*pi;
d2 = 1e-4+1e-9;
C_detector2 = eps0*A2/d2;
Zl2 = 1./(1i.*w*C_detector2);
Zin2 = Z0.*(Zl2.*cosh(gamma.*l)+Z0.*sinh(gamma.*l))./(Z0.*cosh(gamma.*l)+Zl2.*sinh(gamma.*l));
realZin2 = real(Zin2);

[Q_resonator2, wres2, R2, L2, C2, Z_lump2] = lumped(w, realZin2);

g = wres2-wres;

C_coupler = 1/sqrt(R*Z0*wres^2);

kint = 1/R/C;
kout = wres^2*R0*C_coupler^2/C;
kappa = kint+kout;

sideband = 2*g/kappa;
S21 = 2*R/(R+R0-1i/w0/C_coupler);
atR0 = R0/(R0-1i/w0/C_coupler);
output_sideband = abs(sideband*S21*atR0);

Clist = 1e-15:0.1e-15:10e-15;
for jj = 1:length(Clist)

    kint = 1/R/C;
    kout = wres^2*R0*Clist(jj)^2/C;
    kappa = kint+kout;
    
    sideband = 2*g/kappa;
    S21 = 2*R/(R+R0-1i/w0/Clist(jj));
    atR0 = R0/(R0-1i/w0/Clist(jj));
    output_sideband = abs(sideband*S21*atR0);
    temp(jj) = output_sideband;
end
figure;
plot(Clist*1e15, temp,'LineWidth',1.5)
xline(5.035,'LineWidth',1.5)
xlabel('Capacity of coupler (fF)')
ylabel('S11')
fontsize(16,'points')