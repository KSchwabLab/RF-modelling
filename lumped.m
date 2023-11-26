
function [Q, wres, R, L, C, Z_lump] = lumped(w, realZin)
[R,i0] = max(realZin);
wres = w(i0);
wplus = interp1(realZin(i0:end),w(i0:end),R/2);
wminus = interp1(realZin(1:i0),w(1:i0),R/2);
C = 1/R/(wplus-wminus);
L = 1/(wres^2*C);
Z_lump = 1./R./(1/R^2+(w*C-1./(w*L)).^2);
wres = w(i0);
Q = R*sqrt(C/L);

end