% test script for C_EFW_C_EFW_SPINFIT_MX
%
% $Id$

fomega = pi/2;
fnterms = 3;
fitmax = 10;
t=0:255;
at = t*.0156;
az = .1 + .5*sin(fomega*at);

[bad,x,sigma,iter,lim] = c_efw_spinfit_mx(fnterms,fitmax,fomega,at,az);

disp(x(1:fnterms))
size(x)
