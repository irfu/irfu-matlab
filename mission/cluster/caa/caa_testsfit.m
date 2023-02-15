% test script for C_EFW_C_EFW_SPINFIT_MX
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

fomega = pi/2;
fnterms = 3;
fitmax = 10;
t=0:255;
at = t*.0156;
az = .1 + .5*sin(fomega*at);

[bad,x,sigma,iter,lim] = c_efw_spinfit_mx(fnterms,fitmax,fomega,at,az);

disp(x(1:fnterms))
size(x)
