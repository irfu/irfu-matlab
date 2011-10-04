function [apar,aperp]=irf_dec_parperp(b0,a)
%IRF_DEC_PARPERP   Decompose a vector into par/perp to B components
%
% [apar,aperp]=irf_dec_parperp(b0,a)
%
%	Decomposes A to parallel and perpendicular to BO components
%
%	b0,a - martixes A=(t,Ax,Ay,Az) // AV Cluster format
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Yuri Khotyaintsev, 1997
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

btot = irf_abs(b0,1);

ii = find(btot<1e-3);
if ~isempty(ii), btot(ii) = ones(size(ii))*1e-3; end
normb = [b0(:,1) b0(:,2)./btot b0(:,3)./btot b0(:,4)./btot]; 
normb = irf_resamp(normb,a);

apar = irf_dot(normb,a);
aperp = a;
aperp(:,2:4) = a(:,2:4) - normb(:,2:4).*(apar(:,2)*[1 1 1]);

return
