function [apar,aperp]=decomposeParPerp(b0,a)
% function [apar,aperp]=decomposeParPerp(b0,a)
%
%	Decomposes A to parallel and perpendicular to BO components
%
%	b0,a - martixes A=(t,Ax,Ay,Az) // AV Cluster format
%
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Yuri Khotyaintsev, 1997
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndata = length(a(:,1));
normb = zeros(ndata,3);
btot = zeros(ndata,1);
apar = zeros(ndata,1);
aperp = zeros(ndata,3);

btot = av_abs(b0,1);

ii = find(btot<1e-3);
if length(ii)>0, btot(ii) = ones(size(ii))*1e-3; end
normb = [b0(:,1) b0(:,2)./btot b0(:,3)./btot b0(:,4)./btot]; 
normb = c_resamp(normb,a);

apar = av_dot(normb,a);
aperp = a;
aperp(:,2:4) = a(:,2:4) - normb(:,2:4).*(apar(:,2)*[1 1 1]);

return
