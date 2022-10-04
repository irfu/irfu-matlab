function [apar,aperp,alpha]=irf_dec_parperp(b0,a,flagspinplane)
%IRF_DEC_PARPERP   Decompose a vector into par/perp to B components
%
% [Apar,Aperp]=irf_dec_parperp(B0,A)
%
% Decomposes A into parallel and perpendicular to BO components
%
% [Apar,Aperp,Alpha_XY]=irf_dec_parperp(B0,A,1)
%
% Decomposes A into parallel and perpendicular components to the
% projection of B onto the XY plain. Alpha_XY gives the angle between B0
% and the XY plain.

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Yuri Khotyaintsev, 1997
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin temporary fix to convert TS format to older format
rtrnTS = 0;
if isa(b0,'TSeries')
  if isempty(b0), apar = TSeries([]); aperp = TSeries([]); alpha = TSeries([]); return; end
  b0Time = b0.time;
  datatemp = double(b0.data);
  b0 = [b0Time.epochUnix(), double(datatemp)];
end
if isa(a,'TSeries')
  if isempty(a), apar = TSeries([]); aperp = TSeries([]); alpha = TSeries([]); return; end
  aTime = a.time;
  datatemp = double(a.data);
  a = [aTime.epochUnix(), datatemp];
  rtrnTS = 1;
end

% End of temporary fix

if nargin<3 || flagspinplane==0
  btot = irf_abs(b0,1);
  
  ii = find(btot<1e-3);
  if ~isempty(ii), btot(ii) = ones(size(ii))*1e-3; end
  normb = [b0(:,1) b0(:,2)./btot b0(:,3)./btot b0(:,4)./btot];
  normb = irf_resamp(normb,a);
  
  apar = irf_dot(normb,a);
  aperp = a;
  aperp(:,2:4) = a(:,2:4) - normb(:,2:4).*(apar(:,2)*[1 1 1]);
  alpha = [];
  if rtrnTS
    aperp = irf.ts_vec_xyz(aTime,aperp(:,2:4));
  end
else
  irf_log('proc','Decomposing in the XY plane')
  b0 = irf_resamp(b0,a(:,1));
  btot = sqrt(b0(:,2).^2 + b0(:,3).^2);
  alpha = b0(:,1:2);
  alpha(:,2) = atan2d(b0(:,4),btot);
  b0(:,2) = b0(:,2)./btot; b0(:,3) = b0(:,3)./btot;
  apar = a(:,1:2); aperp = apar;
  apar(:,2) = a(:,2).*b0(:,2) + a(:,3).*b0(:,3);
  aperp(:,2) = a(:,2).*b0(:,3) - a(:,3).*b0(:,2);
  if rtrnTS
    aperp = irf.ts_scalar(aTime,aperp(:,2));
    alpha = irf.ts_scalar(aTime,alpha(:,2));
  end
end

if rtrnTS, apar = irf.ts_scalar(aTime,apar(:,2)); end

return
