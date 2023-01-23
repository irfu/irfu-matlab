function [res, ecorr] = c_efw_corrot(cl_id,diEs,diBr,P,R,SAX,diV)
%C_EFW_CORROT  compare EFW, EDI with CORROTATION
%
% [res, diECorr] = C_EFW_CORROT(cl_id,[diEs,diBr,P,R,SAX,diV])
%
% Compare E-fileds measured by EFW with CORROTATION E-filed.
%
% See also: CAA_COMP_EFW_EDI_CORR
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


% Control parameters
DE_LIM = 1;       % Limit on deviation from corrotation in mV/m
SCPOT_LIM = -1.5; % Limit for spacecraft potential
TAV = 60;         % Step in sec
WIN = 10;         % Averaging window WIN*TAV
R_LIM = 6*6378;   % Max R for which we apply the correction in km

if nargin==1
  diEs = c_load('diEs?',cl_id,'var');
  if diEs(1,1) == -157e8, error('No E-field'), end
  diBr = c_load('diBr?',cl_id,'var');
  if diBr(1,1) == -157e8
    diBr = c_load('diBrs?',cl_id,'var');
    if diBr(1,1) == -157e8, error('No B-field'), end
  end
  P = c_load('P?',cl_id,'var');
  R = c_load('R?',cl_id,'var');
  SAX = c_load('SAX?',cl_id,'var');
  diV = c_load('diV?',cl_id,'var');
end

diR = c_gse2dsi(R,SAX);
diRr = irf_resamp(diR,diBr);

if ~any( irf_abs(diRr,1) < R_LIM )
  irf_log('proc','outside the corrotation region')
  res = [];
  ecorr = [];
  return
end


% Earth rotation vector in DSI
geiOM = diBr;
geiOM(:,2:3) = 0;
geiOM(:,4) = 1;
%om = irf_gse2gei(geiOM,-1);
om = irf.geocentric_coordinate_transformation(geiOM,'gse>gei');
om(:,2:4) = om(:,2:4)*2*pi/86400;
diOMr = c_gse2dsi(om,SAX);

% Corrotation E-filed
vCorr = irf_cross(diRr,diOMr);
vCorr(:,2:4) = -vCorr(:,2:4);
idiECorr = irf_e_vxb(vCorr,diBr);

% SC motion induced E-filed in the inertial frame
diEi = irf_tappl(irf_cross(diBr,irf_resamp(diV,diBr)),'*1e-3*(-1)');

% Corrotation E-filed in the SC frame
diECorr = idiECorr;
diECorr(:,2:4) = diECorr(:,2:4) + diEi(:,2:4);

ndata = ceil((diEs(end,1) - diEs(1,1))/TAV);
t = diEs(1,1) + (1:ndata)*TAV - TAV/2; t = t';

diEs(isnan(diEs(:,2)),:) = []; % NaNs can cause nasty problems
if length(diEs) < 2
  irf_log('proc','Not enough E-field data.')
  res = [];
  ecorr = [];
  return
end
diEr = irf_resamp(diEs,t,'fsample',.1/TAV);
diECr = irf_resamp(diECorr,t,'fsample',.1/TAV);
Pr = irf_resamp(P(~isnan(P(:,2)),:),t,'fsample',1/TAV);
diRrr = irf_resamp(diRr,t,'fsample',1/TAV);
diRrr = irf_abs(diRrr,1);
if length(diEr) < 2 || length(diECr) < 2 || length(Pr) < 2
  irf_log('proc','Not enough data after resampling.')
  res = [];
  ecorr = [];
  return
end

dE = sqrt( (diEr(:,2) - diECr(:,2)).^2 + (diEr(:,3) - diECr(:,3)).^2 );
idx = find( Pr(:,2) > SCPOT_LIM & dE > DE_LIM & diRrr < R_LIM);
clear diEr diECr Pr

if ~isempty(idx)
  irf_log('proc',sprintf('max diff from corrotation is %.2f mV/m',max(dE)))
  DT2 = WIN*TAV/2;
  res = t(idx(1)) + [-DT2 DT2];
  idx(1) = [];
  if ~isempty(idx)
    for j=idx'
      if res(end,2)>=t(j)-DT2, res(end,2) = t(j) + DT2;
      else, res = [res; t(j)-DT2 t(j)+DT2];
      end
    end
  end
else
  irf_log('proc','no plasmaspheric wakes')
  res = [];
end

if nargout>1, ecorr = diECorr; end
