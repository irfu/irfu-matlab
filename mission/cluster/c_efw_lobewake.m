function [res,dEx] = c_efw_lobewake(cl_id,diE,diB,Ps,R,diEDI,ecorr)
%C_EFW_LOBEWAKE  detect lobe wakes
%
% [wake,dEx] = C_EFW_LOBEWAKE(cl_id,[diEs,diBrs,Ps,R,diEDI, e_corr_flag])
%
% Detect wakes using conditions on large parallel and penrpendicular
% electric fields. X-tra Ex offset (dEx) is computed from difference
% between EFW and EDI/Ex=0.
%
% e_corr_flag = 1, use dEx(X-tra offset) when detecting wakes
%
% C_EFW_LOBEWAKE(cl_id,'plot')
%
% Do a plot illustrating wake detection

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% Control parameters
TAV = 60;              % averaging window in sec
ANG_LIM = 15;
EPAR_LIM = 1.0;        % limit on E|| in mV/m
EPERP_LIM = 1.5;       % limit on full Eperp (E*B=0) in mV/m
EDEV_LIM = .8;         % limit on s-dev of E
EZ_LIM = 1;            % limit on Ez in mV/m
EPAR_EPERP_RATIO_LIM = .2;
EPERP_RATIO_LIM = 1.8; % Ratio between full Eperp (E*B=0) and
% E along the projecttion of B to spin plane
SCPOT_LIM = -6;        % Limit for spacecraft potential
DEX_MIN = 0.0;         % minimum reasonable dEx
DEX_MAX = 2.0;         % maximum reasonable dEx
EDI_COVERAGE_MIN = 0.5;% minimum EDI coverage

if ischar(diE) && strcmp(diE,'plot'), DOPLOT = 1;
else, DOPLOT = 0;
end

% Load data
if nargin==1 || ( ischar(diE) && strcmp(diE,'plot') )
  ecorr = 0;
  spinFits = caa_sfit_load(cl_id);
  if isempty(spinFits), error('cannot load E'), end
  diE = spinFits.diEs;
  [ok,Ps] = c_load('Ps?',cl_id);
  if ~ok, error('cannot load P'), end
  [Ddsi,Damp] = c_efw_dsi_off(diE(1,1),cl_id,Ps);
  diE = caa_corof_dsi(diE,Ddsi,Damp);
  clear Ddsi Damp
  [ok,diB] = c_load('diBrs?',cl_id);
  if ~ok, error('cannot load B'), end
  [ok,diEDI] = c_load('diEDI?',cl_id);
  if ~ok, disp('cannot load EDI'), end
  [ok,R] = c_load('R?',cl_id);
  if ~ok, disp('cannot load R'), end
elseif nargin<7, ecorr = 0;
end

% In the electric is already corrected, adjust ranges for the offset
if ecorr > 0
  DE_RANGE = (DEX_MAX-DEX_MIN)/2.0;
  % XXX: we do not change this is to help situations with ASPOC on
  %DEX_MAX = DEX_MAX - DE_RANGE;
  DEX_MIN = DEX_MIN - DE_RANGE;
  clear DE_RANGE
end

% Remove NaN
diE=diE(isfinite(diE(:,2)),:);
if size(diE,1) < 1
  ndata = 0;
else
  ndata = ceil((diE(end,1) - diE(1,1))/TAV);
end
if ndata==0
  irf_log('proc','not enough data to compute wakes')
  res = []; dEx = [];
  return
end
t = diE(1,1) + (1:ndata)*TAV - TAV/2; t = t';

diEr = irf_resamp(diE,t,'fsample',1/TAV);

diEDIr = [];
if exist('diEDI','var') && ~isempty(diEDI)
  diEDI = diEDI(~isnan(diEDI(:,2)),:);
  if size(diEDI(:,2),1) > EDI_COVERAGE_MIN*size(diE,1)
    diEDIr = irf_resamp(diEDI,t,'fsample',1/TAV,'thresh',1.3);
  else
    irf_log('dsrc','bad EDI coverage')
  end
end

% Try to empirically correct the sunward offset
dEx = [];
if ~isempty(diEDIr)
  ii = find( ~isnan(diEr(:,5)) & (diEr(:,5) < EDEV_LIM) & ...
    ~isnan(diEDIr(:,2)));
  if length(ii)>2
    dEx = mean(diEr(ii,2) - diEDIr(ii,2));
    irf_log('proc',sprintf('x-tra (EDI) dEx : %.2f mV/m',dEx))
  end
end
if ~isempty(dEx) && ( dEx < DEX_MIN || dEx > DEX_MAX )
  irf_log('proc','not using x-tra dEx')
  dEx = [];
end
% Alternatively we can try zero Ex condition
% Probably this can be applied only in the far tail where one expects
% zero Ex
if isempty(dEx) && ~isempty(R)
  Rr = irf_resamp(R,t,'fsample',1/TAV); Rr = irf_abs(Rr);
  % Look for negative X and distances > than 5 RE, RE = 6378 km
  if isempty(Rr)
    ii=[];
    dEx=[];
    irf_log('proc','Short R data. No dEx.')
  else
    ii = find( Rr(:,2) < 0 & Rr(:,5) > 5*6378 );
  end
  %ii = 1:length(diEr(:,1));
  if ~isempty(ii)
    dEx = mean( diEr(~isnan(diEr(ii,2)),2) );
    s = std( diEr(~isnan(diEr(ii,2)),2) );
    if s>0
      dEx = mean( diEr( diEr(ii,2)>dEx-s & diEr(ii,2)<dEx+s ,2) );
    end
    irf_log('proc',sprintf('x-tra (Ex=0) dEx : %.2f mV/m',dEx))
  end
end
if ~isempty(dEx)
  if dEx > DEX_MIN && dEx < DEX_MAX
    irf_log('proc','correcting x-tra dEx')
    diEr(:,2) = diEr(:,2) - dEx;
    diE(:,2) = diE(:,2) - dEx;
  else
    irf_log('proc','not using x-tra dEx')
    dEx = [];
  end
end

%TODO: the offset needs to be reverified after we found the wakes

Psr = irf_resamp(Ps(~isnan(Ps(:,2)),:),t,'fsample',1/TAV);
diBr = irf_resamp(diB,t,'fsample',1/TAV);
bele = (180/pi)*asin(...
  diBr(:,4)./sqrt(diBr(:,2).^2 + diBr(:,3).^2 + diBr(:,4).^2) );
if DOPLOT
  clf
  nplots = 6;
  h = 1:nplots;
  h(1) = irf_subplot(nplots,1,-1);
  irf_plot(Ps);
  ylabel('Sc pot [-V]')
  st_iso = epoch2iso(diE(end,1));
  title(['Cluster ' num2str(cl_id) ', ' st_iso(1:10)])
  h(2) = irf_subplot(nplots,1,-2);
  irf_plot(diE(:,[1 2]));
  set(gca,'YLim',[-9 9])
  ylabel('Ex [mV/m]')
  h(3) = irf_subplot(nplots,1,-3);
  irf_plot(diE(:,[1 3]));
  set(gca,'YLim',[-9 9])
  ylabel('Ey [mV/m]')
  h(4) = irf_subplot(nplots,1,-4);
  if exist('diEDI','var') && ~isempty(diEDIr)
    axes(h(2)), hold on
    irf_plot(diEDIr(:,[1 2]),'k.');
    axes(h(3)), hold on
    irf_plot(diEDIr(:,[1 3]),'k.');
    diEDIr = diEDIr(:,1:3);
    diEDIr(:,2:3) = diEr(:,2:3) - diEDIr(:,2:3);
    axes(h(4))
    irf_plot(diEDIr)
    ylabel('E diff [mV/m]')
  end
  h(5) = irf_subplot(nplots,1,-5);
  irf_plot([diEr(:,1) bele]);
  ylabel('\theta [deg]')
end

wind = []; ind = [];
% Look for strong apparent Epar if B field is within 15 deg of spin plane
if ~isempty(Psr)
  ind = find( abs(bele) < ANG_LIM & Psr(:,2) < SCPOT_LIM);
  Epar = ( diEr(:,2).*diBr(:,2) + diEr(:,3).*diBr(:,3) )...
    ./sqrt( diBr(:,2).^2 + diBr(:,3).^2 );
  Eperp = abs( diEr(:,2).*diBr(:,3) - diEr(:,3).*diBr(:,2) )...
    ./sqrt( diBr(:,2).^2 + diBr(:,3).^2 );
  wind = ind(abs(Epar(ind)) > EPAR_LIM &...
    abs(Epar(ind)./Eperp(ind)) > EPAR_EPERP_RATIO_LIM);
end

if DOPLOT && ~isempty(wind)
  axes(h(2)), hold on
  irf_plot(diEr(wind,[1 2]),'g.');
  axes(h(3)), hold on
  irf_plot(diEr(wind,[1 3]),'g.');
end


% Look for strong apparent Ez if B field is above 15 deg of spin plane
if ~isempty(Psr)
  ind = find( abs(bele) >= ANG_LIM & Psr(:,2) < SCPOT_LIM);
  if ~isempty(ind)
    diEr(:,4) = -(diEr(:,2).*diBr(:,2)+diEr(:,3).*diBr(:,3))./diBr(:,4);
    diEr = irf_abs(diEr); % diEr(:,6) now contains abs(E)
    ind = ind(abs(diEr(ind,4)) > EZ_LIM & diEr(ind,6) > EPERP_LIM &...
      diEr(ind,6)./Eperp(ind) > EPERP_RATIO_LIM & diEr(ind,5) < EDEV_LIM);
  end
end

if DOPLOT && ~isempty(ind)
  axes(h(2)), hold on
  irf_plot(diEr(ind,[1 2]),'r.');
  axes(h(3)), hold on
  irf_plot(diEr(ind,[1 3]),'r.');
end

wind = sort([wind; ind]);
clear ind

if DOPLOT
  for j=wind'
    %diE( diE(:,1)>=t(j)-TAV/2 & diE(:,1)<=t(j)+TAV/2, 2:end) = NaN;
    diE( diE(:,1)>=t(j)-TAV & diE(:,1)<=t(j)+TAV, 2:end) = NaN;
  end

  h(6) = irf_subplot(nplots,1,-6);
  irf_plot(diE(:,1:3))
  ylabel('E cor [mV/m]')

  irf_zoom(h,'x',[t(1)-TAV/2 t(end)+TAV/2])

  [ok,r] = c_load('R?',cl_id);
  if ok, add_position(h(6),r), end
elseif ~isempty(wind)
  DT2 = TAV;
  res = t(wind(1)) + [-DT2 DT2];
  wind(1) = [];
  if ~isempty(wind)
    for j=wind'
      if res(end,2)>=t(j)-DT2*2, res(end,2) = t(j) + DT2; % throw away one good point in between bad (DT2*2)
      else, res = [res; t(j)-DT2 t(j)+DT2]; %#ok<AGROW>
      end
    end
  end
else
  irf_log('proc','no lobe wakes')
  res = [];
end
