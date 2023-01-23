function data8ord = c_efw_burst_order_signals(data8, varsb, nmdata, pha, ref_probep)
%C_EFW_BURST_ORDER_SIGNALS  Restore correct order of burst signals
%
% DATA8ORD = C_EFW_BURST_ORDER_SIGNALS(DATA8,VARSB,NMDATA,ATWO,REF_PROBEP)
%
% Restore the correct order of the EFW internal burst siglals by comparing
% them to the normal mode data
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

data8ord = data8;

bpha = c_phase(data8(:,1),pha);
if isempty(bpha)
  irf_log('proc','Cannot check burst order: no phase for IB data')
  return;
end

nmdata = irf_tlim(nmdata,data8(1,1),data8(end,1));
if size(nmdata,1)<2
  irf_log('proc','Cannot check burst order: no/short reference NM data')
  return
end
nmpha=c_phase(nmdata(:,1),pha);
if isempty(nmpha)
  irf_log('proc','Cannot check burst order: no phase for NM data')
  return;
end
nmsfit = c_efw_sfit(ref_probep,3,10,20,nmdata(:,1),nmdata(:,2),...
  nmpha(:,1),nmpha(:,2),1,'hx');

nv = length(varsb);

% We always start from the pair which we think is OK
vi = 1;
for i=1:2:nv
  if str2double(varsb{i}(2))==floor(ref_probep/10) && ...
      str2double(varsb{i+1}(2))==rem(ref_probep,10)
    break
  else
    vi = vi+2;
  end
end

% ff gg di
bestguess=[-1 -1 -1];
for di=1:2:nv
  vidx = vi+di-1;
  if vidx+1 > nv, vidx = mod(vi+di-1,nv); end
  irf_log('proc',sprintf('testing %s and %s',...
    varsb{mod(vidx,nv)},varsb{vidx+1}))
  p1(:,1:2)=data8(:,[1 vidx+1]);
  p2(:,1:2)=data8(:,[1 vidx+2]);
  
  % First, check medians
  mp1 = median(p1(:,2)); mp2 = median(p2(:,2));
  if abs(mp1-mp2)/abs(mp1+mp2) > 1
    irf_log('proc','Medians dismatch')
    continue
  end
  
  if ref_probep==32, distance = 62; else, distance = 88;end
  eburst=1000*0.00212*(p2(:,2)-p1(:,2))/distance; % Burst electric field in mV/m
  bsfit=c_efw_sfit(ref_probep,3,10,20,p2(:,1),eburst,bpha(:,1),bpha(:,2),1,'ib');
  if isempty(nmsfit) || isempty(bsfit)
    %        irf_log('proc','sfit empty')
    continue
  end
  [ii1,ii2] = irf_find_comm_idx(nmsfit,bsfit);
  % Second, compare direction of E-fields
  y = abs(((atan2(nmsfit(ii1,2),nmsfit(ii1,3))-atan2(bsfit(ii2,2),bsfit(ii2,3)))/pi)*180);
  y(isnan(y)) = [];
  nbad =  - sum((25<y & y<155) | (205<y & y<335));
  ngood = length(y) - nbad;
  
  if (ngood>bestguess(1) && nbad==0)
    bestguess=[ngood nbad di];
    irf_log('proc','Data match')
  end
end

if bestguess(3)<0
  irf_log('proc','None of the signals match the ref data!!!')
else
  sh = bestguess(3) -1; % Shift
  if sh==0
    irf_log('proc','Order OK')
  else
    newidx = mod((1:nv)+sh-1,nv) + 1;
    data8ord(:,2:nv+1) = data8(:,newidx+1);
    irf_log('proc','reordering')
  end
end