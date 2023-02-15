function [dE_out,dAmp_out,dt_int,weight] = caa_sh_xoff_batch(st,dt,flag_amp)
%CAA_SH_XOFF  sunward offset and amplitude correction in the sh/sw
%
% [dE, dAmp, dt_int, weight] = caa_sh_xoff_batch(st,dt [,flag_amp])
% [dE, dAmp, dt_int, weight] = caa_sh_xoff_batch(iso_st,iso_et [,flag_amp])
%
% Study sunward offset (X GSE) and amplitude correction
% factor by comparing EFW data with CIS HIA
%
% if FLAG_AMP is zero (default), use default amplitude correction
% factor = 1.1
%
% See also CAA_SH_XOFF, CAA_SH_PLAN, CAA_COROF_DSI, CAA_SH_PL_XOFF
%

% Copyright 2007 Yuri Khotyaintsev

if nargin<3, flag_amp = 0; end
if flag_amp~=0, flag_amp = 1; end

if exist('./mPlan.mat','file'), load ./mPlan.mat
else, error('No MPlan.mat found')
end

MP = [];
[st,dt] = irf_stdt(st,dt); tt1 = fromepoch(st); tt2 = fromepoch(st+dt);
for yy=tt1(1):tt2(1)
  v_s = sprintf('MPauseY%d',yy);
  if ~exist(v_s,'var'), error(['Cannot load ' v_s]), end
  if isempty(MP), eval([ 'MP=' v_s ';'])
  else
    eval([ 'MP_tmp=' v_s ';'])
    % Check if years overlap
    if MP(end,2) == MP_tmp(1,2), MP_tmp(1,:) = []; end
    MP = [MP; MP_tmp]; clear MP_tmp %#ok<AGROW>
  end
end
clear tt1 tt2
MP = MP( MP(:,2)>st & MP(:,1)<st+dt ,:);

nint = size(MP,1);
dE = zeros(nint,4);
dAmp = zeros(nint,4);
weight = zeros(nint,1);
dt_int = MP(:,2) - MP(:,1);
for ii=1:nint
  figure(77), clf
  [dE(ii,:), dAmp(ii,:), weight(ii)] = ...
    caa_sh_xoff(MP(ii,1),dt_int(ii),flag_amp);
  orient tall
  print( 77, '-dpdf', [irf_fname(MP(ii,1)) '_XOFF'])
end

dE_out = [MP(:,1) dE];
dAmp_out = [MP(:,1) dAmp];