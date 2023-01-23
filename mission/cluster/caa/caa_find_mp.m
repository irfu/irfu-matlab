function [t_mp_out,t_mp_in] = caa_find_mp(start_time, dt, cl_id, Rin,sc_source)
%CAA_FIND_MP  find model magnetopause crossings using ACE data
%
% [t_mp_out,t_mp_in] = caa_find_mp(start_time, dt, cl_id)
%
% start_time: start time in epoch format
%         dt: duration in seconds or end time in epoch
%      cl_id: which Cluster sat
%        Rin: position of Cluster (optional)
%  sc_source: source for solar wind data, 'omni2' (default) or 'ace'
%
% See also IRF_SHUE_MP
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if dt>toepoch([1996 01 01 00 00 00])
  % et is given
  if dt< start_time, error('STOP_TIME must be larger then START_TIME)'), end
  dt = dt - start_time;
end
if nargin < 5, sc_source='omni2'; end
if nargin < 4, Rin = []; end
if ~strcmp(sc_source,'omni2') && ~strcmp(sc_source,'ace')
  error('Solar wind data source improperly specified. Should be ace or omni2.')
end

t_mp_out = []; t_mp_in = [];

R_E = 6378;
ACE_X_POS = 222*R_E;	% ACE X postition
ACE_VX_DEF = 480;		% Default solar wind speed
ACE_DT_DEF = ACE_X_POS/ACE_VX_DEF;
ACE_N_DEF = 6;			% Default solar wind density
ACE_BZ_DEF = 0;			% Default IMF Bz

irf_log('proc',['orbit : ' epoch2iso(start_time,1) ' -- ' ...
  epoch2iso(start_time+dt,1)])

if isempty(Rin)
  data = getData(ClusterDB, start_time, dt, cl_id, 'r', 'nosave');
  if isempty(data), error('cannot fetch position'), end
  R = data{2};
  clear data
else
  R = irf_tlim(Rin, start_time + [0 dt]);
  if isempty(R), irf_log('proc','empty position'), return, end
end

R = R(R(:,1)>0,:); % we probably cross the MP only for positive X
R = R(irf_abs(R,1)>7*R_E,:); % we probably cross the MP only R > 7 R_E

if isempty(R)
  irf_log('proc','tail season')
  return
end

start_time = R(1,1);
dt = R(end,1) -R(1,1);

irf_log('proc',['X>0, R>7R_E: ' epoch2iso(start_time,1) ' -- ' ...
  epoch2iso(start_time+dt,1)])

% Fetch ACE data
if ismac,  ISTP_PATH = '/Volumes/istp';
else, ISTP_PATH = '/data/istp';
end

ace_B = irf_istp_get(ISTP_PATH, start_time -120*60, dt +240*60, sc_source, 'b');
ace_V = irf_istp_get(ISTP_PATH, start_time -120*60, dt +240*60, sc_source, 'v');
ace_N = irf_istp_get(ISTP_PATH, start_time -120*60, dt +240*60, sc_source, 'n');

% Create new timeline with 30 min step
st_a = fromepoch(start_time);
st = toepoch([st_a(1:4) fix(st_a(5)/30)*30 00]);
dt = ceil((start_time +dt -st)/1800)*1800;

irf_log('proc',['subint: ' epoch2iso(start_time,1) ' -- ' ...
  epoch2iso(start_time+dt,1)])

%v_ttt = []; b_ttt = []; n_ttt = [];
r_prev = [];
for t=st:1800:st+dt
  irf_log('proc',['time: ' epoch2iso(t,1)])
  
  % ACE time shift
  if strcmp(sc_source,'ace')
    if isempty(ace_V), dt_ace = ACE_DT_DEF;
    else
      v_tmp = linear_solve(ace_V, t, ACE_DT_DEF);
      %irf_log('proc',['ace_v_tmp: ' num2str(round(v_tmp)) ' km/s'])
      if isnan(v_tmp)
        dt_ace = ACE_DT_DEF;
        irf_log('proc',['ace_v_tmp: NaN at ' epoch2iso(t,1)])
      else, dt_ace = ACE_X_POS/v_tmp;
      end
    end
    %irf_log('proc',['ace_dt   : ' num2str(round(dt_ace/60)) ' min'])
  else
    dt_ace = 0;
  end
  
  if isempty(ace_V), vx_tmp = ACE_VX_DEF;
  else
    vx_tmp = linear_solve(ace_V, t, dt_ace);
    if isnan(vx_tmp)
      irf_log('proc',['ace_vx: NaN at ' epoch2iso(t,1)])
      vx_tmp = ACE_VX_DEF;
    end
  end
  %irf_log('proc',['ace_vx_tmp: ' num2str(vx_tmp,'%.2f') ' km/s'])
  %v_ttt = [v_ttt; t-dt_ace vx_tmp];
  
  if isempty(ace_N), n_tmp = ACE_N_DEF;
  else
    n_tmp = linear_solve(ace_N, t, dt_ace);
    if isnan(n_tmp)
      irf_log('proc',['ace_n : NaN at ' epoch2iso(t,1)])
      n_tmp = ACE_N_DEF;
    end
  end
  %irf_log('proc',['ace_nn_tmp: ' num2str(n_tmp,'%.2f') ' cc'])
  %n_ttt = [n_ttt; t-dt_ace n_tmp];
  
  if isempty(ace_B), bz_tmp = ACE_BZ_DEF;
  else
    bz_tmp = linear_solve(ace_B(:,[1 4]), t, dt_ace);
    if isnan(bz_tmp)
      irf_log('proc',['ace_bz: NaN at ' epoch2iso(t,1)])
      bz_tmp = ACE_BZ_DEF;
    end
  end
  %irf_log('proc',['ace_vx_tmp: ' num2str(bz_tmp,'%.2f') ' nT'])
  %b_ttt = [b_ttt; t-dt_ace bz_tmp];
  
  r_tmp = linear_solve(R, t, 0);
  if isnan(r_tmp)
    irf_log('proc',['R : NaN at ' epoch2iso(t,1)])
    continue
  end
  
  r_gsm = irf_gse2gsm([t r_tmp]);
  r_gsm(2:4) = r_gsm(2:4)/R_E;
  r_mp = r_shue_mp(r_gsm, bz_tmp, nv2press(n_tmp,vx_tmp^2));
  %irf_log('proc',['r: ' num2str(r_gsm(2:4),'%.2f %.2f %.2f') ...
  %		' mp:' num2str(r_mp,'%.2f') ' Re'])
  
  if isempty(t_mp_out) && ~isempty(r_prev) && (r_prev>0) && (r_mp<0)
    t_mp_out = t -1800;
    irf_log('proc',['FOUND OUTBOUND : ' epoch2iso(t_mp_out,1)])
  end
  if ~isempty(r_prev) && (r_prev<0) && (r_mp>0)
    t_mp_in = t;
    irf_log('proc',['FOUND INBOUND  : ' epoch2iso(t_mp_in,1)])
  end
  r_prev = r_mp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Help functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = linear_solve(f, t, delta_t)
% help function for linear interpolation

t = t -delta_t;
p1 = f( (f(:,1) > t -3600) & (f(:,1) <= t), : ); % data has 1 hour resolutuion
p2 = f( (f(:,1) < t +3600) & (f(:,1) > t), : );
if isempty(p1) && isempty(p2)
  % We have a data gap
  y = NaN;
  return
end
if ~isempty(p1), p1 = p1(1,:); end
if ~isempty(p2)
  p2 = p2(end,:);
  if isempty(p1), y = p2(:,2:end);
  else
    a_tmp = (p1(:,2:end) -p2(:,2:end))/(p1(:,1) -p2(:,1));
    y = a_tmp*t + 0.5*(p1(:,2:end) +p2(:,2:end) ...
      -a_tmp*(p1(:,1) +p2(:,1)));
  end
else, y = p1(:,2:end);
end
return

function res = nv2press(n,v2)
%function res = nv2press(n,v2)
%
% Calculate plasma dynamic pressure in nPa
% n in 1/cc
% v^2 in [km/s]^2

n = n(:);
v2 = v2(:);
% p=nmv^2 ;-)
res = 1.6726*1e-6*v2.*n;
return

function r_mp = r_shue_mp(pos_Re_gsm, bz_nT, swp_nPa)
%IRF_SHUE_MP  estimate distance to model(Shue) magnetopause
%
%  r_mp_Re = irf_shue_mp(pos_Re_gsm, bz_nT, swp_nPa)
%
% Input:
%		pos_Re_gsm - GSM position in Re (3 or 4 components)
%		bz_nT      - IMF Bz in nT
%		swp_nPa    - Solar wind dynamic pressure in nPa
%
% References:
%		Shue et. al., A new functional form to study the solar
% 		wind control of the magnetopause size ans shape,
%		JGR, 102, p.9497, 1997.
%
% See also IRF_GSE2GSM

% Copyright 2006 Yuri Khotyaintsev

if size(pos_Re_gsm,2)>3, pos_Re_gsm = pos_Re_gsm(:,2:4); end

% Shue et. al., Eq. 13
alpha = ( 0.58 -0.01*bz_nT )*( 1.0 +0.01*swp_nPa );

% Shue et. al., Eq. 12
if bz_nT>=0, r0 = ( 11.4 +0.013*bz_nT )*swp_nPa^( -1.0/6.6 );
else,        r0 = ( 11.4 +0.140*bz_nT )*swp_nPa^( -1.0/6.6 );
end

r = irf_abs(pos_Re_gsm,1);
cosTheta = pos_Re_gsm(:,1)./r;
% Shue et. al., Eq. 1
r_mp = r0 *( 2.0./( 1.0 +cosTheta )).^alpha - r;

return

