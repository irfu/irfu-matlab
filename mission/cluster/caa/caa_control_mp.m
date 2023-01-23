function caa_control_mp(st,et)
%CAA_CONTROL_MP  plot predicted MP location
%
% caa_control_mp(st,et)
%
% Plots predicted magnetopause
%
% See also CAA_FIND_MP, CAA_SH_PLAN
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if ischar(st), st = iso2epoch(st); end
if ischar(et), et = iso2epoch(et); end

st_a = fromepoch(st);

v_s = sprintf('ORB1Y%d',st_a(1));

eval(['load mPlan ' v_s])

ORB = [];
if exist(v_s,'var'), eval(['ORB=' v_s ';'])
else
  disp(['cannot load ' v_s])
  return
end

ORB = ORB(ORB(:,1)>=st & ORB(:,1)<=et,:);

v_s = sprintf('MP?Y%d',st_a(1));
c_eval(['load mPlan ' v_s])

for cl_id = 1:4
  if exist(irf_ssub(v_s,cl_id),'var'), c_eval(['MP?=' v_s ';'])
  else
    disp(['cannot load ' v_s])
    return
  end
end

OFF = 2*3600;

for o = 1:length(ORB)
  mp_in = 0; mp_out = 0; cnt = 0;
  c_eval('ii=find(MP?(:,1)>=ORB(o,1) & MP?(:,1)<=ORB(o,1)+ORB(o,2)); if ~isempty(ii), cnt=cnt+1; mp_cur?=MP?(ii,:); mp_out=mp_out+MP?(ii,1); mp_in=mp_in+MP?(ii,2); end')
  if cnt==0
    irf_log('proc',['No MP for orbit: ' epoch2iso(ORB(o,1),1)...
      ' -- ' epoch2iso(ORB(o,1)+ORB(o,2),1)])
    continue
  end
  mp_in = mp_in/cnt +1800; %XXX temporary fix
  mp_out = mp_out/cnt;
  irf_log('proc',['Mp OUT: ' epoch2iso(mp_out,1)...
    ' IN:' epoch2iso(mp_in,1)])
  
  figure(162), clf
  set(gcf,'Position',[520   224   787   876])
  pp = 0;
  for cl_id = 1:4
    subplot(4,2,cl_id*2-1)
    p1 = get(gca, 'Position');
    if ~pp, pp=p1(4)*1.2; end
    set(gca, 'Position', [p1(1:3) pp]);
    irf_log('proc',['C' num2str(cl_id) ' ' epoch2iso(mp_out-OFF,1)...
      ' -- ' epoch2iso(mp_out+OFF,1)])
    Ps = caa_get(mp_out-OFF,2*OFF,cl_id,'Ps?');
    if ~isempty(Ps), irf_plot(Ps), end
    irf_zoom(mp_out +[-OFF OFF],'x',gca)
    c_eval('mp_cur = mp_cur?;', cl_id)
    yl = get(gca,'YLim');
    hold on
    irf_plot([mp_cur(1) mp_cur(1); yl]','ko-')
    hold off
    ylabel(irf_ssub('Cluster ?',cl_id))
    if cl_id==1, title('OUT'), end
    if cl_id~=4, xlabel(''), set(gca,'XTickLabel',[]), end
    
    
    subplot(4,2,cl_id*2)
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p1(2) p(3) pp]);
    irf_log('proc',['C' num2str(cl_id) ' ' epoch2iso(mp_in-OFF,1)...
      ' -- ' epoch2iso(mp_in+OFF,1)])
    Ps = caa_get(mp_in-OFF,2*OFF,cl_id,'Ps?');
    if ~isempty(Ps), irf_plot(Ps), end
    irf_zoom(mp_in +[-OFF OFF],'x',gca)
    c_eval('mp_cur = mp_cur?;', cl_id)
    yl = get(gca,'YLim');
    hold on
    irf_plot([mp_cur(2)+1800 mp_cur(2)+1800; yl]','ko-') %XXX temporary fix
    hold off
    ylabel(irf_ssub('Cluster ?',cl_id))
    if cl_id==1, title('IN'), end
    if cl_id~=4, xlabel(''), set(gca,'XTickLabel',[]), end
  end
  orient tall
  fn = sprintf('C_MP_%s',irf_fname(ORB(o,1)));
  irf_log('save',['saving ' fn])
  print( gcf, '-dpdf', fn)
end
