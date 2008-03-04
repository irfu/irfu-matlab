function caa_update_deltaoff(d,cl_id)
%CAA_UPDATE_DELTAOFF  update delta offsets database
%
% CAA_UPDATE_DELTAOFF(DELTAOFF,CL_ID)
%
% Update database delta offsets for CL_ID by collted data DELTAOFF.
% The database is stored in caa/deltaoff.mat
%
% See also CAA_COROF_DELTA, ClusterProc/getData
% 
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

DT = 86400*6;
SY = 2001;
EY = 2006;

d = real(d) - imag(d); % Old style delta offsets (imag is applied to p34)
ndata = fix((d(end,1)-d(1,1))/DT);
t = d(1,1):DT:( d(1,1) + DT*ndata ) + DT/2;
da = irf_resamp(d,t','thresh',.5); %#ok<NASGU>

figure
h = irf_plot({d,da},'comp','linestyle',{'-','*'});
ylabel(h(1),'\Delta Ex [mV/m]')
ylabel(h(2),'\Delta Ey [mV/m]')
title(h(1),sprintf('Cluster %d',cl_id))

ud = get(gcf,'userdata');
t_st_e = double(ud.t_start_epoch);
tt = zeros(1,EY-SY+1);
ts = cell(size(tt));
for i=1:length(tt)
	tt(i) = toepoch([SY+i-1 1 1 0 0 0]) - t_st_e;
	ts{i} = num2str(SY+i-1);
end
set(h(1),'XTick',tt,'XLim',[tt(1) tt(end)],'Ylim',[-2 2])
set(h(2),'XTick',tt,'XLim',[tt(1) tt(end)],'Ylim',[-2 2],'XTickLabel',ts)

y = irf_ask('Save? yes/no [%]>','y','no');
if strcmpi(y,'y') || strcmpi(y,'yes')
	c_eval('DeltaOff?=da; save deltaoff DeltaOff? -append; disp(''DeltaOff? -> deltaoff.mat'')',cl_id)
end