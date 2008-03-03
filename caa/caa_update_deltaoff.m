function caa_update_deltaoff(d,cl_id)
%CAA_UPDATE_DELTAOFF  update delta offsets database
%
% caa_update_deltaoff(DeltaOff,cl_id)
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

d = real(d) - imag(d); % Old style delta offsets (imag is applied to p34)
ndata = fix((d(end,1)-d(1,1))/DT);
t = d(1,1):DT:( d(1,1) + DT*ndata ) + DT/2;
da = irf_resamp(d,t','thresh',.5); %#ok<NASGU>
c_eval('DeltaOff?=da; save deltaoff DeltaOff? -append; disp(''DeltaOff? -> deltaoff.mat'')',cl_id)
figure
irf_plot({d,da},'comp')
