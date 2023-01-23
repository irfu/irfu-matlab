function caa_ms_pl_xoff(st,et)
%CAA_MS_PL_XOFF  visualize offset in the magnetosphere
%
% caa_ms_pl_xoff(st,et)
%
% See also CAA_SH_PL_XOFF
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

[stepoch,dt] = irf_stdt(st,et);

c_eval('diEs?p34=caa_get(st,et,?,''diEs?p34'');')
c_eval('Ps?=caa_get(st,et,?,''Ps?'');')
c_eval('diEDI?=caa_get(st,et,?,''diEDI?'');')
c_eval('diVCEh?=caa_get(st,et,?,''diVCEh?'');')
R1 = caa_get(st,et,1,'R?');

t = stepoch:600:(stepoch+dt); t = t'; %#ok<NASGU>
c_eval('diEs?r=irf_resamp(diEs?p34(~isnan(diEs?p34(:,2)),1:2),t);')
c_eval('if isempty(diEDI?), diEDI?r=[]; else diEDI?r=irf_resamp(diEDI?(:,1:2),t); end')
c_eval('if isempty(diVCEh?), diVCEh?r=[]; else diVCEh?r=irf_resamp(diVCEh?(:,1:2),t); end')

clf
h(1) = irf_subplot(2,1,-1);
c_pl_tx('diEs?r',2,'.')
hold on
c_pl_tx('diEDI?r',2,'x')
c_pl_tx('diVCEh?r',2,'+')

h(2) = irf_subplot(2,1,-2);
c_pl_tx('Ps?')

irf_zoom(stepoch + [0 dt],'x',h)
add_position(h(2),R1)
