function data = caa_corof_delta(data,probe_p,dof,action)
%CAA_COROF_DELTA  correct delta offsets
%
% NEW_DATA = CAA_COROF_DELTA(DATA,PROBE_P,DELTAOFF,[ACTION])
%
% Correct/remove delta (p12 vs p34) offset DELTAOFF on DATA from probe
% pair PROBE_P.
%
% ACTION:
%     'apply' - apply the offset (default)
%     'undo' - remove already applied offset
%
% See also C_EFW_DELTA_OFF
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if probe_p~=12 && probe_p~=32 && probe_p~=34
  error('bad value for PROBE_P')
end

if probe_p==34 && isreal(dof), return, end
if (probe_p == 12 || probe_p == 32) && ~isreal(dof), return, end

if nargin<4, action = 1;
else
  if ischar(action)
    if strcmp(action,'apply'), action = 1;
    elseif strcmp(action,'undo'), action = -1;
    end
  elseif action~=1 || action~=-1
    error('bad value for ACTION')
  end
end

% Real offset is applied to p12/32,  imaginary to p34.
for comp=1:2
  if (isreal(dof(comp)) && probe_p~=34) || (~isreal(dof(comp)) && probe_p==34)
    if ~isreal(dof(comp)), dof(comp) = imag(dof(comp)); end
    data(:,comp+1) = data(:,comp+1) - action*dof(comp);
    if action==1, do='applying'; else, do='removing'; end
    if comp==1, x='x'; else, x='y'; end
    irf_log('proc',sprintf('%s %.2f mV/m delta E%s on p%d',do,dof(comp),x,probe_p))
  end
end


