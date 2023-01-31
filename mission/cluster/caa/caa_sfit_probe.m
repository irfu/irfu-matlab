function [probeNum,flag_lx,probeStr] = caa_sfit_probe(cl_id,probe)
%CAA_SFIT_PROBE  set/get probe pair for spin resolution data
%
% caa_sfit_probe(cl_id,probe-p)
%    set probe pair for CL_ID (saves to mInfo.mat)
%    PROBE-P is one of 12, 32, 34 (HX) or 420 (42 LX)
%
% [p_res,flag_lx] = caa_sfit_probe(cl_id)
%    get/display probe pair for CL_ID
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2)
if cl_id<1 || cl_id>4, error('CL_ID must be 1..4'), end

if nargin>1
  if probe~=12 && probe~=32 && probe~=34 && probe~=42 && ...
      probe~=120 && probe~=320 && probe~=340 && probe~=420
    error('PROBE must be 12(0), 32(0), 34(0) or 42(0) ')
  end
  c_eval('sfit_probe?=probe;',cl_id)
  if exist('./mInfo.mat','file'), c_eval('save mInfo sfit_probe? -append',cl_id)
  else, c_eval('save mInfo sfit_probe?',cl_id)
  end
  irf_log('save',irf_ssub('sfit_probe? -> mInfo.mat',cl_id))
else
  if exist('./mInfo.mat','file')
    warning off %#ok<WNOFF>
    c_eval('load mInfo sfit_probe?;',cl_id)
    warning on %#ok<WNON>
  end
  if exist(irf_ssub('sfit_probe?',cl_id),'var')
    c_eval('pp=sfit_probe?;',cl_id)
  else
    pp = c_ctl(cl_id,'probe_p');
  end
  
  if nargout>0
    if pp>100
      flag_lx = 1;
      probeNum = pp;
      probeStr = sprintf('%dLX',probeNum/10);
    else
      flag_lx = 0;
      probeNum = pp;
      probeStr = num2str(pp);
    end
  else
    if pp>100, sLX = '(LX)'; pp = pp/10; else, sLX = ''; end
    fprintf('C%d sfit probe-pair : %d %s\n',cl_id,pp,sLX)
  end
end
