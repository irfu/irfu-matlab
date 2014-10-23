function phase = mms_sdc_sdp_comp_phase
%MMS_SDC_SDP_COMP_PHASE  compute phase for datamanager
%
%  phase = mms_sdc_sdp_comp_phase()

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
phase = MMS_CONST.Error;

procId = mms_sdc_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql}
    hk_101 = mms_sdc_sdp_datamanager('hk_101');
    if isnumeric(hk_101) && numel(hk_101)==1 && hk_101==MMS_CONST.Error,
      irf.log('warning','Bad hk_101 input'); return
    end
    dce = mms_sdc_sdp_datamanager('dce');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad dce input'); return
    end
    [dcephase, dcephase_flag] = mms_sdc_sdp_phase_2(hk_101, dce.time);
    phase = struct('data',dcephase,'bitmask',dcephase_flag);
  case MMS_CONST.SDCProc.l2pre
    %XXX: this needs to be replaced with definitive attitude!!!
    hk_101 = mms_sdc_sdp_datamanager('hk_101');
    if isnumeric(hk_101) && numel(hk_101)==1 && hk_101==MMS_CONST.Error,
      irf.log('warning','Bad hk_101 input'); return
    end
    dce = mms_sdc_sdp_datamanager('dce');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad dce input'); return
    end
    [dcephase, dcephase_flag] = mms_sdc_sdp_phase_2(hk_101, dce.time);
    phase = struct('data',dcephase,'bitmask',dcephase_flag);
  case MMS_CONST.Error
    errStr = 'mms_sdc_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end
