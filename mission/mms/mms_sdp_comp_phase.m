function phase = mms_sdp_comp_phase
%MMS_SDP_COMP_PHASE  compute phase for datamanager
%
%  phase = mms_sdp_comp_phase()

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
phase = MMS_CONST.Error;

procId = mms_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql}
    hk_101 = mms_sdp_datamanager('hk_101');
    if mms_is_error(hk_101)
      errStr='Bad HK_101 input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    dce = mms_sdp_datamanager('dce');
    if mms_is_error(dce)
      errStr='Bad DCE input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    [dcephase, dcephase_flag] = mms_sdp_phase_2(hk_101, dce.time);
    phase = struct('data',dcephase,'bitmask',dcephase_flag);
    
  case {MMS_CONST.SDCProc.l2pre,MMS_CONST.SDCProc.l2a}   
    dce = mms_sdp_datamanager('dce');
    if mms_is_error(dce)
      errStr='Bad DCE input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    defatt = mms_sdp_datamanager('defatt');
    if mms_is_error(defatt)
      errStr='Bad DEFATT input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    
    phaseTS = mms_defatt_phase(defatt,dce.time);
    dcephase_flag = zeros(size(phaseTS.data)); % FIXME BETTER FLAG & BITMASKING!
    phase = struct('data', phaseTS.data, 'bitmask', dcephase_flag);

  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end
