function dce_xyz_dsl = mms_sdp_comp_dce_xyz_dsl
%MMS_SDP_COMP_DCE_XYZ_DSL  compute DCE_XYZ_DSL for datamanager
%
%  probe2sc_pot = mms_sdp_comp_dce_xyz_dsl()
%
%  Compute the electric field in DSL (despun)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
dce_xyz_dsl = MMS_CONST.Error;

procId = mms_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl,MMS_CONST.SDCProc.ql,...
      MMS_CONST.SDCProc.l2pre} 
    dce = mms_sdp_datamanager('dce');
    if mms_is_error(dce), irf.log('warning','Bad dce input'); return, end
    phase = mms_sdp_datamanager('phase');
    if mms_is_error(phase), irf.log('warning','Bad phase input'); return, end

    % FIXME: add ADC offsets here

    dE = mms_sdp_despin(dce.e12.data,dce.e34.data,phase.data);

    % FIXME: need to compute from respective bitmasks
    bitmask = dce.e12.bitmask;

    % FIXME: apply DSL offsets here 

    dce_xyz_dsl = struct('time',dce.time,'data',[dE dce.e56.data],...
      'bitmask',bitmask);

  case MMS_CONST.SDCProc.l2a
    % ADC offsets should have already been applied, remainging processing
    % full despin (from spinfits) and DSL offset.
    dce = mms_sdp_datamanager('dce');
    if mms_is_error(dce)
      irf.log('warning','Bad dce input'); return
    end
    spinfits = mms_sdp_datamanager('spinfits');
    if mms_is_error(spinfits)
      irf.log('warning','Bad spinfits input'); return
    end

    irf.log('critical', 'DCE_XYZ_DSL calculation for L2A not performed using spinfits from L2pre. FIXME!');
    bitmask = dce.e12.bitmask;
    % FIXME: apply DSL offsets here
    dce_xyz_dsl = struct('time',dce.time,'data',[dce.e12.data, dce.e34.data, dce.e56.data],...
      'bitmask',bitmask);

  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr);

  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr);
end

end
