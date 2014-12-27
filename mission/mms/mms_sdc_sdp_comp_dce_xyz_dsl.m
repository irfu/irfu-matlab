function dce_xyz_dsl = mms_sdc_sdp_comp_dce_xyz_dsl
%MMS_SDC_SDP_COMP_DCE_XYZ_DSL  compute DCE_XYZ_DSL for datamanager
%
%  probe2sc_pot = mms_sdc_sdp_comp_dce_xyz_dsl()
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

procId = mms_sdc_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl,MMS_CONST.SDCProc.ql,...
      MMS_CONST.SDCProc.l2pre}
    hk_101 = mms_sdc_sdp_datamanager('hk_101');
    if isnumeric(hk_101) && numel(hk_101)==1 && hk_101==MMS_CONST.Error,
      irf.log('warning','Bad hk_101 input'); return
    end
    dce = mms_sdc_sdp_datamanager('dce');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad dce input'); return
    end
    phase = mms_sdc_sdp_datamanager('phase');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad phase input'); return
    end
    
    % FIXME: add ADC offsets here
    
    dE = mms_sdp_despin(dce.e12.data,dce.e34.data,phase.data);
    
    % FIXME: need to compute from respective bitmasks
    bitmask = dce.e12.bitmask;
    
    % FIXME: apply DSL offsets here 

    dce_xyz_dsl = struct('time',dce.time,'data',[dE dce.e56.data],...
      'bitmask',bitmask);
  case MMS_CONST.Error
    errStr = 'mms_sdc_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
