function probe2sc_pot = mms_sdc_sdp_comp_probe2sc_pot
%MMS_SDC_SDP_COMP_PROBE2SC_POT  compute PROBE2SC_POT for datamanager
%
%  probe2sc_pot = mms_sdc_sdp_comp_probe2sc_pot()
%
%  Compute probe-to-spacecraft potential averaged from several probes

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
probe2sc_pot = MMS_CONST.Error;

procId = mms_sdc_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl,MMS_CONST.SDCProc.ql,...
      MMS_CONST.SDCProc.l2pre}
    dcv = mms_sdc_sdp_datamanager('dcv');
    if isnumeric(dcv) && numel(dcv)==1 && dcv==MMS_CONST.Error,
      irf.log('warning','Bad DCV input'); return
    end
    
    % FIXME: see what signals do we acrually have
    % Average the spin plane probes
    avPot = 0.25 * (dcv.v1.data + dcv.v2.data + dcv.v3.data + dcv.v4.data);
    
    % FIXME: need to compute from respective bitmasks
    bitmask = dcv.v1.bitmask;

    probe2sc_pot = struct('time',dcv.time,'data',avPot,'bitmask',bitmask);
    
  case MMS_CONST.Error
    errStr = 'mms_sdc_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
