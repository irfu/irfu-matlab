function sc_pot = mms_sdp_comp_sc_pot
%MMS_SDP_COMP_SC_POT  compute SC_POT for datamanager
%
%  sc_pot = mms_sdp_comp_sc_pot()
%
%  Compute estimate of the spacecraft potential

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
sc_pot = MMS_CONST.Error;

procId = mms_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl,MMS_CONST.SDCProc.ql,...
      MMS_CONST.SDCProc.l2pre}
    probe2sc_pot = mms_sdp_datamanager('probe2sc_pot');
    if isnumeric(probe2sc_pot) && numel(probe2sc_pot)==1 &&...
        probe2sc_pot==MMS_CONST.Error,
      irf.log('warning','Bad PROBE2SC_POT input'); return
    end
    
    % XXX: add a better estimate of the plasma potential
    plasmaPotential = 1;
    
    % XXX: add a better estimate for the shortening factor
    shorteningFactor = 1.1;
    
    % FIXME: need to compute from respective bitmasks
    bitmask = probe2sc_pot.bitmask;
    
    scPot = -probe2sc_pot.data*shorteningFactor + plasmaPotential;

    sc_pot = struct('time',probe2sc_pot.time,'data',scPot,'bitmask',bitmask);
    
  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
