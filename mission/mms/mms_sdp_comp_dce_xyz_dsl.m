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
dce_xyz_dsl = MMS_CONST.Error; %#ok<NASGU>

procId = mms_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot, MMS_CONST.SDCProc.sitl, ...
    MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2a}

    dce = mms_sdp_datamanager('dce');
    if mms_is_error(dce)
      errStr='Bad DCE input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    phase = mms_sdp_datamanager('phase');
    if mms_is_error(phase)
      errStr='Bad PHASE intput, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    adc_off = mms_sdp_datamanager('adc_off');
    if mms_is_error(adc_off)
      errStr='Bad ADC_OFF intput, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end

    sdpProbes = fieldnames(adc_off); % default {'e12', 'e34'}
    Etmp = struct('e12',dce.e12.data,'e34',dce.e34.data);
    for iProbe=1:numel(sdpProbes)
      % Remove ADC offset
      Etmp.(sdpProbes{iProbe}) = ...
        Etmp.(sdpProbes{iProbe}) - adc_off.(sdpProbes{iProbe});
    end
    
    bitmask = uint16(bitor(dce.e12.bitmask,dce.e34.bitmask));
    Etmp.e12 = mask_bits(Etmp.e12, bitmask, MMS_CONST.Bitmask.SWEEP_DATA);
    Etmp.e34 = mask_bits(Etmp.e34, bitmask, MMS_CONST.Bitmask.SWEEP_DATA);

    dE = mms_sdp_despin(Etmp.e12, Etmp.e34, phase.data);

    % FIXME: apply DSL offsets here 

    dce_xyz_dsl = struct('time',dce.time,'data',[dE dce.e56.data],...
      'bitmask',bitmask);

  case MMS_CONST.SDCProc.l2pre
    % L2Pre is a special case should not apply despin and should not output DSL (despun),
    % but output only DCE (in instrument frame).
    errStr='DCE_XYZ_DSL should not be calculated for L2Pre, it should only have spinfits and ADC offset calculated. Should not be here.';
    irf.log('critical', errStr); error(errStr);

  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr);

  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr);
end

end
