function probe2sc_pot = mms_sdp_comp_probe2sc_pot()
%MMS_SDP_COMP_PROBE2SC_POT  compute PROBE2SC_POT for datamanager
%
%  probe2sc_pot = mms_sdp_comp_probe2sc_pot()
%
%  Compute probe-to-spacecraft potential averaged from several probes using
%  the mean value of average filtered data over filterInterval (in
%  seconds) to determine which probe(-s) are possibly bad. For each
%  timestamp either all four, two or one probe(-s) are used.
%  Each datapoint is given a corresponding bitmask where
%   bit 0 = 0, only one probe or no probe at all was used. If no probe at
%              all was available the output is NaN for that point in time.
%   bit 0 = 1, either two or four probes was used.
%   bits 1-16, are a bitor comination of the corresponding bitmasks of the
%              individual probes used for that point in time.

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
probe2sc_pot = MMS_CONST.Error; %#ok<NASGU>

procId = mms_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl,MMS_CONST.SDCProc.ql,...
      MMS_CONST.SDCProc.l2pre}
    Dcv = mms_sdp_datamanager('dcv');
    if mms_is_error(Dcv)
      errStr='Bad DCV input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    
    % FIXME: see what signals do we actually have

    sampleRate = mms_sdp_datamanager('samplerate');
    if mms_is_error(sampleRate)
      errStr='Bad SAMPLERATE input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    
    probe2sc_pot = mms_sdp_dmgr.comp_probe2sc_pot(Dcv,sampleRate,MMS_CONST);
    
  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
