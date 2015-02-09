function probe2sc_pot = mms_sdc_sdp_comp_probe2sc_pot(filterInterval)
%MMS_SDC_SDP_COMP_PROBE2SC_POT  compute PROBE2SC_POT for datamanager
%
%  probe2sc_pot = mms_sdc_sdp_comp_probe2sc_pot(filterInterval)
%
%  Compute probe-to-spacecraft potential averaged from several probes using
%  the mean value of moving average filtered data over filterInterval (in
%  seconds) to determine which probe(-s) are possibly bad.

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% Default to 10 seconds interval, if not specified.
if nargin==0,  filterInterval = 20; end

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
    
    % FIXME: see what signals do we actually have

    tmMode = mms_sdc_sdp_datamanager('tmMode');
    switch tmMode
      case {MMS_CONST.TmMode.slow, MMS_CONST.TmMode.fast, MMS_CONST.TmMode.brst}
        % Filter window size, default 20 s * Samplerate = 160 samples (slow),
        % 640 samples (fast), 163'840 samples (brst).
        windowSize = MMS_CONST.Samplerate.(MMS_CONST.TmModes{tmMode})*filterInterval;
      otherwise
        errStr = 'Unrecognized tmMode';
        irf.log('critical', errStr); error(errStr);
    end
    % Create filter coefficients for moving average filter.
    a = 1; b = (1/windowSize)*ones(1,windowSize);
    % Apply moving average filter (a,b) on spin plane probes 1, 2, 3 & 4.
    MAfilt = filter(b, a, [dcv.v1.data, dcv.v2.data, dcv.v3.data, dcv.v4.data], [], 1);
    % For each timestamp get median value of the moving average.
    MAmedian = median(MAfilt, 2);
    % For each probe check if it is too far off from the median
    for ii=1:4
      % "Bad" filtered values are outside median by absolute value specified in MMS_CONST
      bitsBad.p{ii} = abs(MAfilt(:,ii) - MAmedian) > MMS_CONST.Limit.DIFF_PROBE_TO_SCPOT_MEDIAN;
      if(any(bitsBad.p{ii}>0))
        logStr=sprintf('Removing %i outliers from dcv probe %i before averaging probes to one SC potential.',numel(find(bitsBad.p{ii}>0)), ii);
        irf.log('notice',logStr);
        % Replace "bad" data with NaN
        dcv.(['v', num2str(ii)]).data(bitsBad.p{ii}) = NaN;
      end
    end
    % Compute average of all spin plane probes, ignoring data identified as
    % bad (NaN).
    avPot = irf.nanmean([dcv.v1.data, dcv.v2.data, dcv.v3.data, dcv.v4.data], 2);

    % Compute bitmask of which probes was used so that:
    % bitmask == 0 all probes useful.
    % bitand(0x01, bitmask) probe 1 bad, bitand(0x02, bitmask) probe 2 bad,
    % bitand(0x04, bitmask) probe 3 bad, bitand(0x08, bitmask) probe 4 bad.
    % And any combination thereof, ie. bitand(0x03, bitmask) probe 1 & 2
    % bad and probe 3 & 4 was used.
    bitmask = 1*bitsBad.p{1} + 2*bitsBad.p{2} + 4*bitsBad.p{3} + 8*bitsBad.p{4};

    probe2sc_pot = struct('time',dcv.time,'data',avPot,'bitmask',bitmask);
    
  case MMS_CONST.Error
    errStr = 'mms_sdc_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
