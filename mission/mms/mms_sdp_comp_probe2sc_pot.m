function probe2sc_pot = mms_sdp_comp_probe2sc_pot(filterInterval)
%MMS_SDP_COMP_PROBE2SC_POT  compute PROBE2SC_POT for datamanager
%
%  probe2sc_pot = mms_sdp_comp_probe2sc_pot(filterInterval)
%
%  Compute probe-to-spacecraft potential averaged from several probes using
%  the mean value of moving average filtered data over filterInterval (in
%  seconds) to determine which probe(-s) are possibly bad. For each
%  timestamp either all four, two or one probe(-s) are used.

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% Default to 20 seconds interval (one nominal spin period), if not specified.
if nargin==0,  filterInterval = 20; end

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
probe2sc_pot = MMS_CONST.Error;

procId = mms_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl,MMS_CONST.SDCProc.ql,...
      MMS_CONST.SDCProc.l2pre}
    dcv = mms_sdp_datamanager('dcv');
    if mms_is_error(dcv)
      errStr='Bad DCV input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    
    % FIXME: see what signals do we actually have

    samplerate = mms_sdp_datamanager('samplerate');
    if mms_is_error(samplerate)
      errStr='Bad SAMPLERATE input, cannot proceed.';
      irf.log('critical',errStr); error(errStr);
    end
    % Filter window size, default 20 s * Samplerate = 160 samples (slow),
    % 640 samples (fast), 163'840 samples (brst).
    windowSize = samplerate*filterInterval;
    % Create filter coefficients for moving average filter.
    a = 1; b = (1/windowSize)*ones(1,windowSize);
    % Apply moving average filter (a,b) on spin plane probes 1, 2, 3 & 4.
    MAfilt = filter(b, a, [dcv.v1.data, dcv.v2.data, dcv.v3.data, dcv.v4.data], [], 1);
    % For each timestamp get median value of the moving average.
    MAmedian = median(MAfilt, 2);
    % For each probe check if it is too far off from the median
    absDiff = abs(MAfilt - repmat(MAmedian, [1 4]));
    badBits = absDiff > MMS_CONST.Limit.DIFF_PROBE_TO_SCPOT_MEDIAN;

    % DEBUG: FORCE SOME bad bits
    %badBits(1:25,:)=round(rand(25,4));
    %badBits(20,:)=1;
    %badBits(22,:)=1;

    % Identify times with all four probes marked as bad
    ind_row = ismember(badBits, [1 1 1 1], 'rows');
    if( any(ind_row))
      irf.log('warning', 'Some timestamps show all four probes as outliers. Using the "best" (closest to median) probe.');
      %absDiff = abs(MAfilt(ind_row,:)-repmat(MAmedian(ind_row), [1 4]));
      minAbsDiff = min(absDiff(ind_row,:), [], 2);
      badBitsSeg = absDiff(ind_row,:) > repmat(minAbsDiff, [1 4]);
      badBits(ind_row,:) = badBitsSeg;
    end

    % Identify times with three bad probes
%    ind_row = ismember(badBits, [1 1 1 0; 1 1 0 1; 1 0 1 1; 0 1 1 1], 'rows');
%    if( any(ind_row))
%      irf.log('debug', 'Some timestamps show three probes as outliers. Using only the "best" probe.');
%    end

    % Identify times with two bad probes
%    ind_row = ismember(badBits, [1 0 1 0; 0 1 0 1], 'rows');
%    if( any(ind_row))
%      irf.log('debug','Some timestamps show two (pair 12 or pair 34) as outliers. Using the other pair.');
%    end
    ind_row = ismember(badBits, [1 0 0 1; 0 1 1 0], 'rows');
    if( any(ind_row))
      irf.log('warning','Some timestamps show two (not pairwise) probes as outliers. Using them anyhow...');
      % FIXME: WHAT TO DO? Use the other two? Or only the "Best"?
    end

    % Identify times with one bad probe, set the entire pair as bad.
    ind_row = ismember(badBits, [1 0 0 0 ; 0 1 0 0], 'rows');
    if( any(ind_row))
      irf.log('notice','Some timestamps show one probe as outlier, removing this probe pair (12) for those times.');
      badBits(ind_row, 1:2) = 1;
    end
    ind_row = ismember(badBits, [0 0 1 0 ; 0 0 0 1], 'rows');
    if( any(ind_row))
      irf.log('notice','Some timestamps show one probe as outlier, removing this probe pair (34) for those times.');
      badBits(ind_row, 3:4) = 1;
    end

    % Set all bad bits to NaN in data before calculating the averaged value
    dcv.v1.data(badBits(:,1)) = NaN;
    dcv.v2.data(badBits(:,2)) = NaN;
    dcv.v3.data(badBits(:,3)) = NaN;
    dcv.v4.data(badBits(:,4)) = NaN;

    % Compute average of all spin plane probes, ignoring data identified as
    % bad (NaN).
    avPot = irf.nanmean([dcv.v1.data, dcv.v2.data, dcv.v3.data, dcv.v4.data], 2);

    % Compute bitmask of which probes was used so that:
    % bitmask == 0 all probes useful.
    % bitand(0x01, bitmask) probe 1 bad, bitand(0x02, bitmask) probe 2 bad,
    % bitand(0x04, bitmask) probe 3 bad, bitand(0x08, bitmask) probe 4 bad.
    % And any combination thereof, ie. bitand(0x03, bitmask) probe 1 & 2
    % bad and probe 3 & 4 was used.
    bitmask = 1*badBits(:,1) + 2*badBits(:,2) + 4*badBits(:,3) + 8*badBits(:,4);

    probe2sc_pot = struct('time',dcv.time,'data',avPot,'bitmask',bitmask);
    
  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
