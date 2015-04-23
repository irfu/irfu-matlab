function probe2sc_pot = mms_sdp_comp_probe2sc_pot(filterInterval)
%MMS_SDP_COMP_PROBE2SC_POT  compute PROBE2SC_POT for datamanager
%
%  probe2sc_pot = mms_sdp_comp_probe2sc_pot(filterInterval)
%
%  Compute probe-to-spacecraft potential averaged from several probes using
%  the mean value of moving average filtered data over filterInterval (in
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

% Default to 20 seconds interval (one nominal spin period), if not specified.
if nargin==0,  filterInterval = 20; end

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
probe2sc_pot = MMS_CONST.Error; %#ok<NASGU>

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
    
    % Blank sweeps
    sweepBit = MMS_CONST.Bitmask.SWEEP_DATA;
    dcv.v1.data = mask_bits(dcv.v1.data, dcv.v1.bitmask, sweepBit);
    dcv.v2.data = mask_bits(dcv.v2.data, dcv.v2.bitmask, sweepBit);
    dcv.v3.data = mask_bits(dcv.v3.data, dcv.v3.bitmask, sweepBit);
    dcv.v4.data = mask_bits(dcv.v4.data, dcv.v4.bitmask, sweepBit);
    
    % Compute average of all spin plane probes, ignoring data identified as
    % bad (NaN).
    avPot = irf.nanmean([dcv.v1.data, dcv.v2.data, dcv.v3.data, dcv.v4.data], 2);

    % Combine bitmask so that bit 0 = 0 (either four or two probes was
    % used), bit 0 = 1 (either one probe or no probe (if no probe => NaN
    % output in data). The other bits are a bitor comination of those
    % probes that were used (i.e. bitmask = 2 (dec), would mean at least
    % one probe that was used for that point in time had "bad bias").
    % Start with bit 0
    bitmask = uint16(sum(badBits,2)>=3); % Three or more badBits on each row.
    % Extract probe bitmask, excluding the lowest bit (signal off)
    bits = intmax(class(bitmask)) - MMS_CONST.Bitmask.SIGNAL_OFF;
    vBit = zeros(length(dcv.v1.bitmask),4,'like',MMS_CONST.Bitmask.SIGNAL_OFF);
    vBit(:,1) = bitand(dcv.v1.bitmask, bits);
    vBit(:,2) = bitand(dcv.v2.bitmask, bits);
    vBit(:,3) = bitand(dcv.v3.bitmask, bits);
    vBit(:,4) = bitand(dcv.v4.bitmask, bits);
    % Update badBits due to blanking of sweep.
    badBits(:,1) = isnan(dcv.v1.data); badBits(:,2) = isnan(dcv.v2.data);
    badBits(:,3) = isnan(dcv.v3.data); badBits(:,4) = isnan(dcv.v4.data);
    % Combine bitmasks with bitor of times when probe was used to derive
    % mean. (I.e. not marked by badBits).
    for ii=1:4
      bitmask(~badBits(:,ii)) = bitor(bitmask(~badBits(:,ii)), vBit(~badBits(:,ii),ii));
    end

    probe2sc_pot = struct('time',dcv.time,'data',avPot,'bitmask',bitmask);
    
  case MMS_CONST.Error
    errStr = 'mms_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

end
