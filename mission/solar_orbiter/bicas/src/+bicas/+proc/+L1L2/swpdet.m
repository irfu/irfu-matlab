%
% Collection of functions related to detecting bias sweeps.
%
% NOTE: SBDA and SCDA are intended to be temporary functionality while waiting
% for the long-term solution which is to use the relevant quality bit in L1/L1R
% ZV QUALITY_BITMASK for detecting sweeps. Since it is a (hopefully) short-term
% solution, the functionality is also not as sophisticated (and presumably
% accurate) or configurable as it could be.
%
% NOTE: 2024-11-27: L1 QUALITY_BITMASK in new sample datasets (not yet in
% production) contains bits for sweeps but without adding time margins before
% and after. LESIA has not yet used in "production" but are planned
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef swpdet
  % PROPOSAL: More automatic test code.
  %   PROPOSAL: For separate algorithm functions instead(?) of for consolidated
  %             functions which use multiple algorithms.
  %
  % PROPOSAL: Separate function(s) for detecting sweeps using L1R
  %           QUALITY_BITMASK. Return value in science time.
  % PROPOSAL: Separate time margins t_before, t_after.

  % PROPOSAL: SBDA_wo_margins() should not use BSO as argument.
  %   PRO: Only uses one setting.



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant)

    % The only BDM which sweeps use, but non-sweep data might use this BDM
    % too, i.e.
    % BDM<>BDM_SWEEP_POSSIBLE ==> Not a sweep.
    % BDM==BDM_SWEEP_POSSIBLE <== Sweep
    BDM_SWEEP_POSSIBLE = 4;

  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Detect sweeps using SBDA and without adding time margins.
    %
    % ALGORITHM
    % =========
    % For records before timestamp PROCESSING.L2.SWEEP_DETECTION.SBDA_SCDA_BOUNDARY_UTC:
    % BDN=4 <=> sweep
    %
    function isSweepingSbda = SBDA_wo_margins(tt2000, bdmFpa, Bso)
      sbdaEndTt2000  = spdfcomputett2000(Bso.get_fv('PROCESSING.L2.SWEEP_DETECTION.SBDA_SCDA_BOUNDARY_UTC'));

      irf.assert.sizes(...
        tt2000, [-1, 1], ...
        bdmFpa, [-1, 1]);

      bdm            = bdmFpa.int2doubleNan();
      bSbdaApplies   = tt2000 <= sbdaEndTt2000;

      isSweepingSbda = bSbdaApplies & (bdm == bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE);
    end



    % Detect sweeps using SCDA and without adding time margins.
    %
    % ALGORITHM
    % =========
    % For records after timestamp PROCESSING.L2.SWEEP_DETECTION.SBDA_SCDA_BOUNDARY_UTC:
    % Among records with BDN=4, check sliding windows for the min-max
    % difference for BIAS HK's measured bias current (within the entire
    % window). For windows for which the min-max difference exceeds a
    % threshold, label the entire window (where BDM=4) as sweeping.
    %
    % NOTE: There is some unimportant "imprecision" in the window algorithm
    % which can be seen as an unimportant bug. Records which are BDM<>4 (i.e.
    % definitely not sweep) within a window with currents exceeding the
    % thresholds (for BDM=4) are labelled as sweeps.
    %
    function isSweepingScda = SCDA_wo_margins(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso)
      scdaBeginTt2000          = spdfcomputett2000(Bso.get_fv('PROCESSING.L2.SWEEP_DETECTION.SBDA_SCDA_BOUNDARY_UTC'));
      windowLengthHkCdfRecords =                   Bso.get_fv('PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_LENGTH_HK_CDF_RECORDS');
      % Minimum min-max difference for counting as sweep.
      currentMmDiffMinimumTm   =                   Bso.get_fv('PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM');

      assert(round(windowLengthHkCdfRecords) == windowLengthHkCdfRecords)
      % NOTE: Must be at least two since SCDA requires comparison between
      %       timestamps.
      assert(windowLengthHkCdfRecords        >= 2)
      assert(currentMmDiffMinimumTm >= 0)



      nCdfRecs = irf.assert.sizes(...
        hkTt2000,         [-1, 1], ...
        hkBdmFpa,         [-1, 1], ...
        hkBiasCurrentFpa, [-1, 3]);



      hkBiasCurrent = hkBiasCurrentFpa.int2doubleNan();
      bdm           = hkBdmFpa.int2doubleNan();
      % Whether SCDA applies to (should be used for) records.
      bScdaApplies  = hkTt2000 > scdaBeginTt2000;



      %========================================================================
      % Use a moving window, label those windows for which a condition is true
      %========================================================================
      isSweepingScda = false(size(hkTt2000));    % Preallocate
      % NOTE: Will iterate zero times if window is longer than number of
      %       records.
      for i1 = 1:(nCdfRecs-(windowLengthHkCdfRecords-1))
        i2        = i1 + (windowLengthHkCdfRecords-1);
        iWindowAr = i1:i2;

        % NOTE: Treating all data in window combined, both in time AND over
        % channels/antennas.
        hkBiasCurrentWindow = hkBiasCurrent(i1:i2, :);

        % Exclude bias currents when sweep is impossible.
        % -----------------------------------------------
        % NOTE: This does *not* have the same effect as removing labels based
        % on BDM as done later. This affects windows which cover both BDM=4 and
        % BDM<>4.
        hkBiasCurrentWindow(bdm(iWindowAr) ~= bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE, :) = NaN;

        % NOTE: Calculating min/max and diffs on each bias/antenna channel
        % separately, i.e. detecting whether channels separately vary a lot.
        minWindowAr = min(hkBiasCurrentWindow, [], 1);
        maxWindowAr = max(hkBiasCurrentWindow, [], 1);
        mmDiffAr    = maxWindowAr - minWindowAr;

        % NOTE: Can not reduce only labelling.
        if any(mmDiffAr >= currentMmDiffMinimumTm)
          isSweepingScda(iWindowAr) = true;
        end
      end

      % Remove sweep labelling for selected parts.
      isSweepingScda = isSweepingScda & (bdm == bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE) & bScdaApplies;
    end



    % Try to autodetect sweeps using SBDA (BDM) and SCDA (BIAS HK bias
    % currents).
    %
    %
    % ARGUMENTS
    % =========
    % Bso
    %       NOTE: PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_LENGTH_HK_CDF_RECORDS:
    %       If greater than the number of CDF records/rows of data, then no
    %       record will be labelled as sweeping.
    %
    function isSweepingFpa = SBDA_SCDA_with_margins(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso)
      % TODO-DEC: Does having argument and return value FPAs make sense?
      %           Should caller convert?
      %
      % PROPOSAL: Separate window margins for before and after window.
      %   PRO: Margin after needs to be longer.
      %   NOTE: Looking at sweep for 2024-06-21, margins should be maybe:
      %       before sweep proper:  ~1-2 min
      %       after sweep proper:   ~6-7 min
      % PROPOSAL: Sweep detection algorithm which uses (and labels) the data gaps
      %           before & after the sweep.
      % PROPOSAL: Length of margins should be set in time, not HK CDF records.

      windowMarginSec = Bso.get_fv('PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_MARGIN_SEC');

      % Detect sweeps using SBDA.
      isSweepingSbda = bicas.proc.L1L2.swpdet.SBDA_wo_margins(hkTt2000, hkBdmFpa, Bso);

      % Detect sweeps using SCDA.
      isSweepingScda = bicas.proc.L1L2.swpdet.SCDA_wo_margins(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso);

      % Merge results and add margins.
      isSweeping           = isSweepingSbda | isSweepingScda;
      isSweepingWithMargin = irf.utils.true_with_margin( ...
        hkTt2000, isSweeping, windowMarginSec * 1e9, windowMarginSec * 1e9);

      isSweepingFpa        = bicas.utils.FPArray(isSweepingWithMargin);
    end



  end    % methods(Static)



end
