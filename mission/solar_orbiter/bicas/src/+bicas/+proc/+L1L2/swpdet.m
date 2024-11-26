%
% Collection of functions related to detecting and labelling bias sweeps.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef swpdet
  % PROPOSAL: More automatic test code.
  %
  % PROPOSAL: Rename PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC to reference both
  %           SBDA and SCDA.
  %   PROPOSAL: SBDA_SCDA_BOUNDARY
  %
  % PROPOSAL: Separate function(s) for detecting sweeps, adding margin,
  %           converting to science time(?).
  %   * Detect sweeps using L1R QUALITY_BITMASK. Return value in science time.
  %   * Add time margins to time interval (t_before, t_after).



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



    % Detect
    function isSweepingSbda = SBDA_wo_margins(tt2000, bdmFpa, Bso)
      % PROPOSAL: Use term SBDA in function name.

      sbdaEndTt2000  = spdfcomputett2000(Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC'));

      irf.assert.sizes(...
        tt2000, [-1, 1], ...
        bdmFpa, [-1, 1]);

      bdm            = bdmFpa.int2doubleNan();
      bSbdaApplies   = tt2000 <= sbdaEndTt2000;

      isSweepingSbda = bSbdaApplies & (bdm == bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE);
    end



    function isSweepingScda = SCDA_wo_margins(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso)
      % Time before which a sweep is equivalent to BDM=4. Inclusive threshold.
      sbdaEndTt2000          = spdfcomputett2000(Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC'));
      windowLengthPts        =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS');
      % Minimum min-max difference for counting as sweep.
      currentMmDiffMinimumTm =                   Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM');

      assert(round(windowLengthPts) == windowLengthPts)
      % NOTE: Must be at least two since SCDA requires comparison between
      %       timestamps.
      assert(windowLengthPts        >= 2)
      assert(currentMmDiffMinimumTm >= 0)



      nCdfRecs = irf.assert.sizes(...
        hkTt2000,         [-1, 1], ...
        hkBdmFpa,         [-1, 1], ...
        hkBiasCurrentFpa, [-1, 3]);



      hkBiasCurrent = hkBiasCurrentFpa.int2doubleNan();
      bdm           = hkBdmFpa.int2doubleNan();
      % Whether SCDA applies to (should be used for) records.
      bScdaApplies  = hkTt2000 > sbdaEndTt2000;



      isSweepingScda = false(size(hkTt2000));    % Preallocate
      % NOTE: Will iterate zero times if window is longer than number of
      %       records.
      for i1 = 1:(nCdfRecs-(windowLengthPts-1))
        i2        = i1 + (windowLengthPts-1);
        iWindowAr = i1:i2;

        % NOTE: Treating all data in window combined, both in time AND over
        % channels/antennas.
        hkBiasCurrentWindow = hkBiasCurrent(i1:i2, :);

        hkBiasCurrentWindow(bdm(iWindowAr) ~= bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE, :) = NaN;
        minWindow = min(hkBiasCurrentWindow, [], 1);
        maxWindow = max(hkBiasCurrentWindow, [], 1);
        mmDiff    = maxWindow - minWindow;

        % NOTE: Can not reduce only labelling.
        if any(mmDiff >= currentMmDiffMinimumTm)
          % isSweepingScda(iWindowAr) = isSweepingScda(iWindowAr) | ((bdm(iWindowAr) == bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE) & bScdaApplies(iWindowAr));
          isSweepingScda(iWindowAr) = true;
        end
      end

      isSweepingScda = isSweepingScda & (bdm == bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE) & bScdaApplies;
    end



    % Try to autodetect sweeps using BDM and BIAS HK.
    %
    % NOTE: This function is intended to be temporary functionality while waiting
    % for the long-term solution which is to use the relevant bit in L1/L1R
    % QUALITY_BITMASK for detecting sweeps. Since it is a (hopefully) short-term
    % solution, the functionality is also not as sophisticated (and presumably
    % accurate) or configurable as it could be.
    %
    % NOTE: There is some unimportant "imprecision" in the window algorithm
    % which can be seen as an unimportant bug. Records which are BDM<>4 (i.e.
    % definately not sweep) within a window with currents exceeding the
    % thresholds (for BDM=4) are labelled as sweeps.
    %
    % NOTE: 2024-11-26: It *appear* as if L1 QUALITY_BITMASK in sample datasets
    % (not production) contains bits for sweeps but without time margins before and
    % after. Sweep bits in QUALITY_BITMASK has not been confirmed yet.
    %
    %
    % ALGORITHM
    % =========
    % PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC specifies a timestamp.
    % Records before timestamp:
    %   BDN=4 <=> sweep
    % Records after timestamp:
    %   Among records with BDN=4, check sliding windows for the min-max
    %   difference for BIAS HK's measured bias current (within the entire
    %   window). For windows for which the min-max difference exceeds a
    %   threshold, label entire window (where BDM=4) as sweeping.
    %
    %
    % ARGUMENTS
    % =========
    % Bso
    %       NOTE: PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS: If
    %       greater than the number of CDF records/rows of data, then no
    %       record will be labelled as sweeping.
    %
    function isSweepingFpa = SBDA_SCDA_with_margins(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso)
    % PROPOSAL: Use SBDA, SCDA in function name.
    %
    % TODO-DEC: Does having argument and return value FPAs make sense?
    %           Should caller convert?
    % PROPOSAL: Rename (and negate) bicas.proc.L1L2.swpdet.BDM_SWEEP_POSSIBLE --> BDM_SWEEP_IMPOSSIBLE
    %   CON: Unintuitive for time period when sweep can be deduced frmo BDM.
    % PROPOSAL: Rename "PTS" (unit) for PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS.
    %   PRO: Unclear
    %   PROPOSAL: "HK_CDF_RECORDS", HkCdfRecords
    %     CON: Long
    %
    % PROPOSAL: Separate window margins for before and after window.
    %   PRO: Margin after needs to be longer.
    %   NOTE: Looking at sweep for 2024-06-21, margins should be maybe:
    %       before sweep proper:  ~1-2 min
    %       after sweep proper:   ~6-7 min
    % PROPOSAL: Sweep detection algorithm which uses (and labels) the data gaps
    %           before & after the sweep.
    % PROPOSAL: Length of margins shouls be set in time, not HK CDF records.



    windowMarginSec = Bso.get_fv('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MARGIN_SEC');



    % Detect sweeps using SBDA.
    isSweepingSbda = bicas.proc.L1L2.swpdet.SBDA_wo_margins(hkTt2000, hkBdmFpa, Bso);

    % Detect sweeps using SCDA.
    isSweepingScda = bicas.proc.L1L2.swpdet.SCDA_wo_margins(hkTt2000, hkBdmFpa, hkBiasCurrentFpa, Bso);



    isSweeping           = isSweepingSbda | isSweepingScda;
    isSweepingWithMargin = irf.utils.true_with_margin( ...
      hkTt2000, isSweeping, windowMarginSec * 1e9, windowMarginSec * 1e9);

    isSweepingFpa        = bicas.utils.FPArray(isSweepingWithMargin);
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)
  end    % methods(Static, Access=private)



end
