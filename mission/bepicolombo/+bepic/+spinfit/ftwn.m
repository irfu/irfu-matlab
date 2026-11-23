%
% Miscellaneous internal functions for generating fit windows.
%
%
% IMPLEMENTATION NOTES
% ====================
% The code assumes that fit windows can overlap, and can therefore not assume
% that the end of one fit window is not the beginning of another.
% --
% The code does not (deeply) assume that output time stamps must be the middle
% of the fit window. The code therefore avoids (mostly) using this as a
% variable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ftwn
  % PROPOSAL: Generate fit windows using algorithm which expliticly steps over
  %           samples/time. Avoid working with arrays and CSP.
  %   PRO: Could be more flexible w.r.t. to handling not yet known edge cases.
  %     Ex: Spin phase. Unknown.
  %     Ex: Ignore spin phase if behaving strangely?
  %       PROPOSAL: Detect and then treat as if spin phase unknown.
  %     CON/PROBLEM: Can not extrapolate CSP to timestamps without CSP arrays.
  %     CON: Harder to determine overlapping fit windows?
  %     CON: Still makes sense to split processing into simple intervals where
  %          CSP is well-defined.
  %     NOTE: Any such hacks will (likely) be different from any counterparts
  %           for WPT. ==> Different handling of edge cases, timestamps, data
  %           gaps.
  %
  % PROPOSAL: Generate SAFWs by iterating over time. Not vectorized.
  % PROPOSAL: Generate SAFWs in a way which is more analogous to the
  %           Japanese way/find_spin_zero_time.
  %   PRO: Easier to get agreement with Japan code, initially, and with
  %        future updates/changes.
  %     Ex: Timestamps, removal of FTWNs based on spin phase, number of
  %         samples (in FTWN), (input) timestamps(?!), but not sample values
  %         themselves.
  %   PROPOSAL: Only identify being+end indices (not timestamps). Identity
  %             center timestamp with fit.
  %     NOTE: Currently, bepic.spinfit.fit_array() uses
  %           bepic.spinfit.ftwn.FTWN_time_to_indices() to convert
  %           FTWN timestamps-->indices which all that is truly used anyway.
  %     PRO: Can abolish bepic.spinfit.FTWN_time_to_indices().
  %     PRO: Can still generalize to overlapping FTWNs.
  %     CON: Can not use fit to return (true) non-center timestamp instead.
  %       NOTE: Arguably, using a fit to return a center timestamp does not
  %             truly work (if one is rigorous).
  %       Ex: FTWN boundary timestamps can be inconsistent. The end timestamps
  %           of one FTWN is not automatically the begin timestamp of
  %           the FTWN after.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Convert fit windows, specified as pairs of timestamps (begin+end), to
    % ranges of sample indices.
    %
    %
    % DESIGN NOTE
    % ===========
    % This exists as a separate reusable function since it could be applied to
    % any algorithm for identifying fit windows.
    %
    %
    % RETURN VALUES
    % =============
    % iArCa
    %       Column cell array of column arrays of indices.
    %       iCa{i} = Indices into tt2000Ar, for fit window i.
    %
    function [iArCa] = FTWN_time_to_indices(...
        tt2000Ar, beginTt2000Ar, endTt2000Ar)
      % IMPLEMENTATION NOTE: Not using terminology of "bins" since FTWNs can
      % overlap.

      assert(iscolumn(tt2000Ar)      & isa(tt2000Ar,      "int64"))
      assert(iscolumn(beginTt2000Ar) & isa(beginTt2000Ar, "int64"))
      assert(iscolumn(endTt2000Ar)   & isa(endTt2000Ar,   "int64"))

      assert(issorted(tt2000Ar,      "STRICTASCEND"))    % Required by algorithm.
      % assert(issorted(beginTt2000Ar, "STRICTASCEND"))
      % assert(issorted(endTt2000Ar,   "STRICTASCEND"))
      % assert(all(beginTt2000Ar <= endTt2000Ar))
      assert(numel(beginTt2000Ar) == numel(endTt2000Ar))
      nSamples = numel(tt2000Ar);
      nFtwn    = numel(beginTt2000Ar);

      % -----------------------------------------------------------------------
      % Given timestamps of beginning and end of fit windows, derive the ranges
      % of samples indices for each fit window.
      % -----------------------------------------------------------------------
      % IMPLEMENTATION NOTE: THIS CODE IS PERFORMANCE SENSITIVE FOR LARGE DATA!
      iArCa = cell(nFtwn, 1);
      % IMPLEMENTATION NOTE: Add extra "data points" to makes sure that there
      % are always at least two data points, which is required by interp1().
      % This eliminates the need for explicitly handling special cases.
      xAr      = double([int64(-inf); tt2000Ar; int64(inf)]);
      yAr      = 0 : (nSamples+1);
      iBeginAr = interp1(xAr, yAr, double(beginTt2000Ar), "next");
      iEndAr   = interp1(xAr, yAr, double(endTt2000Ar  ), "next") - 1;

      % Convert the result to the expected return value format.
      for iFtwn = 1:nFtwn
        iArCa{iFtwn} = (iBeginAr(iFtwn) : iEndAr(iFtwn))';
      end
    end



    % Find SAFWs. Data is allowed to contain data gaps (jumps in time).
    %
    %
    % DESIGN NOTE
    % ===========
    % Function is meant to be an example of a function which generates fit
    % windows according to a certain algorithm. Other functions could be made
    % for producing fit windows using other algorithms (e.g. TAFWs). The
    % interface is chosen such that it should be easy to replace this function
    % by any other such function..
    %
    %
    % RETURN VALUES
    % =============
    % FtwnTable
    %       Table of (1) fit window boundaries (TT2000; not CSP (ill-defined) or
    %       indices) and (2) fit window centers (TT2000).
    %       Fit windows can only begin at
    %       ftwnBeginRefRad + m * ftwnPeriodRad, where m=integer.
    %
    function FtwnTable = get_SAFWs(A)
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.ftwnPeriodRad
        A.ftwnLengthRad
        A.ftwnBeginRefRad
        A.dataGapMinNs
      end
      assert(iscolumn(A.tt2000Ar)             & isa(A.tt2000Ar,       "int64"))
      assert(issorted(A.tt2000Ar, "STRICTASCEND"))
      assert(iscolumn(    A.spinPhaseRadAr)   & isa(A.spinPhaseRadAr, "double"))
      assert(all(isfinite(A.spinPhaseRadAr)))
      % IMPLEMENTATION NOTE: Requiring spin phase on interval 0 to 2*pi. This
      % (1) is (analogous to) the spin phase values in L1p CDF files which are
      % also bounded (0 to 360 deg), and (2) is an indirect check on the units
      % used, and (3) is required by the algorithm.
      assert(all((0 <= A.spinPhaseRadAr) & (A.spinPhaseRadAr <= 2*pi)))

      assert(isscalar(A.ftwnPeriodRad)      & isa(A.ftwnPeriodRad,    "double"))
      assert(         A.ftwnPeriodRad > 0)
      assert(isscalar(A.ftwnLengthRad)      & isa(A.ftwnLengthRad,    "double"))
      assert(         A.ftwnLengthRad > 0)
      assert(isscalar(A.ftwnBeginRefRad)    & isa(A.ftwnBeginRefRad,  "double"))
      assert(         A.ftwnBeginRefRad >= 0)

      assert(isscalar(A.dataGapMinNs)       & isa(A.dataGapMinNs,     "int64"))
      assert(         A.dataGapMinNs > 0)

      % =========
      % ALGORITHM
      % =========
      [iBeginAr, iEndAr, nSegments] = bepic.spinfit.utils.split_by_data_gap(...
        A.tt2000Ar, A.dataGapMinNs);

      % Create empty return variable.
      FtwnTable = table();
      FtwnTable.beginTt2000  = int64.empty(0, 1);
      FtwnTable.endTt2000    = int64.empty(0, 1);
      FtwnTable.centerTt2000 = int64.empty(0, 1);

      % Iterate over data gap-separated samples.
      for iSegment = 1:nSegments
        i = iBeginAr(iSegment) : iEndAr(iSegment);

        SegmentFtwnTable = bepic.spinfit.ftwn.get_SAFWs_when_no_data_gaps( ...
          tt2000Ar       = A.tt2000Ar(      i), ...
          spinPhaseRadAr = A.spinPhaseRadAr(i), ...
          ftwnPeriodRad  = A.ftwnPeriodRad, ...
          ftwnLengthRad  = A.ftwnLengthRad, ...
          ftwnBeginRad   = A.ftwnBeginRefRad);

        FtwnTable = [FtwnTable; SegmentFtwnTable];
      end

      % -----------------------------
      % Remove duplicated fit windows
      % -----------------------------
      % IMPLEMENTATION NOTE: Data gaps which
      % (1) are recognized by A.dataGapMinNs, and
      % (2) which are smaller than one revolution, and
      % (3) fit inside the time boundaries of one fit window,
      % lead to bepic.spinfit.ftwn.get_SAFWs_when_no_data_gaps() identifying
      % the same fit window twice. They therefore need to be removed.
      % --
      % PROBLEM: Rounding errors could lead to de facto duplicated fit windows
      %          still being kept?
      FtwnTable = unique(FtwnTable, "STABLE");
    end



    % Find beginning and end of SAFWs on time interval *without data gaps*
    % (jumps in time, equal or greater than one revolution).
    %
    % IMPORTANT NOTE: Assumes that there are no data gaps.
    %
    %
    % RETURN VALUE
    % ============
    % FtwnTable
    %       NOTE: There will be zero fit windows if there are <= 1 samples,
    %       since the algorithm can then not extrapolate spin phase (CMP)
    %       values and can hence not find the beginning and end of the fit
    %       window even in principle.
    %
    function FtwnTable = get_SAFWs_when_no_data_gaps(A)
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.ftwnPeriodRad
        A.ftwnLengthRad
        A.ftwnBeginRad
      end

      nSamples = numel(A.tt2000Ar);

      assert(issorted(A.tt2000Ar))
      assert(all(isfinite(A.spinPhaseRadAr)))    % NOTE: No NaN.
      assert(numel(A.spinPhaseRadAr) == nSamples)

      % IMPORTANT NOTE: Can not calculate CSP if there are data gaps.
      cspRadAr = bepic.spinfit.utils.spin_phase_to_CMP(...
        A.spinPhaseRadAr);

      if nSamples <= 1
        % IMPLEMENTATION NOTE: One can not extrapolate CSP-->TT2000 even in
        % principle when there is zero or one sample. Therefore not returning
        % any fit window at all.

        ftwnBeginTt2000Ar  = int64.empty(0, 1);
        ftwnEndTt2000Ar    = int64.empty(0, 1);
        ftwnCenterTt2000Ar = int64.empty(0, 1);
      else
        % CASE: There are at least TWO samples.
        ftwnBeginCspRadAr = bepic.spinfit.utils.get_incrementing_array(...
          cspRadAr(1), ...
          cspRadAr(end), ...
          A.ftwnPeriodRad, ...
          A.ftwnBeginRad);

        ftwnEndCspRadAr    = ftwnBeginCspRadAr + A.ftwnLengthRad;
        ftwnCenterCspRadAr = (ftwnBeginCspRadAr + ftwnEndCspRadAr) / 2;

        ftwnBeginTt2000Ar  = int64(interp1(cspRadAr, double(A.tt2000Ar), ftwnBeginCspRadAr,  "LINEAR", "extrap"));
        ftwnEndTt2000Ar    = int64(interp1(cspRadAr, double(A.tt2000Ar), ftwnEndCspRadAr,    "LINEAR", "extrap"));
        ftwnCenterTt2000Ar = int64(interp1(cspRadAr, double(A.tt2000Ar), ftwnCenterCspRadAr, "LINEAR", "extrap"));
      end



      % ----------------------
      % Construct return value
      % ----------------------
      % IMPLEMENTATION NOTE: Not returning CSP since it is only well-defined for
      % the data which is used within this function (i.e. locally), and not over
      % large jumps in time (i.e. globally). The CSP arrays in this function can
      % therefore not be meaningfully combined to longer CSP arrays.
      FtwnTable = table();
      FtwnTable.beginTt2000  = ftwnBeginTt2000Ar;
      FtwnTable.endTt2000    = ftwnEndTt2000Ar;
      FtwnTable.centerTt2000 = ftwnCenterTt2000Ar;
    end



    % Get indices to the begin and end of FTWNs, assuming that
    % (1) spin phase is always 0 to 2pi radians,
    % (2) FTWNs are defined as spanning from 0 to 2pi radians (and being 2pi
    %     long).
    % (3) FTWN center timestamps are obtained using a fit.
    %
    % The algorithm is deliberately similar to the Japanese WPT function
    % find_spin_zero_time() (see Shinbori's e-mail 2026-11-02) so
    %
    %
    % DESIGN NOTES
    % ============
    % The design is deliberately similar to Atsuki Shinbori's
    % find_spin_zero_time.m, which should be a translation from IDL(?) to
    % MATLAB of his WPT code. See Shinbori's e-mail 2026-11-02. The advantages
    % are that it can then return similar results, and be updated to return
    % even more similar results if necessary. For this reason, it can not be
    % configured to use arbitrary FTWNs (w.r.t. length, location, overlap) or
    % generate arbitrary FTWN output timestamps (has to be center).
    % It is however not identical, and does not have all the same features.
    %
    % IMPLEMENTATION NOTE: This is not implemented as a separate argument
    % (~nMinSamplesPerFtwn) since such a feature is better implemented
    % elsewhere (when doing fits; also better for testing).
    %
    % IMPLEMENTATION NOTE: The code deliberately does not implement features
    % which are better implemented in other locations, e.g. min number of
    % samples per FTWN (which also is unaware of NaN samples).
    %
    % NOTE: The code can (probably) be generalized to support arbitrary
    % position+length FTWNs, overlapping FTWNs, and arbitrary position output
    % timestamps by using CMPs.
    function FtwnTable = get_SAFWs_limited(A)
      %{
      PROPOSAL: Better name wrt. specifying that it is not configurable wrt.
                FTWN length+positions, output timestamp position.

      PROPOSAL: Exclude FTWNs adjacent to data gaps (similar to Shinbori's
                code).
        PRO: Can not be done outside of this function.

      PROPOSAL: Implement using CMP.
        NOTE: The (current) CMP derivation assumes that all decrements imply
        adding one rev. to CMP. Therefore, the (potential) implementation of
        permitting small decrements affects also the CMP derivation code.
      %}
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.dataGapMinNs
      end

      assert(iscolumn(A.tt2000Ar)        & isa(     A.tt2000Ar,       "int64"))
      assert(issorted(A.tt2000Ar, "STRICTASCEND"))
      assert(iscolumn(A.spinPhaseRadAr)  & isa(     A.spinPhaseRadAr, "double"))
      assert(all((0 <= A.spinPhaseRadAr) & (A.spinPhaseRadAr <= 2*pi)))
      assert(isscalar(A.dataGapMinNs)    & isa(A.dataGapMinNs,        "int64"))
      assert(         A.dataGapMinNs > 0)



      % Normalize 2*pi to 0
      % -------------------
      % IMPLEMENTATION NOTE: Assertion permits both spin phase = 0 and 2*pi,
      % which represent the same information, but influence the algorithm.
      % This normalization
      % (1) permits L1p data with spin phase=360 deg (which might happen; not
      %     sure)
      % (2) simplifies(?) automated tests, since wrapTo2Pi(2*pi)==2*pi (not
      %     zero). The alternative is to use mod(spinPhaseRad, 2*pi) in tests.
      % A.spinPhaseRadAr(A.spinPhaseRadAr == 2*pi) = 0;

      % Minimum number of required samples per FTWN
      % -------------------------------------------
      % NOTE: The linear fit required for obtaining the center timestamp of
      % the FTWN requires at least two samples.
      N_SAMPLES_PER_FTWN_MIN = 2;



      % Naming conventions
      % ------------------
      % b  = Logical index(es) into array of samples.
      % bb = Logical index(es) into array of FTWNs.
      % i  = Linear index(es) into array of samples
      % j  = Linear index(es) into array of FTWNs.

      nSamples = numel(A.tt2000Ar);

      if nSamples == 0
        nFtwn1 = 0;

        iFtwn1BeginAr = double.empty(0, 1);
        iFtwn1EndAr   = double.empty(0, 1);
      else
        % Sample indices before a data gap.
        bDataGapAfterAr       = [diff(A.tt2000Ar) > A.dataGapMinNs; false];
        % Sample indices before the spin phase wraps (decreases).
        bSpinPhaseWrapAfterAr = [diff(A.spinPhaseRadAr) < 0; false];

        % Indices before boundaries *between* preliminary FTWNs.
        % NOTE: (inner) "boundary" is defined as either
        % (1) spin phase wrapping, or
        % (2) data gap
        bInnerBoundaryAfterAr = bDataGapAfterAr | bSpinPhaseWrapAfterAr;
        iInnerBoundaryAfterAr = find(bInnerBoundaryAfterAr);

        % Define preliminary FTWNs: "FTWN1"
        % ---------------------------------
        % NOTE: Some of these may later be rejected and thus excluded from the
        %       return value. Hence the separate name.
        iFtwn1BeginAr = [1; iInnerBoundaryAfterAr+1];
        iFtwn1EndAr   = [iInnerBoundaryAfterAr; nSamples];
        nFtwn1        = numel(iFtwn1BeginAr);

        % Create arrays over FWTNs for whether they are adjacent to (1) data
        % gaps, or (2) spin phase wraps.
        bbDataGapAfterAr    = bDataGapAfterAr(iFtwn1EndAr);
        bbDataGapBeforeAr   = [false; bbDataGapAfterAr(1:end-1)];
        bbDataGapAdjacentAr = bbDataGapBeforeAr | bbDataGapAfterAr;

        % bbSpinPhaseWrapAfterAr    = bSpinPhaseWrapAfterAr(iFtwn1EndAr);
        % bbSpinPhaseWrapBeforeAr   = [false; bbSpinPhaseWrapAfterAr(1:end-1)];
        % bbSpinPhaseWrapAdjacentAr = bbSpinPhaseWrapBeforeAr | bbSpinPhaseWrapAfterAr;
      end

      % Preallocate
      ftwn1CenterTt2000Ar = int64(zeros(nFtwn1, 1));
      bbKeepAr            = true(nFtwn1, 1);

      % Disable warning.
      % Used in find_spin_zero_time.m (Japan).
      statePolyfitWarning = ...
        warning('query','MATLAB:polyfit:RepeatedPointsOrRescale').state;
      warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');

      % ==================
      % Iterate over FTWN1
      % ==================
      for jFtwn1 = 1:nFtwn1
        iFtwnAr       = iFtwn1BeginAr(jFtwn1) : iFtwn1EndAr(jFtwn1);
        nSamplesFtwn1 = length(iFtwnAr);

        % Eliminate FTWNs adjacent to data gaps.
        if bbDataGapAdjacentAr(jFtwn1)
          bbKeepAr(jFtwn1) = false;
          continue
        end

        % Eliminate FTWNs with too few samples (as required by algorithm; not
        % configurable).
        if nSamplesFtwn1 < N_SAMPLES_PER_FTWN_MIN
          bbKeepAr(jFtwn1) = false;
          continue
        end

        ftwnTt2000Ar       = A.tt2000Ar(      iFtwnAr);
        ftwnSpinPhaseRadAr = A.spinPhaseRadAr(iFtwnAr);

        % ==============================================
        % Derive the timestamp of the center of the FTWN
        % ==============================================
        % Linear fit on data points (timestamp, spin phase).
        % NOTE: May trigger warning "MATLAB:polyfit:RepeatedPointsOrRescale"
        %       though the data seems legit.
        [cAr, ~] = polyfit(double(ftwnTt2000Ar), ftwnSpinPhaseRadAr, 1);
        radPerNs  = cAr(1);
        offsetRad = cAr(2);

        % pi = offsetRad + radPerNs*centerTt2000
        % ==> centerTt2000 = (pi - offsetRad) / radPerNs
        ftwn1CenterTt2000Ar(jFtwn1) = int64((pi - offsetRad) / radPerNs);
      end
      warning(statePolyfitWarning, 'MATLAB:polyfit:RepeatedPointsOrRescale');

      FtwnTable = table();
      FtwnTable.iBegin       = iFtwn1BeginAr(      bbKeepAr, 1);
      FtwnTable.iEnd         = iFtwn1EndAr(        bbKeepAr, 1);
      FtwnTable.centerTt2000 = ftwn1CenterTt2000Ar(bbKeepAr, 1);
    end



  end    % methods(Static)



end
