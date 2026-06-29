%
% Miscellaneous internal "utility functions" used for implementing spin fitting
% (bepic.spinfit).
%
%
% IMPLEMENTATION NOTES
% ====================
% The code assumes that fit windows can overlap, and can therefore not assume
% that the end of one fit window is not the beginning of another.
% --
% The code does not (deeply) assume that output time stamps must be the middle
% of the fit window. The code therefore avoids using this as a variable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils
  % PROPOSAL: Handle absence of spin phase.
  %     Ex: In eclipse?
  %   PROPOSAL: Use default spin.
  %   PROPOSAL: Extrapolate spin.
  %     NOTE: Only works if spin phase data ceases (or returns)
  %             within the time interval currently be processed.
  %     CON/PROBLEM: Spin rate may change in eclipse.
  %   PROPOSAL: Remove samples without spin phase.
  % PROPOSAL: Generera fit windows mha algorithm which expliticly steps over
  %           samples/time. Avoid working with arrays and CSP.
  %   PRO: Could be more flexible w.r.t. to handling not yet known edge cases.
  %     Ex: Spin phase. Unknown.
  %     Ex: Ignore spin phase if behaving strangely?
  %       PROPOSAL: Detect and then treat as if spin hase unknown.
  %     CON/PROBLEM: Can not extrapolate CSP to timestamps without CSP arrays.
  %     CON: Harder to determine overlapping fit windows?
  %     CON: Still makes sense to split processing into simple intervals where
  %          CSP is well-defined.
  %     NOTE: Any such hacks will (likely) be different from any counterparts
  %           for WPT. ==> Different handling of edge cases, timestamps, data
  %           gaps.




  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Convert fit windows specified as pairs of timestamps (boundaries) to
    % ranges of sample indices.
    %
    % DESIGN NOTE
    % ===========
    % This exists as a separate reusable function since it could be applied to
    % any algorithm for identifying fit windows.
    %
    function [iBeginAr, iEndAr] = fit_window_time_to_indices(...
        tt2000Ar, beginTt2000Ar, endTt2000Ar)

      assert(iscolumn(tt2000Ar)      & isa(tt2000Ar,      "int64"))
      assert(iscolumn(beginTt2000Ar) & isa(beginTt2000Ar, "int64"))
      assert(iscolumn(endTt2000Ar)   & isa(endTt2000Ar,   "int64"))

      assert(issorted(tt2000Ar, "STRICTASCEND"))    % Required by algorithm.
      nSamples = numel(tt2000Ar);

      % ----------------------------------------------------------------
      % Given timestamps of beginning and end of fit windows, derive the
      % ranges of samples indices for each fit window.
      % ----------------------------------------------------------------
      % IMPLEMENTATION NOTE: Could have been implemented by iterating over fit
      % windows and comparing time variables with each other, e.g.
      %     b = (tt2000Begin <= tt2000Ar) & (tt2000Ar <= tt2000End))
      % . This might however scale badly with the size of the data. Using
      % interp1() instead should hopefully avoid bad performance.
      iBeginEnd = interp1(...
        double(tt2000Ar), 1:nSamples, ...
        double([beginTt2000Ar, endTt2000Ar]), ...
        "linear", "extrap");

      iBeginAr = max(ceil( iBeginEnd(:, 1)), 1);
      iEndAr   = min(floor(iBeginEnd(:, 2)), nSamples);
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
    % Table of fit window boundaries (TT2000; not CSP (ill-defined) or indices).
    % Fit windows can only begin at
    % fitWindowBeginRefRad + m * fitWindowPeriodRad, where m=integer.
    %
    function FitWindowTable = get_SAFWs(A)
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowBeginRefRad
        A.dataGapMinNs
      end

      % =========
      % ALGORITHM
      % =========
      [iBeginAr, iEndAr, nSegments] = bepic.spinfit.split_by_data_gap(...
        A.tt2000Ar, A.dataGapMinNs);

      % Create empty return variable.
      FitWindowTable = table();
      FitWindowTable.beginTt2000 = int64.empty( 0, 1);
      FitWindowTable.endTt2000   = int64.empty( 0, 1);

      for iSegment = 1:nSegments
        i = iBeginAr(iSegment) : iEndAr(iSegment);

        SegmentFitWindowTable = bepic.spinfit.utils.get_segment_SAFWs( ...
          tt2000Ar           = A.tt2000Ar(      i), ...
          spinPhaseRadAr     = A.spinPhaseRadAr(i), ...
          fitWindowPeriodRad = A.fitWindowPeriodRad, ...
          fitWindowLengthRad = A.fitWindowLengthRad, ...
          fitWindowBeginRad  = A.fitWindowBeginRefRad);

        FitWindowTable = [FitWindowTable; SegmentFitWindowTable];
      end

      % -----------------------------
      % Remove duplicated fit windows
      % -----------------------------
      % IMPLEMENTATION NOTE: Data gaps which
      % (1) are recognized by A.dataGapMinNs, and
      % (2) which are smaller than one revolution, and
      % (3) fit inside the time boundaries of one fit window,
      % lead to bepic.spinfit.utils.get_segment_SAFWs() identifying the same fit
      % window twice. They therefore need to be removed.
      % --
      % PROBLEM: Rounding errors could lead to de facto duplicated fit windows
      % still being kept?
      FitWindowTable = unique(FitWindowTable, "STABLE");
    end



    % Find beginning and end of SAFWs on time interval *without data gaps*
    % (jumps in time equal or greater than one revolution). (This is what
    % "segment" in the function name refers to.).
    %
    % IMPORTANT NOTE: Assumes that there are no data gaps.
    %
    %
    % RETURN VALUE
    % ============
    % Table
    %       NOTE: Zero fit windows if there are <= 1 samples, since the
    %       algorithm can then not extrapolate spin phase (CMP) values and can
    %       hence not find the beginning and end of the window even in
    %       principle.
    %
    function [FitWindowTable] = get_segment_SAFWs(A)
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowBeginRad
      end
      % PROPOSAL: Better name.
      %   ~(no) data gap
      %   get_SAFWs_no_data_gap

      nSamples = numel(A.tt2000Ar);

      assert(issorted(A.tt2000Ar))
      assert(all(isfinite(A.spinPhaseRadAr)))
      assert(numel(A.spinPhaseRadAr) == nSamples)

      % IMPORTANT NOTE: Can not calculate CSP if there are data gaps.
      cspRadAr = bepic.spinfit.spin_phase_to_cumulative_spin_phase(...
        A.spinPhaseRadAr);

      if nSamples <= 1
        % IMPLEMENTATION NOTE: Can not extrapolate CSP-->TT2000 even in
        % principle when there is zero or one sample. Therefore not returning
        % any fit window at all.

        fitWindowBeginTt2000Ar = int64.empty( 0, 1);
        fitWindowEndTt2000Ar   = int64.empty( 0, 1);
      else
        fitWindowBeginCspRadAr = bepic.spinfit.utils.get_incrementing_array(...
          cspRadAr(1), ...
          cspRadAr(end), ...
          A.fitWindowPeriodRad, ...
          A.fitWindowBeginRad);

        fitWindowEndCspRadAr   = fitWindowBeginCspRadAr + A.fitWindowLengthRad;

        fitWindowBeginTt2000Ar = int64(interp1(cspRadAr, double(A.tt2000Ar), fitWindowBeginCspRadAr, "LINEAR", "extrap"));
        fitWindowEndTt2000Ar   = int64(interp1(cspRadAr, double(A.tt2000Ar), fitWindowEndCspRadAr,   "LINEAR", "extrap"));
      end

      % ----------------------
      % Construct return value
      % ----------------------
      % IMPLEMENTATION NOTE: Not returning CSP since it is only well-defined for
      % the data which is used within this function (i.e. locally), and not over
      % large jumps in time (i.e. globally). The CSP arrays in this function can
      % therefore not be meaningfully combined to longer CSP arrays.
      FitWindowTable = table();
      FitWindowTable.beginTt2000 = fitWindowBeginTt2000Ar;
      FitWindowTable.endTt2000   = fitWindowEndTt2000Ar;
    end



    % Generate array of incrementing values x = xRef + n * xPeriod (n=integer),
    % for an interval, such that the lowest value is the highest legal
    % x <= xBegin.
    %
    %
    % IMPLEMENTATION NOTE
    % ===================
    % Designed to work for both integers and float for maximum flexibility, so
    % that it can support both e.g. TT2000 and cumulative sphin phase. This puts
    % constraints on the implementation though.
    %
    %
    % ARGUMENTS
    % =========
    % All arguments must be scalars, have the same MATLAB class, and may only be
    % floats or integers.
    %
    %
    % RETURN VALUE
    % ============
    % xArray
    %       Column array of all values x = xRef + n * xPeriod (n=integer) such
    %       that
    %       (1) the lowest value is the highest possible value which satisfies
    %       x <= xBegin, and
    %       (2) the highest value is the highest possible value which satisfies
    %       x <= xEnd.
    %
    function xAr = get_incrementing_array(xBegin, xEnd, xPeriod, xRef)
      % ==========
      % ASSERTIONS
      % ==========
      mc = class(xPeriod);
      assert(isa(xRef,   mc))
      assert(isa(xBegin, mc))
      assert(isa(xEnd,   mc))

      assert(isscalar(xPeriod) & (xPeriod > 0))
      assert(isscalar(xRef   ))
      assert(isscalar(xBegin ))
      assert(isscalar(xEnd   ) & (xBegin <= xEnd))

      % =========
      % ALGORITHM
      % =========
      % Derive the highest x such that
      % (1) x = xRef + m*xPeriod, and
      % (2) x<=xBegin.
      % NOTE: mod() works for both floats and integers.
      xFirst = xBegin - mod(xBegin-xRef, xPeriod);

      % NOTE: x:y:z can not produce a value outside [x,z].
      xAr = (xFirst : xPeriod : xEnd)';
    end



  end    % methods(Static)



end
