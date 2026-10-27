%
% Miscellaneous internal functions for generating fit windows ("fw").
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
classdef fw
  % PROPOSAL: Generate fit windows using algorithm which expliticly steps over
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
    function [iArCa] = fit_window_time_to_indices(...
        tt2000Ar, beginTt2000Ar, endTt2000Ar)
      % PROPOSAL: Change terminology "fit windows" --> "bins".
      %   CON: "Fit windows" can overlap.

      assert(iscolumn(tt2000Ar)      & isa(tt2000Ar,      "int64"))
      assert(iscolumn(beginTt2000Ar) & isa(beginTt2000Ar, "int64"))
      assert(iscolumn(endTt2000Ar)   & isa(endTt2000Ar,   "int64"))

      assert(issorted(tt2000Ar,      "STRICTASCEND"))    % Required by algorithm.
      % assert(issorted(beginTt2000Ar, "STRICTASCEND"))
      % assert(issorted(endTt2000Ar,   "STRICTASCEND"))
      % assert(all(beginTt2000Ar <= endTt2000Ar))
      assert(numel(beginTt2000Ar) == numel(endTt2000Ar))
      nSamples    = numel(tt2000Ar);
      nFitWindows = numel(beginTt2000Ar);

      % -----------------------------------------------------------------------
      % Given timestamps of beginning and end of fit windows, derive the ranges
      % of samples indices for each fit window.
      % -----------------------------------------------------------------------
      % IMPLEMENTATION NOTE: THIS CODE IS PERFORMANCE SENSITIVE FOR LARGE DATA!
      iArCa = cell(nFitWindows, 1);
      % IMPLEMENTATION NOTE: Add extra "data points" to makes sure that there
      % are always at least two data points, which is required by interp1().
      % This eliminates the need for explicitly handling special cases.
      xAr      = double([int64(-inf); tt2000Ar; int64(inf)]);
      yAr      = 0 : (nSamples+1);
      iBeginAr = interp1(xAr, yAr, double(beginTt2000Ar), "next");
      iEndAr   = interp1(xAr, yAr, double(endTt2000Ar  ), "next") - 1;

      % Convert the result to the expected return value format.
      for iFitWindow = 1:nFitWindows
        iArCa{iFitWindow} = (iBeginAr(iFitWindow) : iEndAr(iFitWindow))';
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
    % FitWindowTable
    %       Table of (1) fit window boundaries (TT2000; not CSP (ill-defined) or
    %       indices) and (2) fit window centers (TT2000).
    %       Fit windows can only begin at
    %       fitWindowBeginRefRad + m * fitWindowPeriodRad, where m=integer.
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
      assert(iscolumn(A.tt2000Ar)              & isa(A.tt2000Ar,              "int64"))
      assert(issorted(A.tt2000Ar, "STRICTASCEND"))
      assert(iscolumn(    A.spinPhaseRadAr)    & isa(A.spinPhaseRadAr,        "double"))
      assert(all(isfinite(A.spinPhaseRadAr)))
      % IMPLEMENTATION NOTE: Requiring spin phase on interval 0 to 2*pi. This
      % (1) is (analogous to) the spin phase values in L1pre CDF files which are
      % also bounded (0 to 360 deg), and (2) is an indirect check on the units
      % used, and (3) is required by the algorithm.
      assert(all((0 <= A.spinPhaseRadAr) & (A.spinPhaseRadAr <= 2*pi)))

      assert(isscalar(A.fitWindowPeriodRad)    & isa(A.fitWindowPeriodRad,    "double"))
      assert(         A.fitWindowPeriodRad > 0)
      assert(isscalar(A.fitWindowLengthRad)    & isa(A.fitWindowLengthRad,    "double"))
      assert(         A.fitWindowLengthRad > 0)
      assert(isscalar(A.fitWindowBeginRefRad)  & isa(A.fitWindowBeginRefRad,  "double"))
      assert(         A.fitWindowBeginRefRad >= 0)

      assert(isscalar(A.dataGapMinNs)          & isa(A.dataGapMinNs,          "int64"))
      assert(         A.dataGapMinNs > 0)

      % =========
      % ALGORITHM
      % =========
      [iBeginAr, iEndAr, nSegments] = bepic.spinfit.utils.split_by_data_gap(...
        A.tt2000Ar, A.dataGapMinNs);

      % Create empty return variable.
      FitWindowTable = table();
      FitWindowTable.beginTt2000  = int64.empty(0, 1);
      FitWindowTable.endTt2000    = int64.empty(0, 1);
      FitWindowTable.centerTt2000 = int64.empty(0, 1);

      % Iterate over data gap-separated samples.
      for iSegment = 1:nSegments
        i = iBeginAr(iSegment) : iEndAr(iSegment);

        SegmentFitWindowTable = bepic.spinfit.fw.get_SAFWs_when_no_data_gaps( ...
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
      % lead to bepic.spinfit.fw.get_SAFWs_when_no_data_gaps() identifying the same fit
      % window twice. They therefore need to be removed.
      % --
      % PROBLEM: Rounding errors could lead to de facto duplicated fit windows
      % still being kept?
      FitWindowTable = unique(FitWindowTable, "STABLE");
    end



    % Find beginning and end of SAFWs on time interval *without data gaps*
    % (jumps in time, equal or greater than one revolution).
    %
    % IMPORTANT NOTE: Assumes that there are no data gaps.
    %
    %
    % RETURN VALUE
    % ============
    % FitWindowTable
    %       NOTE: There will be zero fit windows if there are <= 1 samples,
    %       since the algorithm can then not extrapolate spin phase (CMP) values
    %       and can hence not find the beginning and end of the fit window even
    %       in principle.
    %
    function FitWindowTable = get_SAFWs_when_no_data_gaps(A)
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowBeginRad
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

        fitWindowBeginTt2000Ar  = int64.empty(0, 1);
        fitWindowEndTt2000Ar    = int64.empty(0, 1);
        fitWindowCenterTt2000Ar = int64.empty(0, 1);
      else
        fitWindowBeginCspRadAr = bepic.spinfit.utils.get_incrementing_array(...
          cspRadAr(1), ...
          cspRadAr(end), ...
          A.fitWindowPeriodRad, ...
          A.fitWindowBeginRad);

        fitWindowEndCspRadAr    = fitWindowBeginCspRadAr + A.fitWindowLengthRad;
        fitWindowCenterCspRadAr = (fitWindowBeginCspRadAr + fitWindowEndCspRadAr) / 2;

        fitWindowBeginTt2000Ar  = int64(interp1(cspRadAr, double(A.tt2000Ar), fitWindowBeginCspRadAr,  "LINEAR", "extrap"));
        fitWindowEndTt2000Ar    = int64(interp1(cspRadAr, double(A.tt2000Ar), fitWindowEndCspRadAr,    "LINEAR", "extrap"));
        fitWindowCenterTt2000Ar = int64(interp1(cspRadAr, double(A.tt2000Ar), fitWindowCenterCspRadAr, "LINEAR", "extrap"));
      end

      % ----------------------
      % Construct return value
      % ----------------------
      % IMPLEMENTATION NOTE: Not returning CSP since it is only well-defined for
      % the data which is used within this function (i.e. locally), and not over
      % large jumps in time (i.e. globally). The CSP arrays in this function can
      % therefore not be meaningfully combined to longer CSP arrays.
      FitWindowTable = table();
      FitWindowTable.beginTt2000  = fitWindowBeginTt2000Ar;
      FitWindowTable.endTt2000    = fitWindowEndTt2000Ar;
      FitWindowTable.centerTt2000 = fitWindowCenterTt2000Ar;
    end



  end    % methods(Static)



end
