%
% Functions which implement the core processing for MEFISTO spin fits, when
% processing L1p CDFs --> L2pre CDFs.
%
% The purpose of this file is to isolate many of the technical details
% surrounding the details of the spin fit, e.g. selecting the fit windows (time
% windows) used for spin fitting, and which exact timestamps should represent
% fit windows.
%
% The implementation uses mms_spinfit_m() for the core processing. In the long
% term, this function should probably be replaced by a copy dedicated to
% BepiColombo processing (in case it needs to be modified).
%
%
% TECHNICALITIES WHICH ARE BEYOND THE SCOPE OF THIS CLASS
% =======================================================
% * Rotating the coordinate system.
% * Converting first-order sin/cos fit coefficients to an E field vector.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit
  % PROPOSAL: Better class name
  %   ~spin fit
  %   MEFISTO, Mio
  %     CON: Mission is already implied by the parent package "bepic".
  %   ~spin fit "engine", "core"
  %   fit
  %   process
  %
  % PROPOSAL: Convert class to package of functions dedicated to official
  %           processing.
  %   CON: Can not put constants in class.
  %     CON: Can create dedicated class for constants.
  %
  % PROBLEM: How handle spin phase if values are unknown during eclipse, or if
  %          spin phase values jump when exiting eclipse?
  %     """"
  %     * Spin rate is 'unknown' during the eclipse, but do not take care.
  %     * Spin rate is jumped & changed after the eclipse, but do not take care.
  %     """"
  %     /Yasumasa Kasaba, e-mail 2026-06-16, 23:16



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Internal helper function. NOT INTENDED TO BE USED OUTSIDE THIS CLASS.
    %
    % Convert spin phase (0 to 2*pi) to cumulative spin phase (which always
    % increases).
    %
    % NOTE: Assumes that every decrement implies that 2*pi should be added.
    % NOTE: The function assumes that there are no data gaps.
    %
    function cumulSpinPhaseRadAr = spin_phase_to_cumulative_spin_phase(...
        spinPhaseRadAr)

      assert(iscolumn(spinPhaseRadAr) & isa(spinPhaseRadAr, "double"))
      assert(all((0 <= spinPhaseRadAr) & (spinPhaseRadAr <= 2*pi)))
      n = numel(spinPhaseRadAr);

      % IMPLEMENTATION NOTE: unwrap() decrements cumulative spin phase if the
      % spin phase jumps are longer than pi. Therefore not using unwrap().
      cumulSpinPhaseRadAr = NaN(n, 1);
      if n >= 1
        nRevol = 0;
        cumulSpinPhaseRadAr(1) = spinPhaseRadAr(1);
        for i = 2:n

          if spinPhaseRadAr(i-1) > spinPhaseRadAr(i)
            nRevol = nRevol + 1;
          end

          cumulSpinPhaseRadAr(i) = spinPhaseRadAr(i) + 2*pi*nRevol;
        end
      end

      assert(issorted(cumulSpinPhaseRadAr))
    end



    % Internal helper function.  NOT INTENDED TO BE USED OUTSIDE THIS CLASS.
    %
    % Convert cumulative spin phase to (true) TT2000 using interpolation from
    % tabulated known conversions of cumulative spin phase to/from TT2000.
    %
    %
    % ARGUMENTS
    % =========
    % dataTt2000Ar
    %       Column array of TT2000 values
    % dataCumulSpinPhaseRadAr
    %       Column array of known cumulative spin phase values for the
    %       dataTt2000Ar values.
    % inCumulSpinPhaseRadAr
    %       Column array of cumulative spin phase values for which TT2000 shall
    %       be derived.
    %
    function outTt2000Ar = cumulative_spin_phase_to_TT2000( ...
        dataTt2000Ar, dataCumulSpinPhaseRadAr, inCumulSpinPhaseRadAr)

      assert(iscolumn(dataTt2000Ar)            & isa(dataTt2000Ar,            "int64" ))
      assert(iscolumn(dataCumulSpinPhaseRadAr) & isa(dataCumulSpinPhaseRadAr, "double"))
      assert(iscolumn(inCumulSpinPhaseRadAr)   & isa(inCumulSpinPhaseRadAr,   "double"))

      assert(numel(dataTt2000Ar) == numel(dataCumulSpinPhaseRadAr))
      n = numel(dataTt2000Ar);

      % NOTE: Technically, arrays do not need to be sorted (ascending), but the
      % data must still describe a monotonic function for interpolation to work,
      % i.e. if one permutes the elements the same way for both arrays, and so
      % that one of the arrays is sorted, then the other must also become
      % sorted. Otherwise interpolation does not work.
      assert(issorted(dataTt2000Ar,            "STRICTASCEND"))
      assert(issorted(dataCumulSpinPhaseRadAr, "STRICTASCEND"))

      if n >= 2
        % --------------------------------------------------------------
        % CASE: There is enough data for interpolation and extrapolation
        % --------------------------------------------------------------
        % NOTE: interp1() returns double. It returns and NaN if it can not
        % interpolate.
        y = interp1(...
          dataCumulSpinPhaseRadAr, double(dataTt2000Ar), inCumulSpinPhaseRadAr, ...
          "LINEAR", "extrap");
        assert(all(~isnan(y)))
        outTt2000Ar = int64(y);
      else
        % ----------------------------------------------------------
        % CASE: There is NO data for interpolation and extrapolation
        % ----------------------------------------------------------
        % Only permit execution if no actual interpolation/extrapolation is
        % requested.
        assert(...
          isempty(inCumulSpinPhaseRadAr), ...
          "Trying to interpolate/extrapolate when there are fewer than two data points.")
        outTt2000Ar = int64.empty(0, 1);
      end
    end



    % Internal helper function.  NOT INTENDED TO BE USED OUTSIDE THIS CLASS.
    %
    % Given timestamps, identify segments of timestamps which do not increase
    % more than a specified threshold.
    %
    function [iBeginAr, iEndAr, nSegments] = split_by_data_gap(...
        tt2000Ar, dataGapMinNs)

      % PROPOSAL: Move to some "utils" package.
      % PROPOSAL: Relax argument assertions which are not really needed.
      %   Ex: int64
      %     NOTE: Must then relax units/variable types in variable names.
      %       Ex: TT2000, nanoseconds
      % PROPOSAL: Return table.

      assert(iscolumn(tt2000Ar)     & isa(tt2000Ar, 'int64') & issorted(tt2000Ar, "STRICTASCEND"))
      assert(isscalar(dataGapMinNs) & isa(tt2000Ar, 'int64') & dataGapMinNs >= 0)

      if isempty(tt2000Ar)
        iBeginAr   = double.empty(0, 1);
        iEndAr     = double.empty(0, 1);
      else
        boundaryAr = find(diff(tt2000Ar) >= dataGapMinNs);
        iBeginAr   = [1; boundaryAr+1];
        iEndAr     = [boundaryAr; numel(tt2000Ar)];
      end
      nSegments  = numel(iBeginAr);
    end



    % Do spin fit assuming that fit windows used should have constant length and
    % period in spin phase (SAFW).
    %
    % Implemented as a wrapper around bepic.spinfit.fit_TAFW() which adds
    % functionality for splitting processing into smaller time segments based on
    % time increments exceeding a threshold.
    %
    %
    % IMPLEMENTATION NOTE
    % ===================
    % Implements spin-aligned fit windows using mms_spin_fit() by internally
    % using fake TT2000 timestamps derived from cumulative spin phase. This in
    % turn causes problems with handling data gaps when converting from spin
    % phase (0 to 2*pi) to cumulative spin phase, and when then constructing the
    % output timestamps in the presence of data gaps. It therefore makes sense
    % for the implementation to split by data gap before using this approach.
    %
    %
    % POTENTIAL PROBLEM / BUG
    % =======================
    % The function could possibly (theoretically) generate the same timestamps
    % twice if identifying a data gap within a fit window. It is unclear what
    % is the best way to handle such a situation. However,
    % fit_TAFW()'s functionality for removing output timestamps
    % outside the range of the input timestamps should eliminate this
    % possibility().
    %
    %
    % NAME-VALUE ARGUMENTS
    % ====================
    % tt2000Ar
    %       Column array of TT2000 timestamps for every sample.
    % spinPhaseRadAr
    %       Column array of spin phase values. Radians (0 to 2*pi). Must be
    %       finite.
    % samplesAr
    %       Column array of sample values.
    % fitWindowPeriodRad
    %       Length of time between the beginning of each fit window. Radians.
    % fitWindowLengthRad
    %       Length of fit window. Radians.
    % fitWindowCenterRad
    %       Scalar value. Describes where the center of fit windows (output
    %       timestamps) should be in cumulative spin phase. Any time
    %       window center will be located at
    %       fitWindowCenterRad + n * fitWindowPeriodRad, n=integer.
    % nMinFitSamples
    %       Minimum number of samples required for a fit.
    % nFitCoefficients
    %       Number of fit coefficients to use: 3 or 5.
    % dataGapMinNs
    %       Threshold for when a jump in tt2000Ar should count as a data gap.
    %
    %
    % RETURN VALUES
    % =============
    % R ("Result")
    %       Table with ~self-explanatory column names. One spin fit per row.
    %       NOTE: Will never return timestamps outside the interval of input
    %             timestamps.
    %
    function R = fit_SAFW(A)
      % PROPOSAL: Remove data if the same output timestamp is generate twice.
      %   PROBLEM: The timestamps might only be approximately equal?

      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowCenterRad
        A.nMinFitSamples
        A.nFitCoefficients
        A.dataGapMinNs
      end

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(A.tt2000Ar)           & isa(A.tt2000Ar,           "int64"))
      assert(iscolumn(A.spinPhaseRadAr)     & isa(A.spinPhaseRadAr,     "double"))
      assert(iscolumn(A.samplesAr)          & isa(A.samplesAr,          "double"))
      assert(isscalar(A.fitWindowPeriodRad) & isa(A.fitWindowPeriodRad, "double"))
      assert(isscalar(A.fitWindowLengthRad) & isa(A.fitWindowLengthRad, "double"))
      assert(isscalar(A.fitWindowCenterRad) & isa(A.fitWindowCenterRad, "double"))
      assert(isscalar(A.dataGapMinNs)       & isa(A.dataGapMinNs,       "int64"))
      %
      nIn = numel(A.tt2000Ar);
      assert(nIn == numel(A.spinPhaseRadAr))
      assert(nIn == numel(A.samplesAr))
      %
      assert(issorted(A.tt2000Ar, "STRICTASCEND"))
      assert(all(isfinite(A.spinPhaseRadAr)))
      % IMPLEMENTATION NOTE: Requiring spin phase on interval 0 to 2*pi. This is
      % not required by mms_spin_fit(), but (1) is (analogous to) the spin phase
      % values in L1p CDF files which are also bounded (0 to 360 deg), and (2)
      % is an indirect check on the units used.
      assert(all((0 <= A.spinPhaseRadAr) & (A.spinPhaseRadAr <= 2*pi)))
      assert(A.fitWindowPeriodRad > 0)
      assert(A.fitWindowLengthRad > 0)
      %
      assert(A.dataGapMinNs > 0)

      % ========================================================================
      % Convert from "spin phase radians" to fake TT2000/duration so that values
      % can be fed to bepic.spinfit.fit_TAFW()
      % ========================================================================
      % Fake nanoseconds per radian when converting to/from fake TT2000.
      N = 4e9 / (2*pi);

      % IMPLEMENTATION NOTE: Cumulative spin phase values will not increment
      % correctly for time jumps (error n*2*pi) but that does not matter, since
      % the processing will be split by data gaps anyway.
      cumulSpinPhaseRadAr = bepic.spinfit.spin_phase_to_cumulative_spin_phase(...
        A.spinPhaseRadAr);
      fakeTt2000Ar              = int64(cumulSpinPhaseRadAr  * N);
      fakeFitWindowPeriodNs     = int64(A.fitWindowPeriodRad * N);
      fakeFitWindowLengthNs     = int64(A.fitWindowLengthRad * N);
      fakeFitWindowCenterTt2000 = int64(A.fitWindowCenterRad * N);

      nSamples = numel(A.tt2000Ar);
      if nSamples == 0
        % IMPLEMENTATION NOTE: Call bepic.spinfit.fit_TAFW() with
        % empty data, just to create a consistent return value.
        R = bepic.spinfit.fit_TAFW( ...
          tt2000Ar              = fakeTt2000Ar, ...
          spinPhaseRadAr        = A.spinPhaseRadAr, ...
          samplesAr             = A.samplesAr, ...
          fitWindowPeriodNs     = fakeFitWindowPeriodNs, ...
          fitWindowLengthNs     = fakeFitWindowLengthNs, ...
          fitWindowCenterTt2000 = fakeFitWindowCenterTt2000, ...
          nMinFitSamples        = A.nMinFitSamples, ...
          nFitCoefficients      = A.nFitCoefficients);
      else
        % =========================================================
        % Identify indices defining the beginning and end of a time
        % jump-separated segment.
        % =========================================================
        [iBeginAr, iEndAr, nSegments] = bepic.spinfit.split_by_data_gap(...
          A.tt2000Ar, A.dataGapMinNs);

        rCa = cell(nSegments, 1);
        for i = 1:nSegments
          iAr = iBeginAr(i):iEndAr(i);

          rSegment = bepic.spinfit.fit_TAFW( ...
            tt2000Ar              = fakeTt2000Ar    (iAr), ...
            spinPhaseRadAr        = A.spinPhaseRadAr(iAr), ...
            samplesAr             = A.samplesAr     (iAr), ...
            fitWindowPeriodNs     = fakeFitWindowPeriodNs, ...
            fitWindowLengthNs     = fakeFitWindowLengthNs, ...
            fitWindowCenterTt2000 = fakeFitWindowCenterTt2000, ...
            nMinFitSamples        = A.nMinFitSamples, ...
            nFitCoefficients      = A.nFitCoefficients);

          % ======================================================
          % Modify the timestamps, from fake TT2000 to true TT2000
          % ======================================================
          % IMPLEMENTATION NOTE: Must do this separately for every call to
          % bepic.spinfit.fit_TAFW() in to correctly handle spin
          % phase values which are (legitimately) identical just before and
          % after a data gap.
          outCumulSpinPhaseRadAr = double(rSegment.fitWindowCenterTt2000Ar) / N;
          rSegment.fitWindowCenterTt2000Ar = bepic.spinfit.cumulative_spin_phase_to_TT2000(...
            A.tt2000Ar         (iAr), ...
            cumulSpinPhaseRadAr(iAr), ...
            outCumulSpinPhaseRadAr);

          rCa{i, 1} = rSegment;
        end

        R = vertcat(rCa{:});
      end

    end



    % Do spin fit assuming that fit windows should have constant length and
    % period in time (TAFW).
    %
    % Imlemented as reusable wrapper around mms_spinfit_m() for
    % (1) enforcing strict input arguments and return values:
    %   (1a) MATLAB classes (data types),
    %   (1b) array sizes
    % (2) automatically deriving "t0" so that a constant argument can be used
    %     for this function.
    % (3) better naming of input arguments and return values.
    % (4) hardcoding certain argument(s).
    %
    %
    % NAME-VALUE ARGUMENTS
    % ====================
    % tt2000Ar
    %       Column array of TT2000 timestamps for every sample.
    % spinPhaseRadAr
    %       Column array of spin phase values. Radians (0 to 2*pi). Must be
    %       finite.
    % samplesAr
    %       Column array of sample values.
    % fitWindowPeriodNs.
    %       Length of time between the beginning of each fit window.
    %       Nanoseconds.
    % fitWindowLengthNs.
    %       Length of fit window. Nanoseconds.
    % fitWindowCenterTt2000
    %       Scalar value. Describes where the center of fit windows (output
    %       timestamps) should be in time. Any fit window center will be
    %       located at fitWindowCenterTt2000 + n * fitWindowPeriodNs, n=integer.
    % nMinFitSamples
    %       Minimum number of samples required for a fit.
    % nFitCoefficients
    %       Number of fit coefficients to use: 3 or 5.
    %
    %
    % RETURN VALUES
    % =============
    % R ("Result")
    %       Table with ~self-explanatory column names. One spin fit per row.
    %       NOTE: Will never return timestamps outside the interval of input
    %             timestamps.
    %
    function R = fit_TAFW(A)
      % PROPOSAL: Expose constants as arguments?
      %   Ex: N_MIN_REQUIRED_FIT_SAMPLES.
      %   CON: Need to write more tests.
      %     PRO: Might need to test for behaviour which is not used.
      %   PRO: Easier to change constants (update implementation)

      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.fitWindowPeriodNs
        A.fitWindowLengthNs
        A.fitWindowCenterTt2000
        A.nMinFitSamples
        A.nFitCoefficients
      end

      N_MAX_FIT_ITERATIONS = 5;              % TODO: Determine proper value.

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(A.tt2000Ar)              & isa(A.tt2000Ar,              "int64"))
      assert(iscolumn(A.spinPhaseRadAr)        & isa(A.spinPhaseRadAr,        "double"))
      assert(iscolumn(A.samplesAr)             & isa(A.samplesAr,             "double"))
      assert(isscalar(A.fitWindowPeriodNs)     & isa(A.fitWindowPeriodNs,     "int64"))
      assert(isscalar(A.fitWindowLengthNs)     & isa(A.fitWindowLengthNs,     "int64"))
      assert(isscalar(A.fitWindowCenterTt2000) & isa(A.fitWindowCenterTt2000, "int64"))
      assert(isscalar(A.nMinFitSamples))
      assert(isscalar(A.nFitCoefficients))
      %
      nIn = numel(A.tt2000Ar);
      assert(nIn == numel(A.spinPhaseRadAr))
      assert(nIn == numel(A.samplesAr))
      %
      assert(issorted(A.tt2000Ar, "STRICTASCEND"))
      assert(all(isfinite(A.spinPhaseRadAr)))
      % IMPLEMENTATION NOTE: Requiring spin phase on interval 0 to 2*pi. This is
      % not required by mms_spin_fit(), but (1) is (analogous to) the spin phase
      % values in L1p CDF files which are also bounded (0 to 360 deg), and (2)
      % is an indirect check on the units used.
      assert(all((0 <= A.spinPhaseRadAr) & (A.spinPhaseRadAr <= 2*pi)))
      %
      assert(A.fitWindowPeriodNs > 0)
      assert(A.fitWindowLengthNs > 0)
      assert(ismember(A.nFitCoefficients, [3, 5]))   % NOTE: Does not support 7!

      % -----------------------------------------
      % DOCUMENTATION COPIED FROM mms_spinfit_m()
      % -----------------------------------------
      % function [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(maxIt, minPts, nTerms, timeData, data, phase, fitEvery, fitInterv, t0)
      %  Compute spinfit coefficients to spinning data. Data is fitted to
      %  function y = A + Bcos(phase) + Csin(phase) + (Dcos(2*phase) +
      %  Esin(2*phase) + Fcos(3*phase) + Gsin(3*phase)). According to the number
      %  of terms specified.
      %
      % Input: (all required)
      %   maxIt    - maximum of iterations for each fit
      %   minPts   - minimum number of points required for each fit
      %   nTerms   - number of terms to fit, must be odd (3, 5, 7)
      %   timeData - time of measurement (int64 TT2000)
      %   data     - data (samples) to be fitted
      %   phase    - phase of instrument at corresponding time of measurement
      %              (radians)
      %   fitEvery - one spinfit every X ns (default every 5*10^9 ns)
      %   fitInter - spinfit is fitted to data during this interval
      %              (default 20*10^9 ns)
      %   t0       - the first time inside timeData which is evenly divisable
      %              with fitEvery, accounting for leap seconds and such.
      %              (With default, each fit should line up with times 00:00:05,
      %              00:00:10, 00:00:15 etc), (int64 TT2000).
      % Output: (all required)
      %   timeFit  - middle of each spinfit ( 00:00:05, 00:00:10 etc). (int64 TT2000)
      %   sfit     - matrix with each fit coefficents
      %   sdev     - standard deviation of each fit
      %   iter     - number of iterations used for each fit
      %   nBad     - number of bad points, outliers for each fit
      %
      % Bad fits will have value NaN.
      %
      % This is an interface function used by Matlab to display help and/or
      % hints, the real processing occurs in mms_spinfit_mx (mex file).
      % ------------------------------------------------------------------------
      if nIn == 0
        % ========================
        % CASE: Empty input arrays
        % ========================
        % IMPLEMENTATION NOTE: Using empty arrays when calling mms_spinfit_m() can
        % crash MATLAB!!! ("25.2.0.3042426 (R2025b) Update 1", Linux). Can therefore
        % not call it for this case.
        % Ex: [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 5+1, 5, int64.empty(0, 1), double.empty(0, 1), double.empty(0, 1), 4e9, 4e9, 0)

        % Create the equivalent of return values from mms_spinfit_m(), but
        % empty.
        timeFit = int64.empty( 0, 1);
        sfit    = double.empty(0, A.nFitCoefficients);
        sdev    = double.empty(0, 1);
        iter    = double.empty(0, 1);
        nBad    = double.empty(0, 1);
      else
        % ============================
        % CASE: Non-empty input arrays
        % ============================

        % ------------------------------------------------------------
        % Modify A.fitWindowCenterTt2000 to work with mms_spinfit_m()
        % ------------------------------------------------------------
        % IMPLEMENTATION NOTE: mms_spinfit_m() requires "t0" to be within or
        % close to the submitted timestamps but is unclear what this exactly
        % means. If it is not, it might crash or add (not NaN) or omit return
        % values for timestamps for fit windows which there are no samples.
        m = idivide(...
          A.tt2000Ar(1) - A.fitWindowCenterTt2000, ...
          A.fitWindowPeriodNs, "FLOOR");
        modifFitWindowCenterTt2000 = A.fitWindowCenterTt2000 + m * A.fitWindowPeriodNs;

        % --------------------
        % CALL mms_spinfit_m()
        % --------------------
        % Manual test code for mms_spinfit_m(): One-liner which can be
        % copy-pasted to the MATLAB command-line for manual experimentation.
        % n=9; i=1:n; tt2000Ar=int64(linspace(0.1e9, 7.8e9, n)); spinRad=double(tt2000Ar)/4e9; ; samplesAr=3+4*cos(spinRad)+5*sin(spinRad); [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(50, 3+1, 3, tt2000Ar(i), samplesAr(i), spinRad(i), 4e9, 4e9, 2e9)
        %
        % IMPORTANT NOTE: mms_spinfit_m() behaves strangely for
        % N_MIN_REQUIRED_FIT_SAMPLES=4, N_FIT_TERMS=3. In practice, 4 samples
        % per fit window does not appear to be enough to produce non-NaN fits
        % coefficients.
        % Ex: n=25*4+24; i=1:n; tt2000Ar=int64(linspace(0.1e9, 99.8e9, n)); spinRad=double(tt2000Ar)/4e9; ; samplesAr=3+4*cos(spinRad)+5*sin(spinRad); [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 3+1, 3, tt2000Ar(i), samplesAr(i), spinRad(i), 4e9, 4e9, 2e9)
        [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(...
          N_MAX_FIT_ITERATIONS, A.nMinFitSamples, A.nFitCoefficients, ...
          A.tt2000Ar, A.samplesAr, A.spinPhaseRadAr, ...
          A.fitWindowPeriodNs, A.fitWindowLengthNs, ...
          modifFitWindowCenterTt2000);
      end

      % =========================================
      % Assertions on mms_spinfit_m return values
      % =========================================
      % IMPLEMENTATION NOTE: Useful for (1) testing the understanding/behaviour
      % of mms_spinfit_m(), and (2) ensuring that emulated mms_spinfit_m()
      % return values are consistent with the non-emulated ones.
      assert(isa(timeFit, "int64"))
      nOut = numel(timeFit);
      assert(isequal(size(timeFit), [nOut, 1]))
      assert(isequal(size(sfit), [nOut, A.nFitCoefficients]))
      assert(isequal(size(sdev), [nOut, 1]))
      assert(isequal(size(iter), [nOut, 1]))
      assert(isequal(size(nBad), [nOut, 1]))

      % ==========================================
      % Determine which output timestamps too keep
      % ==========================================
      % IMPLEMENTATION NOTE: It is unclear exactly which output timestamps
      % mms_spinfit_m() will return. It can vary depending on the choice of
      % "t0", and timestamps can be both added and removed. This code tries to
      % handle this problem by removing timestamps outside of input data and
      % hopefully make the function more "deterministic".
      if nOut >= 1
        bKeep = (A.tt2000Ar(1) <= timeFit) & (timeFit <= A.tt2000Ar(end));
      else
        bKeep = logical.empty(0, 1);
      end

      % =========================
      % Construct return value(s)
      % =========================
      R = table();
      R.fitWindowCenterTt2000Ar = timeFit(bKeep, 1);
      R.stdDeviationAr          = sdev(   bKeep);
      R.nBadPoints              = nBad(   bKeep);
      R.offsetAr                = sfit(   bKeep, 1);
      R.coefficientCos1Ar       = sfit(   bKeep, 2);
      R.coefficientSin1Ar       = sfit(   bKeep, 3);
      if A.nFitCoefficients >= 5
        R.coefficientCos2Ar     = sfit(   bKeep, 4);
        R.coefficientSin2Ar     = sfit(   bKeep, 5);
      end
    end



  end    % methods(Static)



end
