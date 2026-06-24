%
% Functions which implement the core processing for MEFISTO spin fits, when
% processing L1p CDFs --> L2pre CDFs.
%
% The purpose of this file is to isolate many of the technical details
% surrounding the details of the spin fit, e.g. the time windows used for spin
% fitting, and which exact timestamps should represent spin fit time windows.
%
% The implementation uses mms_spinfit_m for the core processing. In the long
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
  %   ~spin fit "engine", "core"
  %
  % PROPOSAL: Convert class to package dedicated to official processing.
  %   CON: Can not put constants in class.
  %     CON: Can create dedicated class for constants.
  %
  % TODO-DEC: Naming convention for spin fitting functions
  %   ~diff, E field
  %   Type of time windows: constant-length, spin phase-aligned
  %
  % PROBLEM: Unclear how one should handle data gaps, in particular for
  %          spin-aligned time windows?
  %   Ex: Large jums in time indicate that spin phase values (wrapped) can not
  %       be used for determining cumulative spin phase across the time jump.
  %   Ex: If a time window contains a jump in time, contains
  %       sufficiently many samples for making a spin fit, but the center spin
  %       phase point is on the other side of the time jump, how should one then
  %       calculate the center spin phase point as TT2000?
  %     PROPOSAL: Do fit on spin phase and extrapolate.
  %       CON: mms_spin_fit() does the grouping of samples and spin phases.
  %            Therefore does not have direct access to groups of values on
  %            which to make the fit.
  %       CON: mms_spin_fit() was built because doing so many fits was
  %            time-consuming. Therefore, this might be time consuming too.
  %     PROPOSAL: Should detect data gaps and split up the processing
  %               into separate chunks before calling mms_spin_fit().



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant)
    % Constants for when using constant-length time windows.
    TIME_WINDOW_PERIOD_NS     = int64(4e9);
    TIME_WINDOW_LENGTH_NS     = int64(4e9);
    TIME_WINDOW_BEGIN_REFERENCE_TT2000 = int64(0e9);
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Internal helper function. Convert spin phase (0 to 2*pi) to cumulative
    % spin phase.
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



    % Internal helper function. Convert cumulative spin phase to spin phase (0
    % to 2*pi radians).
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
      % that one of the arrays is sorted, then the other must also become sorted.
      % Otherwise interpolation does not work.
      assert(issorted(dataTt2000Ar,            "strictascend"))
      assert(issorted(dataCumulSpinPhaseRadAr, "strictascend"))

      if n >= 2
        % --------------------------------------------------------------
        % CASE: There is enough data for interpolation and extrapolation
        % --------------------------------------------------------------
        % NOTE: interp1() returns double. It returns and NaN if it can not
        % interpolate.
        y = interp1(...
          dataCumulSpinPhaseRadAr, double(dataTt2000Ar), inCumulSpinPhaseRadAr, ...
          "linear", "extrap");
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



    % Wrapper around bepic.spinfit.mms_spinfit_wrapper() which
    % (1) adds functionality for splitting processing into smaller time segments
    %     based on time increments exceeding a threshold.
    % (2) works with spin-aligned time windows.
    %
    %
    % IMPLEMENTATION NOTE
    % ===================
    % Implements spin-aligned time windows using mms_spin_fit() by using fake
    % TT2000 timestamps based on cumulative spin phase. This in causes problems
    % with handling data gaps when converting from spin phase (0 to 2*pi) to
    % cumulative spin phase, and when constructing the output timestamps in the
    % presence of data gaps. It therefore makes sense for the implementation to
    % split by data gap before using this approach.
    %
    % POTENTIAL PROBLEM / BUG
    % =======================
    % Function should be able to generate the same timestamps twice if choosing
    % a dataGapMinNs which is smaller than the time window (in time). It is
    % unclear what is the best way to handle such a situation.
    %
    function R = spin_aligned(A)
      % PROPOSAL: Better name
      %
      % TODO: Automated tests.
      %
      % PROPOSAL: Remove data when the same output timestamp is generate twice.
      %   PROBLEM: The timestamps might only be approximately equal?

      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.timeWindowPeriodRad
        A.timeWindowLengthRad
        A.timeWindowReferenceRad
        A.dataGapMinNs
      end

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(A.tt2000Ar)        & isa(A.tt2000Ar,       "int64")  & issorted(A.tt2000Ar, 'strictascend'))
      assert(iscolumn(A.spinPhaseRadAr)  & isa(A.spinPhaseRadAr, "double") & all(isfinite(A.spinPhaseRadAr)))
      % IMPLEMENTATION NOTE: Requiring spin phase on interval 0 to 2*pi. This is
      % not required by mms_spin_fit(), but (1) is (analogous to) the spin phase
      % values in L1p CDF files which are also bounded (0 to 360 deg), and (2)
      % is an indirect check on the units used.
      assert(all((0 <= A.spinPhaseRadAr) & (A.spinPhaseRadAr <= 2*pi)))
      %
      assert(iscolumn(A.samplesAr)       & isa(A.samplesAr,  "double"))
      nIn = numel(A.tt2000Ar);
      assert(nIn == numel(A.spinPhaseRadAr))
      assert(nIn == numel(A.samplesAr))
      %
      assert(isscalar(A.timeWindowPeriodRad) & (A.timeWindowPeriodRad > 0) & isa(A.timeWindowPeriodRad,    "double"))
      assert(isscalar(A.timeWindowLengthRad) & (A.timeWindowLengthRad > 0) & isa(A.timeWindowLengthRad,    "double"))
      assert(isscalar(A.timeWindowReferenceRad)                            & isa(A.timeWindowReferenceRad, "double"))
      %
      assert(isscalar(A.dataGapMinNs) & (A.dataGapMinNs > 0) & isa(A.dataGapMinNs, 'int64'))

      % Fake nanoseconds per radian when converting to/from fake TT2000.
      N = 4e9 / 2*pi;


      % ========================================================================
      % Convert from "spin phase radians" to fake TT2000/duration so that values
      % can be fed to bepic.spinfit.mms_spinfit_wrapper()
      % ========================================================================
      % IMPLEMENTATION NOTE: Cumulative spin phase values will not increment
      % correctly for time jumps (error n*2*pi) but that does not matter, since
      % the processing will be split by data gaps anyway.
      cumulSpinPhaseRadAr = bepic.spinfit.spin_phase_to_cumulative_spin_phase(...
        A.spinPhaseRadAr);
      fakeTt2000Ar              = int64(cumulSpinPhaseRadAr      * N);
      timeWindowPeriodNs        = int64(A.timeWindowPeriodRad    * N);
      timeWindowLengthNs        = int64(A.timeWindowLengthRad    * N);
      timeWindowReferenceTt2000 = int64(A.timeWindowReferenceRad * N);

      nSamples = numel(A.tt2000Ar);
      if nSamples == 0
        % IMPLEMENTATION NOTE: Call bepic.spinfit.mms_spinfit_wrapper() with empty
        % data, just to create a consistent return value.
        R = bepic.spinfit.mms_spinfit_wrapper( ...
          tt2000Ar                  = fakeTt2000Ar, ...
          spinPhaseRadAr            = A.spinPhaseRadAr, ...
          samplesAr                 = A.samplesAr, ...
          timeWindowPeriodNs        = timeWindowPeriodNs, ...
          timeWindowLengthNs        = timeWindowLengthNs, ...
          timeWindowReferenceTt2000 = timeWindowReferenceTt2000);
      else
        % =========================================================
        % Identify indices defining the beginning and end of a time
        % jump-separated segment.
        % =========================================================
        boundaryAr = diff(A.tt2000Ar);
        iBeginAr   = [1; boundaryAr];
        iEndAr     = [boundaryAr+1; nSamples];
        nSegments  = numel(iBeginAr);

        rCa = cell(nSegments, 1);
        for i = 1:nSegments
          iAr = iBeginAr(i):iEndAr(i);

          rCa{i, 1} = bepic.spinfit.mms_spinfit_wrapper( ...
            tt2000Ar                  = A.tt2000Ar                 (iAr), ...
            spinPhaseRadAr            = A.spinPhaseRadAr           (iAr), ...
            samplesAr                 = A.samplesAr                (iAr), ...
            timeWindowPeriodNs        = A.timeWindowPeriodNs       (iAr), ...
            timeWindowLengthNs        = A.timeWindowLengthNs       (iAr), ...
            timeWindowReferenceTt2000 = A.timeWindowReferenceTt2000(iAr));
        end

        R = vertcat(rCa{:});
      end

      % ======================================================
      % Modify the timestamps, from fake TT2000 to true TT2000
      % ======================================================
      outCumulSpinPhaseRadAr = double(R.tt2000Ar) / N;
      R.tt2000Ar = bepic.spinfit.cumulative_spin_phase_to_TT2000(...
        A.tt2000Ar, cumulSpinPhaseRadAr, outCumulSpinPhaseRadAr);
    end



    % Reusable wrapper around mms_spinfit_m() for
    % (1) enforcing strict input arguments and return values.
    %   (1a) MATLAB classes (data types),
    %   (1b) array sizes
    % (2) better naming of input arguments and return values.
    %
    %
    % ARGUMENTS
    % =========
    % tt2000Ar
    %       TT2000 timestamps for samples.
    % spinPhaseRadAr
    %       Spin phase values
    % samplesAr
    %       Sample values.
    % --
    % All arguments are same-length column arrays.
    %
    %
    % RETURN VALUES
    % =============
    % Struct with self-explanatory fieldnames. One spin fit per row.
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function R = mms_spinfit_wrapper(A)
      % PROPOSAL: Expose constants as arguments?
      %   Ex: N_MIN_REQUIRED_FIT_SAMPLES.
      %   CON: Need to write more tests.
      %     PRO: Might need to test for behaviour which is not used.
      %   PRO: Easier to change constants (update implementation)
      %
      % TODO-DEC: Range of spin phase? Infinite? 0 to 2*pi.
      %   PROPOSAL: Assert spin phase interval: 0 to 2*pi; -pi to +pi
      %     PRO: Checks radians vs degrees

      % Manual test code for mms_spinfit_m():
      % n=10; spinRad=linspace(0, 0.1*pi, n)'; [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(5, 5+1, 3+2, int64([1:n]'), 3+cos(spinRad)+2*sin(spinRad), spinRad, 4e9, 4e9, 0)

      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.timeWindowPeriodNs
        A.timeWindowLengthNs
        A.timeWindowReferenceTt2000
      end

      N_FIT_TERMS                = 3+2;
      N_MAX_FIT_ITERATIONS       = 5;              % TODO: Determine proper value.
      N_MIN_REQUIRED_FIT_SAMPLES = N_FIT_TERMS+3;  % TODO: Determine proper value.

      tt2000Ar                  = A.tt2000Ar;
      spinPhaseRadAr            = A.spinPhaseRadAr;
      samplesAr                 = A.samplesAr;
      timeWindowPeriodNs        = A.timeWindowPeriodNs;
      timeWindowLengthNs        = A.timeWindowLengthNs;
      timeWindowReferenceTt2000 = A.timeWindowReferenceTt2000;

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(tt2000Ar)        & isa(tt2000Ar,       "int64")  & issorted(tt2000Ar, 'strictascend'))
      assert(iscolumn(spinPhaseRadAr)  & isa(spinPhaseRadAr, "double") & all(isfinite(spinPhaseRadAr)))
      % IMPLEMENTATION NOTE: Requiring spin phase on interval 0 to 2*pi. This is
      % not required by mms_spin_fit(), but (1) is (analogous to) the spin phase
      % values in L1p CDF files which are also bounded (0 to 360 deg), and (2)
      % is an indirect check on the units used.
      assert(all((0 <= spinPhaseRadAr) & (spinPhaseRadAr <= 2*pi)))
      %
      assert(iscolumn(samplesAr)       & isa(samplesAr,      "double"))
      nIn = numel(tt2000Ar);
      assert(nIn == numel(spinPhaseRadAr))
      assert(nIn == numel(samplesAr))
      %
      assert(isscalar(timeWindowPeriodNs) & (timeWindowPeriodNs > 0) & isa(timeWindowPeriodNs,        "int64"))
      assert(isscalar(timeWindowLengthNs) & (timeWindowLengthNs > 0) & isa(timeWindowLengthNs,        "int64"))
      assert(isscalar(timeWindowReferenceTt2000)                     & isa(timeWindowReferenceTt2000, "int64"))

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

        timeFit = int64.empty( 0, 1);
        sfit    = double.empty(0, 5);
        sdev    = double.empty(0, 1);
        iter    = double.empty(0, 1);
        nBad    = double.empty(0, 1);
      else
        % ============================
        % CASE: Non-empty input arrays
        % ============================
        [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(...
          N_MAX_FIT_ITERATIONS, N_MIN_REQUIRED_FIT_SAMPLES, N_FIT_TERMS, ...
          tt2000Ar, samplesAr, spinPhaseRadAr, ...
          timeWindowPeriodNs, timeWindowLengthNs, ...
          timeWindowReferenceTt2000);
      end

      % =========================================
      % Assertions on mms_spinfit_m return values
      % =========================================
      % IMPLEMENTATION NOTE: Useful for (1) testing the understanding/behaviour
      % of mms_spinfit_m(), and (2) ensuring that emulated mms_spinfit_m()
      % return values are consistent with the non-emulated ones.
      assert(isa(timeFit, "int64"))
      nOut = numel(timeFit);
      assert(iscolumn(timeFit) & (numel(timeFit) == nOut))
      assert(isequal(size(sfit), [nOut, N_FIT_TERMS]))
      assert(iscolumn(sdev)    & (numel(sdev)    == nOut))
      assert(iscolumn(iter)    & (numel(iter)    == nOut))
      assert(iscolumn(nBad)    & (numel(nBad)    == nOut))

      % =========================
      % Construct return value(s)
      % =========================
      R = struct();
      R.tt2000Ar          = timeFit;
      R.offsetAr          = sfit(:, 1);
      R.coefficientCos1Ar = sfit(:, 2);
      R.coefficientSin1Ar = sfit(:, 3);
      R.coefficientCos2Ar = sfit(:, 4);
      R.coefficientSin2Ar = sfit(:, 5);
    end



  end    % methods(Static)



end
