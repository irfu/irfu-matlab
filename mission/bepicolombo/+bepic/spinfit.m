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
  %     CON: The mission is already implied by the parent package "bepic".
  %       CON: The instrument is not implied.
  %   ~spin fit "engine", "core"
  %   fit
  %   processing
  %
  % PROPOSAL: Separate classes for functions (1) based on, and (2) not based on
  %           mms_spinfit_m().
  %
  % PROBLEM: How handle spin phase if values are unknown during eclipse, or if
  %          spin phase values jump when exiting eclipse?
  %     """"
  %     * Spin rate is 'unknown' during the eclipse, but do not take care.
  %     * Spin rate is jumped & changed after the eclipse, but do not take care.
  %     """"
  %     /Yasumasa Kasaba, e-mail 2026-06-16, 23:16
  %
  % PROPOSAL: Optional key-value argument in "official functions" specifies
  %           whether debugging data should be returned in additional return
  %           values (default=false).



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    %###################################
    %###################################
    % FUNCTIONS BASED ON mms_spinfit_m()
    %###################################
    %###################################



    % Do spin fit assuming that fit windows used should have both constant
    % length and constant period in units of spin phase (SAFW).
    %
    % Implemented as a wrapper around bepic.spinfit.fit_TAFW_MMS() which adds
    % functionality for splitting processing into smaller time segments based on
    % time increments exceeding a threshold.
    %
    %
    % IMPLEMENTATION NOTE
    % ===================
    % Implements spin-aligned fit windows using mms_spin_fit() by internally
    % generating fake TT2000 timestamps derived from cumulative spin phase. This
    % in turn causes problems with handling data gaps when converting from spin
    % phase (0 to 2*pi) to cumulative spin phase, and when then constructing the
    % output timestamps in the presence of data gaps. It therefore makes sense
    % for the implementation to split by data gap before using this approach.
    %
    %
    % POTENTIAL PROBLEM / BUG
    % =======================
    % The function could possibly (theoretically) generate the same timestamps
    % twice if identifying a data gap within a fit window. It is unclear what is
    % the best way to handle such a situation. However, fit_TAFW_MMS()'s
    % functionality for removing output timestamps outside the range of the
    % input timestamps should eliminate this possibility().
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
    % fitWindowCenterRefRad
    %       Scalar value. Describes where the center of fit windows (output
    %       timestamps) should be in cumulative spin phase. Any time
    %       window center will be located at a phase
    %       fitWindowCenterRefRad + n * fitWindowPeriodRad,
    %       where n=integer.
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
    function R = fit_SAFW_MMS(A)
      % PROPOSAL: Remove data if the same output timestamp is generated twice.
      %   PROBLEM: The timestamps might only be approximately equal?

      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowCenterRefRad
        A.nMinFitSamples
        A.nFitCoefficients
        A.dataGapMinNs
      end

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(A.tt2000Ar)              & isa(A.tt2000Ar,              "int64"))
      assert(iscolumn(A.spinPhaseRadAr)        & isa(A.spinPhaseRadAr,        "double"))
      assert(iscolumn(A.samplesAr)             & isa(A.samplesAr,             "double"))
      assert(isscalar(A.fitWindowPeriodRad)    & isa(A.fitWindowPeriodRad,    "double"))
      assert(isscalar(A.fitWindowLengthRad)    & isa(A.fitWindowLengthRad,    "double"))
      assert(isscalar(A.fitWindowCenterRefRad) & isa(A.fitWindowCenterRefRad, "double"))
      assert(isscalar(A.dataGapMinNs)          & isa(A.dataGapMinNs,          "int64"))
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
      % can be fed to bepic.spinfit.fit_TAFW_MMS()
      % ========================================================================
      % Fake nanoseconds per radian when converting to/from fake TT2000. The
      % exact value should not matter.
      FAKE_NS_PER_RAD = 4e9 / (2*pi);
      % IMPLEMENTATION NOTE: CPS values will not increment correctly for time
      % jumps (error n*2*pi) but that does not matter, since the processing will
      % be split by data gaps anyway.
      cspRadAr = bepic.spinfit.utils.spin_phase_to_CMP(...
        A.spinPhaseRadAr);
      fakeTt2000Ar                 = int64(cspRadAr                * FAKE_NS_PER_RAD);
      fakeFitWindowPeriodNs        = int64(A.fitWindowPeriodRad    * FAKE_NS_PER_RAD);
      fakeFitWindowLengthNs        = int64(A.fitWindowLengthRad    * FAKE_NS_PER_RAD);
      fakeFitWindowCenterRefTt2000 = int64(A.fitWindowCenterRefRad * FAKE_NS_PER_RAD);

      nSamples = numel(A.tt2000Ar);
      if nSamples == 0
        % IMPLEMENTATION NOTE: Call bepic.spinfit.fit_TAFW_MMS() with
        % empty data, just to create a consistent return value.
        R = bepic.spinfit.fit_TAFW_MMS( ...
          tt2000Ar                 = fakeTt2000Ar, ...
          spinPhaseRadAr           = A.spinPhaseRadAr, ...
          samplesAr                = A.samplesAr, ...
          fitWindowPeriodNs        = fakeFitWindowPeriodNs, ...
          fitWindowLengthNs        = fakeFitWindowLengthNs, ...
          fitWindowCenterRefTt2000 = fakeFitWindowCenterRefTt2000, ...
          nMinFitSamples           = A.nMinFitSamples, ...
          nFitCoefficients         = A.nFitCoefficients);
      else
        % =========================================================
        % Identify indices defining the beginning and end of a time
        % jump-separated segment.
        % =========================================================
        [iBeginAr, iEndAr, nSegments] = bepic.spinfit.utils.split_by_data_gap(...
          A.tt2000Ar, A.dataGapMinNs);

        % ====================================
        % Iterate over jump-separated segments
        % ====================================
        rCa = cell(nSegments, 1);
        for i = 1:nSegments
          iAr = iBeginAr(i):iEndAr(i);

          rSegment = bepic.spinfit.fit_TAFW_MMS( ...
            tt2000Ar                 = fakeTt2000Ar    (iAr), ...
            spinPhaseRadAr           = A.spinPhaseRadAr(iAr), ...
            samplesAr                = A.samplesAr     (iAr), ...
            fitWindowPeriodNs        = fakeFitWindowPeriodNs, ...
            fitWindowLengthNs        = fakeFitWindowLengthNs, ...
            fitWindowCenterRefTt2000 = fakeFitWindowCenterRefTt2000, ...
            nMinFitSamples           = A.nMinFitSamples, ...
            nFitCoefficients         = A.nFitCoefficients);

          % ======================================================
          % Modify the timestamps, from fake TT2000 to true TT2000
          % ======================================================
          % IMPLEMENTATION NOTE: Must do this separately for every call to
          % bepic.spinfit.fit_TAFW_MMS() in order to correctly handle spin
          % phase values which are (legitimately) identical just before and
          % after a data gap.
          outCspRadAr = double(rSegment.fitWindowCenterTt2000) / FAKE_NS_PER_RAD;
          rSegment.fitWindowCenterTt2000 = bepic.spinfit.utils.CMP_to_TT2000(...
            A.tt2000Ar(iAr), ...
            cspRadAr  (iAr), ...
            outCspRadAr);

          rCa{i, 1} = rSegment;
        end

        R = vertcat(rCa{:});
      end

    end



    % Do spin fit assuming that fit windows should have constant length and
    % period in units of time (TAFW).
    %
    % Implemented as a reusable wrapper around mms_spinfit_m() which adds
    % functionality for:
    % (1) enforcing (asserting) strict input arguments and return values:
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
    %       Column array of incrementing TT2000 timestamps for every sample.
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
    % fitWindowCenterRefTt2000
    %       Scalar value. Describes where the center of fit windows (output
    %       timestamps) should be in time. Any fit window center will be
    %       located at a time
    %       fitWindowCenterRefTt2000 + n * fitWindowPeriodNs,
    %       where n=integer.
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
    function R = fit_TAFW_MMS(A)
      % PROPOSAL: Expose constant as argument?
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
        A.fitWindowCenterRefTt2000
        A.nMinFitSamples
        A.nFitCoefficients
      end

      N_MAX_FIT_ITERATIONS = 5;    % TODO: Determine proper value.

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(A.tt2000Ar)                 & isa(A.tt2000Ar,                 "int64"))
      assert(iscolumn(A.spinPhaseRadAr)           & isa(A.spinPhaseRadAr,           "double"))
      assert(iscolumn(A.samplesAr)                & isa(A.samplesAr,                "double"))
      assert(isscalar(A.fitWindowPeriodNs)        & isa(A.fitWindowPeriodNs,        "int64"))
      assert(isscalar(A.fitWindowLengthNs)        & isa(A.fitWindowLengthNs,        "int64"))
      assert(isscalar(A.fitWindowCenterRefTt2000) & isa(A.fitWindowCenterRefTt2000, "int64"))
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
      % function [timeFit, sfit, sdev, iter, nBad] = mms_spinfit_m(
      %     maxIt, minPts, nTerms, timeData, data, phase, fitEvery, fitInterv, t0)
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
        % IMPLEMENTATION NOTE: Using empty arrays when calling mms_spinfit_m()
        % can crash MATLAB!!! ("25.2.0.3042426 (R2025b) Update 1", Linux). Can
        % therefore not call it for this case.
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

        % --------------------------------------------------------------
        % Modify A.fitWindowCenterRefTt2000 to work with mms_spinfit_m()
        % --------------------------------------------------------------
        % IMPLEMENTATION NOTE: mms_spinfit_m() requires "t0" to be within or
        % close to the submitted timestamps but is unclear what this exactly
        % means. If it is not, it might (a) crash, or (b) add (not NaN) or omit
        % return values for timestamps for fit windows which there are no
        % samples.
        m = idivide(...
          A.tt2000Ar(1) - A.fitWindowCenterRefTt2000, ...
          A.fitWindowPeriodNs, "FLOOR");
        modifFitWindowCenterTt2000 = A.fitWindowCenterRefTt2000 + m * A.fitWindowPeriodNs;

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
      assert(isequal(size(sfit),    [nOut, A.nFitCoefficients]))
      assert(isequal(size(sdev),    [nOut, 1]))
      assert(isequal(size(iter),    [nOut, 1]))
      assert(isequal(size(nBad),    [nOut, 1]))

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
      R.fitWindowCenterTt2000 = timeFit(bKeep, 1);
      R.stdDeviation          = sdev(   bKeep);
      R.nBadPoints            = nBad(   bKeep);
      R.offset                = sfit(   bKeep, 1);
      R.coefficientCos1       = sfit(   bKeep, 2);
      R.coefficientSin1       = sfit(   bKeep, 3);
      if A.nFitCoefficients >= 5
        R.coefficientCos2     = sfit(   bKeep, 4);
        R.coefficientSin2     = sfit(   bKeep, 5);
      end
    end



    %#######################################
    %#######################################
    % FUNCTIONS NOT BASED ON mms_spinfit_m()
    %#######################################
    %#######################################



    % Given plasma potential samples for MEF1 and MEF2 (from e.g. an L1pre
    % "potential" CDF), derive averages for fit windows (time windows).
    %
    % NOTE: Derives fits for both MEFISTO probes in one function call since
    %       the same timestamps and fit windows should be used for both probes.
    %
    %
    % NAME-VALUE ARGUMENTS
    % ====================
    % tt2000Ar
    %       Column array of incrementing TT2000 timestamps for every sample.
    % spinPhaseRadAr
    %       Column array of spin phase values. Radians (0 to 2*pi). Must be
    %       finite.
    % samplesAr
    %       Nx2 array for MEF1 and MEF2 sample values. Finite or NaN.
    % fitWindowPeriodRad
    %       Length of time between the beginning of each fit window. Radians.
    % fitWindowLengthRad
    %       Length of fit window. Radians.
    % fitWindowCenterRefRad
    %       Scalar value. Describes where the center of fit windows (output
    %       timestamps) should be in cumulative spin phase. Any time
    %       window center will be located at a phase
    %       fitWindowCenterRefRad + n * fitWindowPeriodRad,
    %       where n=integer.
    % dataGapMinNs
    %       Threshold for when a jump in tt2000Ar should count as a data gap.
    % nMinSamplesPerFitWindow
    %       Minimum number of finite sample values per fit window.
    %
    %
    % RETURN VALUE
    % ============
    % FitTable
    %       Table with partly self-explanatory column names. One fit per row.
    %       FitTable.A
    %         Contains the fit coefficients. See implementation for the exact
    %         meaning. Meaning of indices:
    %         FitTable.A(i_timestamp, i_coefficient).
    %
    function FitTable = fit_pot_SAFW(A)
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowCenterRefRad
        A.dataGapMinNs
        A.nMinSamplesPerFitWindow
      end

      % ==========
      % ASSERTIONS
      % ==========
      nSamples = numel(A.tt2000Ar);
      assert(isequal(size(A.samplesAr), [nSamples, 2]) & isa(A.samplesAr,             "double"))
      assert(all(isfinite(A.samplesAr) | isnan(A.samplesAr), "all"))
      assert(isscalar(A.fitWindowCenterRefRad)         & isa(A.fitWindowCenterRefRad, "double"))

      % =========
      % ALGORITHM
      % =========
      FtwnTable = bepic.spinfit.ftwn.get_SAFWs( ...
        tt2000Ar        = A.tt2000Ar, ...
        spinPhaseRadAr  = A.spinPhaseRadAr, ...
        ftwnPeriodRad   = A.fitWindowPeriodRad, ...
        ftwnLengthRad   = A.fitWindowLengthRad, ...
        ftwnBeginRefRad = A.fitWindowCenterRefRad - A.fitWindowLengthRad/2, ...
        dataGapMinNs    = A.dataGapMinNs);

      % Fit according to "Japan" wiki 2026-10-27.
      % https://chs.isee.nagoya-u.ac.jp/scwiki/doku.php?id=miosc:scteam:instrument:efd:ifdefinitiontomefisto:start
      % Y = A[0]*sin(X + A[1]) + A[2]
      f   = @(C, x) (C(1)*sin(x + C(2)) + C(3));
      C_0 = [0, 0, 0];

      [FitTable1, ~] = bepic.spinfit.fit_array( ...
        tt2000Ar           = A.tt2000Ar, ...
        samplesAr          = A.samplesAr(:, 1), ...
        spinPhaseRadAr     = A.spinPhaseRadAr, ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = f, ...
        fittingStartValue  = C_0, ...
        nMinSamplesPerFtwn = A.nMinSamplesPerFitWindow);

      [FitTable2, ~] = bepic.spinfit.fit_array( ...
        tt2000Ar           = A.tt2000Ar, ...
        samplesAr          = A.samplesAr(:, 2), ...
        spinPhaseRadAr     = A.spinPhaseRadAr, ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = f, ...
        fittingStartValue  = C_0, ...
        nMinSamplesPerFtwn = A.nMinSamplesPerFitWindow);

      FitTable = table();
      FitTable.fitWindowCenterTt2000 = FtwnTable.centerTt2000;
      FitTable.A                     = cat(3, ...
        FitTable1.coefficients, ...
        FitTable2.coefficients);
    end



    % Given E field samples (from e.g. an L1pre "E field" CDF), derive spin fit
    % for fit windows (time windows) and return the coefficients of the fits.
    %
    %
    % NAME-VALUE ARGUMENTS
    % ====================
    % tt2000Ar
    %       Column array of incrementing TT2000 timestamps for every sample.
    % spinPhaseRadAr
    %       Column array of spin phase values. Radians (0 to 2*pi). Must be
    %       finite.
    % samplesAr
    %       Nx1 array of electric field samples (in rotating frame).
    %       This should be the difference between MEF1 and MEF2 (or MEF2 and
    %       MEF1).
    % fitWindowPeriodRad
    %       Length of time between the beginning of each fit window. Radians.
    % fitWindowLengthRad
    %       Length of fit window. Radians.
    % fitWindowCenterRefRad
    %       Scalar value. Describes where the center of fit windows (output
    %       timestamps) should be in cumulative spin phase. Any time
    %       window center will be located at a phase
    %       fitWindowCenterRefRad + n * fitWindowPeriodRad,
    %       where n=integer.
    % dataGapMinNs
    %       Threshold for when a jump in tt2000Ar should count as a data gap.
    % nMinSamplesPerFitWindow
    %       Minimum number of finite sample values per fit window.
    %
    %
    % RETURN VALUE
    % ============
    % FitTable
    %       Table with partly self-explanatory column names. One fit per row.
    %       FitTable.A
    %         Contains the fit coefficients. See implementation for the exact
    %         meaning. Meaning of indices:
    %         FitTable.A(i_timestamp, i_coefficient, i_mefisto_probe).
    %
    function [FitTable] = fit_Efield_SAFW(A)
    %function [FitTable, fitAr] = fit_Efield_SAFW(A)    % DEBUG
      arguments
        A.tt2000Ar
        A.spinPhaseRadAr
        A.samplesAr
        A.fitWindowPeriodRad
        A.fitWindowLengthRad
        A.fitWindowCenterRefRad
        A.dataGapMinNs
        A.nMinSamplesPerFitWindow
      end
      % PROPOSAL: for-->parfor (iteration over fit windows).
      %
      % Code for manually experimenting with lsqcurvefit():
      % clear; f=@(A,x) (A(1)+A(2)*x+A(3)*x.^2); xdata=[-10:10]; ydata=[sin(xdata/5)]; [A,resnorm] = lsqcurvefit(f,[0,0,0],xdata,ydata); close all; plot(xdata,ydata,'o'); hold on; x=linspace(-12,12); plot(x,f(A,x))

      % ==========
      % ASSERTIONS
      % ==========
      nSamples = numel(A.tt2000Ar);
      assert(isscalar(A.dataGapMinNs)                  & isa(A.dataGapMinNs,          "int64"))
      assert(isequal(size(A.samplesAr), [nSamples, 1]) & isa(A.samplesAr,             "double"))
      assert(isscalar(A.fitWindowCenterRefRad)         & isa(A.fitWindowCenterRefRad, "double"))

      % %%%%%%%%%
      % ALGORITHM
      % %%%%%%%%%
      FtwnTable = bepic.spinfit.ftwn.get_SAFWs( ...
        tt2000Ar        = A.tt2000Ar, ...
        spinPhaseRadAr  = A.spinPhaseRadAr, ...
        ftwnPeriodRad   = A.fitWindowPeriodRad, ...
        ftwnLengthRad   = A.fitWindowLengthRad, ...
        ftwnBeginRefRad = A.fitWindowCenterRefRad - A.fitWindowLengthRad/2, ...
        dataGapMinNs    = A.dataGapMinNs);

      % =======================================
      % Function to fit the spinning E field to
      % =======================================
      % if 0
      %   % According to mms_spinfit_m().
      %   f = @(C, x) (C(1) + C(2)*cos(x) + C(3)*sin(x) + C(4)*cos(2*x) + C(5)*sin(2*x));
      %   C_0 = [0, 0, 0, 0, 0];
      % end
      % if 0
      %   % According to "Japan" in "BepiColombo Ground Segment (GS) Development
      %   % and Planning".
      %   % NOTE: Originally found on wiki, but is no longer there. Obsolete?!!
      %   f = @(C, x) (C(1)*cos(x-C(2)) + C(3)*x + C(4));
      %   C_0 = [0, 0, 0, 0];    % Start value for iterations.
      % end
      if 1
        % Fit according to "Japan" on wiki 2026-10-27.
        % https://chs.isee.nagoya-u.ac.jp/scwiki/doku.php?id=miosc:scteam:instrument:efd:ifdefinitiontomefisto:start
        % Y = A[0] * cos(X - A[1]) + A[2] * sin(2*X) + A[3] * cos(2*X) + A[4] * X + A[5]
        % NOTE: Resembles MMS spinfit, but is not the same!!!
        f = @(C, x) (C(1)*cos(x - C(2)) + C(3)*sin(2*x) + C(4) * cos(2*x) + C(5) * x + C(6));
        C_0 = [0, 0, 0, 0, 0, 0];    % Start value for iterations.
      end

      [FitTable, fitAr] = bepic.spinfit.fit_array( ...
        tt2000Ar           = A.tt2000Ar, ...
        samplesAr          = A.samplesAr, ...
        spinPhaseRadAr     = A.spinPhaseRadAr, ...
        FtwnTable          = FtwnTable, ...
        fittingFunction    = f, ...
        fittingStartValue  = C_0, ...
        nMinSamplesPerFtwn = A.nMinSamplesPerFitWindow);

      FitTable = renamevars(FitTable, "coefficients", "A");
      FitTable.fitWindowCenterTt2000 = FtwnTable.centerTt2000;

      %FitTable.fitWindowBeginTt2000  = FtwnTable.beginTt2000;   % DEBUG
      %FitTable.fitWindowEndTt2000    = FtwnTable.endTt2000;     % DEBUG
    end



    % Given fit windows and a fitting function, do fit on 1D array of samples,
    % sample=f(spinPhase).
    %
    % RETURN VALUES
    % =============
    % FitTable
    %       Table of fitting coefficients.
    % fitAr
    %       Sample values corresponding to a perfect fit (associated with
    %       the input timestamps).
    %
    %
    % DESIGN INTENT
    % =============
    % This function is intended for sharing the implementation between other
    % functions and for simplifying testing.
    %
    function [FitTable, fitAr] = fit_array(A)
      arguments
        A.tt2000Ar
        A.samplesAr
        A.spinPhaseRadAr
        A.FtwnTable
        A.fittingFunction
        A.fittingStartValue
        A.nMinSamplesPerFtwn
      end
      % IMPLEMENTATION NOTE: This function should not calculate FtwnTable since
      % it needs to be shared for potential data (same timestamps for two 1D
      % arrays of samples).
      % NOTE: FtwnTable as produced by bepic.spinfit.ftwn.get_SAFWs() also
      %       generates the FTWN centers, which is unnecessary for this
      %       function.
      %
      % PROPOSAL: Use previous fit as start value to speed up+help
      %           convergence?!
      % PROPOSAL: for-->parfor (iteration over fit windows).
      % PROPOSAL: Better name for "fitAr". Should be disting from fit
      %           coefficients.
      %   fittedSamplesAr

      % ==========
      % ASSERTIONS
      % ==========
      assert(iscolumn(A.tt2000Ar)           & isa(     A.tt2000Ar,       "int64"))
      assert(issorted(A.tt2000Ar, "STRICTASCEND"))
      assert(iscolumn(A.spinPhaseRadAr)     & isa(     A.spinPhaseRadAr, "double"))
      assert(iscolumn(A.samplesAr)          & isa(     A.samplesAr,      "double"))
      assert(isa(     A.fittingFunction, "function_handle"))
      assert(isrow(   A.fittingStartValue)  & isfloat(  A.fittingStartValue))
      assert(isscalar(A.nMinSamplesPerFtwn) & isnumeric(A.nMinSamplesPerFtwn))
      %
      assert(A.nMinSamplesPerFtwn >= 0)

      % %%%%%%%%%
      % ALGORITHM
      % %%%%%%%%%
      nCoeff = numel( A.fittingStartValue);
      nFtwn  = height(A.FtwnTable);

      FitTable = table();
      FitTable.coefficients = NaN(nFtwn, nCoeff);     % Pre-allocate
      fitAr          = NaN(size(A.samplesAr));

      iFtwnArCa = bepic.spinfit.ftwn.FTWN_time_to_indices(...
        A.tt2000Ar, ...
        A.FtwnTable.beginTt2000, ...
        A.FtwnTable.endTt2000);

      %N_MAX_ITERATIONS = 10;
      %OPTS = optimoptions("lsqcurvefit", "Display", "off", "MaxIterations", N_MAX_ITERATIONS);
      OPTS = optimoptions("lsqcurvefit", "Display", "off");

      for i = 1:nFtwn
        jAr = iFtwnArCa{i};

        ftwnSpinPhaseRadAr1 = A.spinPhaseRadAr(jAr);
        ftwnSamplesAr1      = A.samplesAr     (jAr);

        % Remove NaN from FTWN.
        bFinite             = isfinite(ftwnSamplesAr1);
        ftwnSpinPhaseRadAr2 = ftwnSpinPhaseRadAr1(bFinite);
        ftwnSamplesAr2      = ftwnSamplesAr1(     bFinite);

        nFtwnSamples2       = numel(ftwnSamplesAr2);

        if max(A.nMinSamplesPerFtwn, nCoeff) <= nFtwnSamples2
          % =============
          % Calculate fit
          % =============
          [ftwnCoefficientAr, ~] = lsqcurvefit(...
            A.fittingFunction, A.fittingStartValue, ...
            ftwnSpinPhaseRadAr2, ftwnSamplesAr2, [], [], OPTS);
        else
          ftwnCoefficientAr = NaN(nCoeff, 1);
        end

        % NOTE: Must use spin phase from before removing fill value samples.
        ftwnFitAr = A.fittingFunction(ftwnCoefficientAr, ftwnSpinPhaseRadAr1);

        FitTable.coefficients(i, :) = ftwnCoefficientAr;
        fitAr         (jAr)  = ftwnFitAr;

        % Interpretation of "Averaged residual of sine fit of Ev".
        % https://chs.isee.nagoya-u.ac.jp/scwiki/doku.php?id=miosc:scteam:instrument:efd:ifdefinitiontomefisto:start
        % NOTE: Not the expression I would expect, mathematically.
        %       /Erik P G Johansson 2026-10-26.
        % FitTable.residualAbsDiff(i, :) = mean(abs(ftwnFitAr-ftwnSamplesAr));
      end    % for
    end    % function



  end    % methods(Static)



end
