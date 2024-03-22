%
% Code relating to generating test data for tests.
%
%
% NOTE: Including "MANUAL TESTS" for generating quicklooks based on test data
% (not CDFs). This speeds up the generation of plots which is useful for
% manually inspecting quicklooks. This code will not be called by runtests().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef testdata
  % PROPOSAL: Move manual tests to separate file(s).
  %   PROPOSAL: One class
  %   manual tests
  %   ~mantests, ~manual_tests



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    % zVariable "Energy" in solo_L2_swa-pas-eflux which appears to be constant
    % metadata between datasets.
    %
    % zVariable attributes:
    % CATDESC  = "Center of energy bins",
    % VAR_TYPE = "support_data".
    % No DEPEND_0, so not time-dependent.
    SWA_PAS_EFLUX_ENERGY = [
      18450.500000, 17398.822266, 16407.089844, 15471.887695, 14589.988281, ...
      13758.356445, 12974.128906, 12234.612305, 11537.236328, 10879.615234, ...
      10259.477539, 9674.684570, 9123.227539, 8603.206055, 8112.821777, ...
      7650.391113, 7214.318359, 6803.105957, 6415.328125, 6049.655273, ...
      5704.823730, 5379.648926, 5073.010742, 4783.850098, 4511.167480, ...
      4254.028320, 4011.552734, 3782.895996, 3567.269043, 3363.936768, ...
      3172.194092, 2991.376953, 2820.865479, 2660.074707, 2508.449707, ...
      2365.471436, 2230.639404, 2103.489258, 1983.590942, 1870.529663, ...
      1763.910400, 1663.363770, 1568.549561, 1479.143066, 1394.834473, ...
      1315.328613, 1240.356079, 1169.656250, 1102.985107, 1040.117065, ...
      980.832703, 924.921570, 872.199524, 822.485840, 775.600952, 731.390076, ...
      689.703064, 650.389893, 613.315918, 578.361084, 545.395081, 514.303528, ...
      484.991272, 457.347992, 431.279358, 406.694763, 383.509888, 361.649261, ...
      341.038269, 321.601471, 303.269409, 285.981812, 269.678864, 254.310608, ...
      239.816696, 226.142548, 213.253281, 201.098328, 189.633362, 178.828094, ...
      168.637360, 159.021439, 149.955460, 141.409119, 133.347549, 125.745872, ...
      118.578949, 111.821907, 105.449593, 99.437134, 93.769684, 88.426933, ...
      83.384018, 78.631226, 74.153107, 69.924812 ...
      ]';

    % Return struct field "f" from read_TNR().
    % Effectively zVariable "TNR_BAND_FREQ", below 100 kHz. (solo.read_TNR()
    % filters out higher frequencies/values).
    % Ex: solo_L2_rpw-tnr-surv-cdag_20230101_V08.cdf
    TNR_F = uint32([
      3992, 4169, 4353, 4546, 4747, 4957, 5177, 5406, 5645, 5895, 6156, 6429, ...
      6713, 7011, 7321, 7645, 7984, 8337, 8706, 9092, 9494, 9914, 10353, 10812, ...
      11290, 11790, 12312, 12857, 13427, 14021, 14642, 15290, 15967, 16674, ...
      17412, 18183, 18988, 19829, 20707, 21624, 22581, 23581, 24625, 25715, ...
      26853, 28042, 29284, 30580, 31934, 33348, 34825, 36366, 37976, 39658, ...
      41414, 43247, 45162, 47161, 49249, 51430, 53707, 56085, 58568, 61161, ...
      63869, 66696, 69649, 72733, 75953, 79316, 82827, 86494, 90324, 94323, ...
      98499 ...
      ])'
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Generate plots for manual inspection from test data.
    function manual_24h_6h_2h_test(outputDir)
      close all
      tic

      QuicklooksTint = irf.tint(...
        '2024-01-01T00:00:00.00Z', ...
        '2024-01-02T00:00:00.00Z');
      SpacePosTint = irf.tint(...
        '2023-12-01T00:00:00.00Z', ...
        '2024-02-02T00:00:00.00Z');

      % Shorter time interval to speed up test, but not too short so that it is
      % hard to inspect the relevant part of the plot (spectrum).
      BTint = irf.tint(...
        '2024-01-01T01:00:00', ...
        '2024-01-01T01:30:00');

      OutputPaths = solo.qli.utils.create_output_directories(outputDir);
      Data        = solo.qli.testdata.generate_test_data(...
        QuicklooksTint, SpacePosTint, BTint);

      irfLogoPath = solo.qli.testdata.get_test_logo_path();
      %irfLogoPath = '/nonhome_data/work_files/SOLAR_ORBITER/irfu-matlab_qli/mission/solar_orbiter/+solo/+qli/+offgen/irf_logo.png';

      solo.qli.generate_quicklooks_24h_6h_2h(Data, OutputPaths, QuicklooksTint, irfLogoPath)
      toc
    end



    % Generate plot for manual inspection from test data.
    function manual_7days_test(outputDir)
      close all
      tic

      QuicklooksTint = irf.tint(...
        '2023-12-27T00:00:00.00Z', ...
        '2024-01-03T00:00:00.00Z');
      SpacePosTint = irf.tint(...
        '2023-12-01T00:00:00.00Z', ...
        '2024-02-02T00:00:00.00Z');
      % Shorter time interval to speed up test, but not too short so that it is
      % hard to inspect the relevant part of the plot (spectrum).
      BTint = QuicklooksTint;

      Data = solo.qli.testdata.generate_test_data(...
        QuicklooksTint, SpacePosTint, BTint);

      irfLogoPath = solo.qli.testdata.get_test_logo_path();

      solo.qli.generate_quicklook_7days(Data, outputDir, QuicklooksTint, irfLogoPath)
      toc
    end



    % Generate test data which can be used as input to
    % solo.qli.generate_quicklooks_24h_6h_2h() and
    % solo.qli.generate_quicklook_7days().
    function Data = generate_empty_test_data(SpacePosTint)
      % NOTE: Data.Vrpw is loaded as a TSeries from .mat files and time interval
      % is then selected with TSeries.tlim(). ==> Empty data is represented as a
      % length-zero TSeries.

      Data = struct();
      Data.Vrpw     = TSeries(EpochTT(zeros(0, 1)), zeros(0, 1));
      Data.E        = [];
      Data.Ne       = [];
      Data.B        = [];
      Data.Tpas     = [];
      Data.Vpas     = [];
      Data.Npas     = [];
      Data.ieflux   = [];
      Data.tnrBand  = [];

      % NOTE: Tnr can be
      % [],     if there are no datasets for the specified date)
      % 0,      if some condition is met in solo.read_TNR()
      % struct, if there is data (half-truth; struct can represent absence of data too)
      Data.Tnr      = [];

      Data.swaEnergyMetadata = solo.qli.testdata.SWA_PAS_EFLUX_ENERGY;

      Data.soloPos  = solo.qli.testdata.get_space_position(SpacePosTint);
      Data.earthPos = solo.qli.testdata.get_space_position(SpacePosTint);
    end



    % Generate test data which can be used as input to
    % solo.qli.generate_quicklooks_24h_6h_2h() and
    % solo.qli.generate_quicklook_7days().
    %
    % ARGUMENTS
    % =========
    % QuicklooksTint
    %       Time interval for data in general (not data referred to in other
    %       arguments).
    % SpacePosTint
    %       Time interval which will be filled with simulated SPICE data for
    %       SolO and Earth positions.
    % BTint
    %       Time interval which will be filled with B data. This is useful for
    %       specifying the amount of B data and thereby how slow the plotting
    %       should be.
    %
    function Data = generate_test_data(QuicklooksTint, SpacePosTint, BTint)

      % Variables for supporting 1 h time resolution data.
      % n1h = Number of timestamps.
      [T1h, n1h] = solo.qli.testdata.get_time_array(QuicklooksTint, 60*60);
      %f1h = linspace(0, 1, n1h)';

      % Variables for supporting 1 min. time resolution data
      [T1m, n1m] = solo.qli.testdata.get_time_array(QuicklooksTint, 60);
      f1m = linspace(0, 1, n1m)';


      %====================
      % Generate test data
      %====================
      Data = struct();



      % Vrpw = VHT data
      % ---------------
      % NOTE: There is only VHT data up until
      %       2020-12-05T22:00:00.000000000Z (1h) and
      %       2020-12-05T18:00:00.000000000Z (6h).
      %       ==> Vrpw is often an empty (zero elements) TSeries.
      % NOTE: Only returning data with one time resolution for simplicity. In
      % principle, 24h/6h/2h and 7day quicklooks use 2h and 6h time resolution
      % respectively.
      % NOTE: Individual values are plotted as circles, i.e. one wants FEW
      % values.
      % Ex: Day with Vrpw data 2020-06-03.
      V_RPW_DATA = linspace(-779.9284, -36.0898, n1h)';
      Data.Vrpw  = irf.ts_scalar(T1h, double(V_RPW_DATA));
      %Data.Vrpw  = irf.ts_scalar(T1h([], 1), zeros(0, 1));



      Data.E  = irf.ts_vec_xyz(T1m, single([NaN(n1m, 1), f1m, -f1m]));
      Data.Ne = irf.ts_scalar( T1m, single(0.1 + f1m*10));



      % NOTE: Having a large amount of B data slows down the plotting a lot in
      %       the call to irf_ebsp(). Therefore only want to use a short time
      %       period with B data.
      % NOTE: B is also used for calculating the proton gyration frequency which is
      %       plotted in the same panel as the B spectrum.
      % NOTE: A badly chosen sampling rate for B can generate an error: "F_MAX
      %       must be lower than the Nyquist frequency". Must therefore generate
      %       data with an accepted sampling frequency.
      % NOTE: The sampling frequency is taken from B data for 2023-01-01
      %       (1/7.9997=median(diff(Data.B.time.tts))).
      B_SAMPLING_RATE_HZ = 1/7.9997;
      [TB, nB] = solo.qli.testdata.get_time_array(BTint, B_SAMPLING_RATE_HZ);
      fB       = linspace(0, 1, nB)';   % f=fraction
      B_DATA   = [sin(fB*5), cos(fB*4), 2*fB] .* (1+10*fB);
      % This test data works, but it makes |B| fluctuate a lot.
      % B_DATA   = solo.qli.testdata.get_2D_array(nB, 3, -10, 10, 'LIN');
      Data.B   = irf.ts_vec_xyz(TB, single(B_DATA));



      Data.Tpas = irf.ts_scalar( T1m, single(f1m*100));
      Data.Vpas = irf.ts_vec_xyz(T1m, single([f1m*1000, f1m*100, 90*f1m]));



      NPAS_DATA = linspace(7.2675, 23.5056, n1m)';
      Npas      = irf.ts_scalar(T1m, single(NPAS_DATA));
      Data.Npas = Npas;



      % solo_L2_swa-pas-eflux_20230101_V02.cdf:
      %   min(Data.ieflux.data(:)) = 0
      %   max(Data.ieflux.data(:)) = 3.1242e+09
      %   numel(find(Data.ieflux.data(:) == 0)) = 1386744
      %   numel(find(Data.ieflux.data(:) ~= 0)) = 663336
      %   NOTE: Many zeros for something that is plotted on a logarithmic scale.
      N_IEFLUX_COLUMNS = 96;
      IEFLUX_MAX       = 3.1242e+09;
      % NOTE: Deliberately including zeros for something that is plotted on a
      % logarithmic scale.
      IEFLUX_DATA = solo.qli.testdata.get_2D_array(n1m, N_IEFLUX_COLUMNS, 0, IEFLUX_MAX, 'LIN');
      % TODO-NI: Is there some better irf.* function for initializing a TSeries?
      Data.ieflux = TSeries(T1m, single(IEFLUX_DATA));



      % NOTE: Should only be non-empty when "Tnr" is non-empty.
      % Should ideally be higher time resolution.
      % Values are arbitrary (except range), but may in principle indicate some
      % setting. The value is only checked for emptiness anyway.
      Data.tnrBand = irf.ts_scalar(...
        T1m, ...
        uint8(solo.qli.testdata.get_2D_array(n1m, 1, 0, 3, 'LIN')));



      % Ex: solo_L2_rpw-tnr-surv-cdag_20230101_V08.cdf
      % t: [21505×1 double]
      % f: [75×1 uint32]
      % p: [21505×75 double]
      % p_label: {'dB'}
      % min(Data.Tnr.p, [], 'all') = 1.3639e-154
      % max(Data.Tnr.p, [], 'all') = 2.4172e-90
      Tnr   = struct();
      Tnr.t = T1m.epochUnix;    % EpochTT.epochUnix
      Tnr.f = solo.qli.testdata.TNR_F;
      % NOTE: p values run over many orders of magnitude and are plotted as a
      % spectrum. Therefore important to spread out values "exponentially" or
      % else the spectrum will visually appear to be constant value (even if it
      % is not).
      P_MIN = 1.3639e-154;
      P_MAX = 2.4172e-90;
      Tnr.p = solo.qli.testdata.get_2D_array(...
        n1m, length(solo.qli.testdata.TNR_F), ...
        P_MIN, P_MAX, 'LOG');
      Tnr.p_label = {'dB'};
      Data.Tnr    = Tnr;



      Data.swaEnergyMetadata = solo.qli.testdata.SWA_PAS_EFLUX_ENERGY;



      % Data.Etnr    = TSeries(Times, uint8(mod([1:n]', 4)));
      Data.soloPos  = solo.qli.testdata.get_space_position(SpacePosTint);
      Data.earthPos = solo.qli.testdata.get_space_position(SpacePosTint);
    end



    % Return synthetic TSeries for the position of SolO/Earth for a specified
    % time interval.
    %
    % NOTE: The tested code extracts the correct time interval of data itself,
    % as opposed to for science data itself. Tests should thus ideally NOT
    % specify the same time range as the science data.
    function SoloPosTSeries = get_space_position(TotalTint)
      assert(isa(TotalTint, 'GenericTimeArray'))

      DT_SEC = 60*60;
      UNITS  = irf_units;
      AU_KM  = UNITS.AU / UNITS.km;    % Astronomical unit [km]

      LengthSec = TotalTint.stop - TotalTint.start;    % double, seconds
      nPts      = round(LengthSec / DT_SEC) + 1;

      SoloPosTint = irf.time_array(TotalTint.start, double(linspace(0, LengthSec, nPts)));
      f = linspace(0, 1, nPts)';
      % Columns: [soloSunDistance, soloEclLongitude, soloEclLatitude]
      solPosArray = [AU_KM*ones(size(f)), 180*f, 180*f];

      SoloPosTSeries = TSeries(SoloPosTint, solPosArray);
    end



    % Return 2D array with test data. Uses sin() for variation.
    % Is designed to behave nicely also for 1D arrays.
    %
    % ARGUMENTS
    % =========
    function A = get_2D_array(nRows, nCols, aMin, aMax, scale)
      % For debugging.
      %       function mm(a)
      %         m1 = min(a, [], 'all');
      %         m2 = max(a, [], 'all');
      %         [m1, m2]
      %       end

      % IMPLEMENTATION NOTE: Using multiplication cos()*cos() to make sure to
      % always have full variation (from min to max) on the edges. This gives
      % good behaviour for 1D arrays.

      assert(aMin < aMax)
      assert(0 <= nRows)
      assert(0 <= nCols)

      N_PERIODS_1 = 2;
      N_PERIODS_2 = 2;

      % Values in range 0-1.
      f1 = linspace(0, 1, nRows);
      f2 = linspace(0, 1, nCols);

      % Values in range -0.5 -- +0.5
      g1 = cos(f1*N_PERIODS_1*2*pi);
      g2 = cos(f2*N_PERIODS_2*2*pi);

      % Values in range -0.5 -- +0.5
      F1 = repmat(g1', 1, nCols);
      F2 = repmat(g2,  nRows, 1);

      % Value in range 0-1.
      F = 0.5 + 0.5 * (F1.*F2);

      switch(scale)
        case 'LIN'
          A = aMin + F * (aMax-aMin);
        case 'LOG'
          assert(aMin > 0)
          A = exp( log(aMin) + F * (log(aMax)-log(aMin)) );
        otherwise
          error('Illegal argument scale=%s', scale)
      end

      % Rounding errors can make values go sligtly outside limits. Remove such
      % "overflow".
      A = max(A, aMin);
      A = min(A, aMax);
    end



    % Create image file which can be used as test logo.
    function create_test_logo(filePath)
      % NOTE: Argument determines the file type, which is not perfect.

      % NOTE: irf_logo.png is 1220x1226 pixels.
      TEST_LOGO = [
        0,1,0;
        1,0,1;
        0,1,0
        ];
      imwrite(TEST_LOGO, filePath)
    end



    function testLogoPath = get_test_logo_path()
      testLogoPath = fullfile(tempdir, 'test_logo.png');
      solo.qli.testdata.create_test_logo(testLogoPath)
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Return GenericTimeArray for time interval. Useful for generating series of
    % timestamps for synthetic data.
    %
    % ARGUMENTS
    % =========
    % Syntax 1: utsStr1, utcStr2, dtSec
    % Syntax 2: Tint, dtSec
    function [T, n] = get_time_array(varargin)
      switch(nargin)
        case 2
          [Tint, dtSec]             = varargin{:};
        case 3
          [utsStr1, utcStr2, dtSec] = varargin{:};
          Tint                      = irf.tint(utsStr1, utcStr2);
        otherwise
          error('Illegal number of arguments')
      end

      lengthSec  = Tint.stop - Tint.start;
      dtSecArray = 0:dtSec:lengthSec;
      T          = Tint.start + dtSecArray;
      n          = length(T);
    end



  end    % methods(Static, Access=private)



end
