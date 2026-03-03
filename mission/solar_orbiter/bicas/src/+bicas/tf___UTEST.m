%
% matlab.unittest automatic test code for bicas.tf.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-10
%
classdef tf___UTEST < matlab.unittest.TestCase



  properties(TestParameter)
    METHOD = {'FFT', 'KERNEL'}
  end



  properties (Constant)
    NON_FV_SPLIT_SETTINGS = struct(...
      'snfEnabled',          false, ...
      'snfSubseqMinSamples', 1 ...
      )
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Enable RE-trending without DE-trending. ==> Error
    %
    function test_apply_TF___illegal_detrending(T, METHOD)

      N  = 100;
      dt = 0.1;
      y1 = 5 * ones(N, 1);
      tf = bicas.tf.utest_utils.get_TF_constant(29, 29);    % Constant TF.

      T.verifyError(...
        @() (bicas.tf.apply_TF(...
        dt, y1, tf, ...
        snfEnabled          = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
        snfSubseqMinSamples = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
        method              = METHOD, ...
        detrendingDegreeOf  = -10, ...
        retrendingEnabled   = 1 ...   % Should be bool/logical.
        )), ...
        ?MException)
    end



    % Zero-order detrending. Constant signal ==> Output=0
    %
    function test_apply_TF___detrending0(T, METHOD)

      N  = 100;
      dt = 0.1;
      y1 = 5 * ones(N, 1);
      tf = bicas.tf.utest_utils.get_TF_constant(29, 29);    % Constant TF.

      [y2, D] = bicas.tf.apply_TF(...
        dt, y1, tf, ...
        snfEnabled              = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
        snfSubseqMinSamples     = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
        method                  = METHOD, ...
        detrendingDegreeOf      = 0, ...
        retrendingEnabled       = false ...
        );

      % All three tests: No tolerance works on irony but fails in GitHub CI.
      T.verifyEqual(numel(D.y1ModifCa), 1)
      T.verifyEqual(D.y1ModifCa{1}, 0*y1, 'AbsTol', 1e-14)
      T.verifyEqual(D.y2ModifCa{1}, 0*y1, 'AbsTol', 1e-13)
      T.verifyEqual(y2,             0*y1, 'AbsTol', 1e-13)
    end



    % Test one signal with different parts being removed depending of degree
    % of detrending, and re-trending enabled/disabled.
    %
    function test_apply_TF___detrending_parts(T, METHOD)

      N  = 100;
      dt = 0.1;
      x  = linspace(-1, 1, N)';   % "Normalized" time.

      A  = 5;
      B  = 2;
      % NOTE: y1_0 removed exactly by zero-order de-trending.
      y1_0 = A * ones(size(x));
      y1_1 = B * x;
      y1   = y1_0 + y1_1;
      tf   = bicas.tf.utest_utils.get_TF_constant(29, 29);    % Constant TF.



      % No detrending
      [y2, D] = bicas.tf.apply_TF(...
        dt, y1, tf, ...
        snfEnabled              = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
        snfSubseqMinSamples     = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
        method                  = METHOD, ...
        detrendingDegreeOf      = -1, ...
        retrendingEnabled       = false);
      T.verifyEqual(D.y1ModifCa{1}, y1,    'AbsTol', 1e-14)
      T.verifyEqual(D.y2ModifCa{1}, y1*29, 'AbsTol', 1e-13)
      T.verifyEqual(y2,             y1*29, 'AbsTol', 1e-13)



      % Zero order detrending.
      [y2, D] = bicas.tf.apply_TF(...
        dt, y1, tf, ...
        snfEnabled              = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
        snfSubseqMinSamples     = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
        method                  = METHOD, ...
        detrendingDegreeOf      = 0, ...
        retrendingEnabled       = false);

      T.verifyEqual(D.y1ModifCa{1}, y1_1,    'AbsTol', 1e-14)
      T.verifyEqual(D.y2ModifCa{1}, y1_1*29, 'AbsTol', 1e-13)
      T.verifyEqual(y2,             y1_1*29, 'AbsTol', 1e-13)



      % Second order detrending.
      [y2, D] = bicas.tf.apply_TF(...
        dt, y1, tf, ...
        snfEnabled              = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
        snfSubseqMinSamples     = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
        method                  = METHOD, ...
        detrendingDegreeOf      = 2, ...
        retrendingEnabled       = false);

      T.verifyEqual(D.y1ModifCa{1}, zeros(size(y1_1)), 'AbsTol', 1e-14)
      T.verifyEqual(D.y2ModifCa{1}, zeros(size(y1_1)), 'AbsTol', 1e-13)
      T.verifyEqual(y2,             zeros(size(y1_1)), 'AbsTol', 1e-13)



      % Second order detrending + RE-trending
      [y2, D] = bicas.tf.apply_TF(...
        dt, y1, tf, ...
        snfEnabled              = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
        snfSubseqMinSamples     = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
        method                  = METHOD, ...
        detrendingDegreeOf      = 2, ...
        retrendingEnabled       = true);

      T.verifyEqual(D.y1ModifCa{1}, zeros(size(y1_1)), 'AbsTol', 1e-14)
      T.verifyEqual(D.y2ModifCa{1}, zeros(size(y1_1)), 'AbsTol', 1e-13)
      T.verifyEqual(y2,  29*y1,             'AbsTol', 1e-13)
    end



    % Test (Nyquist) frequency cutoff.
    %
    % NOTE: Obsolete since argument "tfHighFreqLimitFraction" and corresponding
    % functionality has been removed.
    %
    % function test_apply_TF___freq_cutoff(T)
    %
    %   N  = 2^7;
    %   dt = 0.1;
    %   t  = [0:N-1]' * dt;
    %
    %   nyquistOmegaRps = pi/dt;
    %   omega1 = nyquistOmegaRps*0.25;   % Survives   tfHighFreqLimitFraction.
    %   omega2 = nyquistOmegaRps*0.50;   % Removed by tfHighFreqLimitFraction.
    %   y1_1   = sin(omega1*t);
    %   y1_2   = sin(omega2*t);
    %   y1     = y1_1 + y1_2;
    %
    %   tf     = bicas.tf.utest_utils.get_TF_delay(1*dt);
    %
    %   if 1
    %
    %     [y2, D] = bicas.tf.apply_TF(...
    %       dt, y1, tf, ...
    %       snfEnabled              = T.NON_FV_SPLIT_SETTINGS.snfEnabled, ...
    %       snfSubseqMinSamples     = T.NON_FV_SPLIT_SETTINGS.snfSubseqMinSamples, ...
    %       method                  = 'FFT', ...
    %       detrendingDegreeOf      = -1, ...
    %       retrendingEnabled       = false, ...
    %       tfHighFreqLimitFraction = 0.4 ...
    %       );
    %
    %     y2_exp = circshift(y1_1, 1);   % Requires FFT method.
    %     %bicas.tf.apply_TF___UTEST.plot_test(y1, y2, y2_exp)
    %
    %     T.verifyEqual(abs(D.tfModif(omega1)), 1)
    %     T.verifyEqual(abs(D.tfModif(omega2)), 0)
    %     T.verifyEqual(y2, y2_exp, 'AbsTol', 1e-13)
    %   end
    % end



    function test_apply_TF___SNF(T)
      % NOTE: Tests wrt. snfEnabled are overkill, but exist for historical
      % reasons. This is mostly tested via the corresponding function
      % bicas.tf.split_samples_by_nonfinite(). Testing wrt. snfSubseqMinSamples
      % is still relevant.

      C = 3;

      function test(S)
        % IMPLEMENTATION NOTE: Using struct "S" as argument to get around
        % MATLAB fobidding name-value arguments in nested functions.

        % expY2ModifCa = cell(1,0);
        % for j = 1:numel(S.expY1ModifCa)
        %   expY2ModifCa{j, 1} = S.expY1ModifCa{j} * C;
        % end

        dt = 0.1;    % Value should be irrelevant.
        tf = bicas.tf.utest_utils.get_TF_constant(C, C);

        [actY2, ActD] = bicas.tf.apply_TF(...
          dt, S.y1, tf, ...
          method              = 'KERNEL', ...
          detrendingDegreeOf  = -1, ...
          retrendingEnabled   = false, ...
          snfEnabled          = true, ...
          snfSubseqMinSamples = S.snfSubseqMinSamples ...
          );

        % T.assertEqual(ActD.y1ModifCa, S.expY1ModifCa, 'RelTol', 1e-14)
        % T.assertEqual(ActD.y2ModifCa,   expY2ModifCa, 'RelTol', 1e-14)
        T.assertEqual(actY2,          S.expY2,        'RelTol', 1e-14)
      end

      test(...
        struct( ...
        'snfSubseqMinSamples', 1, ...
        'y1',           [1:100]', ...
        'expY2',        [1:100]' * C ...
        ));

      test(...
        struct( ...
        'snfSubseqMinSamples', 1, ...
        'y1',           [1:10 NaN 12:100]', ...
        'expY2',        [1:10 NaN 12:100]' * C ...
        ));

      test(...
        struct( ...
        'snfSubseqMinSamples', 1, ...
        'y1',           [NaN 2:10 Inf 12:100]', ...
        'expY2',        [NaN 2:10 NaN 12:100]' * C ...
        ));

      test(...
        struct( ...
        'snfSubseqMinSamples', 1, ...
        'y1',           [1:90 NaN 92:99 NaN]', ...
        'expY2',        [1:90 NaN 92:99 NaN]' * C ...
        ));

      % No splitting. Subsequence too short.
      test(...
        struct( ...
        'snfSubseqMinSamples', 1000, ...
        'y1',           [1:100]', ...
        'expY2',        [NaN(1, 100)]' * C ...
        ));

      % Splitting. 1/2 subsequences is too short.
      test(struct( ...
        'snfSubseqMinSamples', 20, ...
        'y1',           [NaN 2:10 NaN 12:100]', ...
        'expY2',        [NaN(1,11)    12:100]' * C ...
        ))
    end



    function test_split_samples_by_nonfinite(T)
      function test(y, expYCa)
        %============================
        % Argument consistency check
        %============================
        % Implicitly consistency check on the tested function.
        % --
        % IMPLEMENTATION NOTE: cell2mat() is ~inconsistent and must therefore
        % be handles using a special case. cell2mat(cell(0, 1)) ==> 0x0
        % NOTE: The empty cell2mat() return value can not contain any MATLAB
        % class corresponding to that of "y".
        if isempty(y)
          y2 = zeros(0, 1);
        else
          y2 = cell2mat(expYCa);
        end
        assert(isequaln(y, y2))



        actYCa = bicas.tf.split_samples_by_nonfinite(y);

        % Actual return value format check.
        assert(iscolumn(actYCa))
        for i = 1:numel(actYCa)
          assert(iscolumn(actYCa{i}));
        end

        T.verifyEqual(actYCa, expYCa)
      end

      test(zeros(0, 1), cell(0, 1))

      test([3]', {3}')
      test( ...
        [ 3 4 5]', ...
        {[3 4 5]'}')

      test( ...
        [ inf nan -inf]', ...
        {[inf nan -inf]'}')

      test(...
        [3 inf 5 -inf 7 nan 9]', ...
        {3 inf 5 -inf 7 nan 9}')
      test(...
        [inf  2 3    nan inf -inf     4 5 6 7]', ...
        {inf [2 3]' [nan inf -inf]', [4 5 6 7]'}')
    end



    function test_make_hard_low_pass_TF(T)
      tf       = bicas.tf.utest_utils.get_TF_constant(3, 3);
      dtSec    = 0.1;
      fraction = 0.7;

      actTf = bicas.tf.make_hard_low_pass_TF(tf, fraction, dtSec);

      nyquistFreqRps = pi / dtSec;   % 0.5/dtSec [Hz] = pi/dtSec [rad/s]
      omegaRps       = linspace(0, 100);
      expZ           = tf(omegaRps);
      expZ(omegaRps > fraction*nyquistFreqRps) = 0;

      % Verify configuration
      assert(ismember(0, expZ))
      assert(ismember(3, expZ))

      actZ = actTf(omegaRps);

      T.assertEqual(actZ, expZ)
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function plot_test(y1, y2_act, y2_exp)
      figure
      plot([y1, y2_act, y2_exp+0.01], '*-')
      legend({'y1', 'y2_{act}', 'y2_{exp}'})
    end



  end    % methods(Static, Access=private)



end
