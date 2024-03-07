%
% matlab.unittest automatic test code for solo.qli.context_info_strings().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef context_info_strings___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      % Arbitrary number output variables.
      function test(...
          soloPosTSeries, earthPosTSeries, Tint, ...
          expSoloStr, expEarthStr)

        [actSoloStr, actEarthStr] = solo.qli.context_info_strings(soloPosTSeries, earthPosTSeries, Tint);
        testCase.verifyEqual(actSoloStr,  expSoloStr)
        testCase.verifyEqual(actEarthStr, expEarthStr)
      end

      %===================================================================

      Units = irf_units;
      AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

      ETT = EpochTT( ...
        [ ...
        '2024-01-10T00:00:00.000000000Z'; ...
        '2024-01-11T00:00:00.000000000Z'; ...
        '2024-01-12T00:00:00.000000000Z' ...
        ] ...
        );
      POSITION_TS = TSeries( ...
        ETT, ...
        [ ...
        1*AU_KM, 3, 4; ...
        2*AU_KM, 5, 6; ...
        3*AU_KM, 7, 8; ...
        ] ...
        );

      % In-range time interval.
      TI_1 = EpochTT(['2024-01-09T00:00:00.000000000Z'; '2024-03-13T00:00:00.000000000Z']);
      test(...
        POSITION_TS, POSITION_TS, TI_1, ...
        'SolO:  1.00 AU,  EcLat 229\circ,  EcLon 172\circ', ...
        'Earth:  EcLon 172\circ')

      % Out-of-range time interval.
      TI_2 = EpochTT(['2024-01-01T00:00:00.000000000Z'; '2024-01-02T00:00:00.000000000Z']);
      test(...
        POSITION_TS, POSITION_TS, TI_2, ...
        '', ...
        '')
    end



  end    % methods(Test)



end
