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

%             % Arbitrary number output variables.
            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));

                [actOutputs{:}] = solo.qli.context_info_strings(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end

            %===================================================================

            Units = irf_units;
            AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

            ett = EpochTT( ...
                [ ...
                    '2024-01-10T00:00:00.000000000Z'; ...
                    '2024-01-11T00:00:00.000000000Z'; ...
                    '2024-01-12T00:00:00.000000000Z' ...
                ] ...
            );
            Ts = TSeries( ...
                ett, ...
                [ ...
                    1*AU_KM, 3, 4; ...
                    2*AU_KM, 5, 6; ...
                    3*AU_KM, 7, 8; ...
                ] ...
            );

            % In-range time interval.
            ti1 = EpochTT(['2024-01-09T00:00:00.000000000Z'; '2024-03-13T00:00:00.000000000Z']);
            test({Ts, Ts, ti1}, {'SolO:  1.00 AU,  EcLat 229\circ,  EcLon 172\circ', 'Earth:  EcLon 172\circ'})

            % Out-of-range time interval.
            ti2 = EpochTT(['2024-01-01T00:00:00.000000000Z'; '2024-01-02T00:00:00.000000000Z']);
            test({Ts, Ts, ti2}, {'', ''})

        end



    end    % methods(Test)

end
