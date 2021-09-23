%
% matlab.unittest automatic test code for bicas.tf.apply_TF().
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-10
%
classdef apply_TF___UTEST < matlab.unittest.TestCase
    
    
    
    properties(TestParameter)
        METHOD = {'FFT', 'kernel'}
    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        % Enable RE-trending without DE-trending. ==> Error
        %
        function test_Illegal_detrending(testCase, METHOD)
            
            N  = 100;
            dt = 0.1;
            y1 = 5 * ones(N, 1);
            tf = @(omegaRps) (29);
            
            testCase.verifyError(...
                @() (bicas.tf.apply_TF(...
                    dt, y1, tf, ...
                    'method',             METHOD, ...
                    'detrendingDegreeOf', -10, ...
                    'retrendingEnabled',  1 ...
                )), ...
                ?MException)
        end
        
        
        
        % Zero-order detrending. Constant signal ==> Output=0
        %
        function test_detrending0(testCase, METHOD)
            
            N  = 100;
            dt = 0.1;
            y1 = 5 * ones(N, 1);
            tf = @(omegaRps) (29);
            
            [y2, y1B, y2B, tfB] = bicas.tf.apply_TF(...
                dt, y1, tf, ...
                'method',                  METHOD, ...
                'detrendingDegreeOf',      0, ...
                'retrendingEnabled',       false, ...
                'tfHighFreqLimitFraction', Inf);
            
            testCase.verifyEqual(y1B, 0*y1)
            testCase.verifyEqual(y2B, 0*y1)
            testCase.verifyEqual(y2,  0*y1)
        end
        
        
        
        % Test one signal with different parts being removed depending of degree
        % of detrending, and re-trending enabled/disabled.
        %
        function test_detrending_parts(testCase, METHOD)
            
            N  = 100;
            dt = 0.1;
            x  = linspace(-1, 1, N)';   % "Normalized" time.
            
            A  = 5;
            B  = 2;
            % NOTE: y1_0 removed exactly by zero-order de-trending.
            y1_0 = A * ones(size(x));
            y1_1 = B * x;
            y1   = y1_0 + y1_1;
            tf   = @(omegaRps) (29);   % Constant TF.
            
            
            
            % No detrending
            [y2, y1B, y2B, tfB] = bicas.tf.apply_TF(...
                dt, y1, tf, ...
                'method',                  METHOD, ...
                'detrendingDegreeOf',      -1, ...
                'retrendingEnabled',       false, ...
                'tfHighFreqLimitFraction', Inf);
            
            testCase.verifyEqual(y1B, y1,    'AbsTol', 1e-14)
            testCase.verifyEqual(y2B, y1*29, 'AbsTol', 1e-13)
            testCase.verifyEqual(y2,  y1*29, 'AbsTol', 1e-13)
            
            
            
            % Zero order detrending.
            [y2, y1B, y2B, tfB] = bicas.tf.apply_TF(...
                dt, y1, tf, ...
                'method',                  METHOD, ...
                'detrendingDegreeOf',      0, ...
                'retrendingEnabled',       false, ...
                'tfHighFreqLimitFraction', Inf);
            
            testCase.verifyEqual(y1B, y1_1,    'AbsTol', 1e-14)
            testCase.verifyEqual(y2B, y1_1*29, 'AbsTol', 1e-13)
            testCase.verifyEqual(y2,  y1_1*29, 'AbsTol', 1e-13)
            
            
            
            % Second order detrending.
            [y2, y1B, y2B, tfB] = bicas.tf.apply_TF(...
                dt, y1, tf, ...
                'method',                  METHOD, ...
                'detrendingDegreeOf',      2, ...
                'retrendingEnabled',       false, ...
                'tfHighFreqLimitFraction', Inf);
            
            testCase.verifyEqual(y1B, zeros(size(y1_1)), 'AbsTol', 1e-14)
            testCase.verifyEqual(y2B, zeros(size(y1_1)), 'AbsTol', 1e-13)
            testCase.verifyEqual(y2,  zeros(size(y1_1)), 'AbsTol', 1e-13)
            
            
            
            % Second order detrending + RE-trending            
            [y2, y1B, y2B, tfB] = bicas.tf.apply_TF(...
                dt, y1, tf, ...
                'method',                  METHOD, ...
                'detrendingDegreeOf',      2, ...
                'retrendingEnabled',       true, ...
                'tfHighFreqLimitFraction', Inf);
            
            testCase.verifyEqual(y1B, zeros(size(y1_1)), 'AbsTol', 1e-14)
            testCase.verifyEqual(y2B, zeros(size(y1_1)), 'AbsTol', 1e-13)
            testCase.verifyEqual(y2,  29*y1,             'AbsTol', 1e-13)
        end
        
        
        
        % Test (Nyquist) frequency cutoff.
        %
        function test_Freq_cutoff(testCase)
            %close all
            
            N  = 2^7;
            dt = 0.1;
            t  = [0:N-1]' * dt;
            
            nyquistOmegaRps = pi/dt;
            omega1 = nyquistOmegaRps*0.25;   % Survives   tfHighFreqLimitFraction.
            omega2 = nyquistOmegaRps*0.50;   % Removed by tfHighFreqLimitFraction.
            y1_1   = sin(omega1*t);
            y1_2   = sin(omega2*t);
            y1     = y1_1 + y1_2;
                        
            tf     = bicas.tf.utest_utils.get_tf_delay(1*dt);
            
            if 1

                [y2, y1B, y2B, tfB] = bicas.tf.apply_TF(...
                    dt, y1, tf, ...
                    'method',                  'FFT', ...
                    'detrendingDegreeOf',      -1, ...
                    'retrendingEnabled',       false, ...
                    'tfHighFreqLimitFraction', 0.4);

                y2_exp = circshift(y1_1, 1);   % Requires FFT method.
                %bicas.tf.apply_TF___UTEST.plot_test(y1, y2, y2_exp)

                testCase.verifyEqual(abs(tfB(omega1)), 1)
                testCase.verifyEqual(abs(tfB(omega2)), 0)
                testCase.verifyEqual(y2, y2_exp, 'AbsTol', 1e-13)
            end
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
