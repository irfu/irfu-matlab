%
% matlab.unittest automatic test code for bicas.proc.L1L2.cal.rct.findread.
%
% NOTE: Only tests one function.
% NOTE: Tests for "bicas.proc.L1L2.cal.rct.findread.find_RCT_regexp()" has partly been
% written in order to try out functionality for testing code with file
% operations.
% NOTE: Creates temporary directory and files for every test, separately.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-16
%
classdef findread___UTEST < matlab.unittest.TestCase
    % PROPOSAL: Tests for bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log() for BIAS
    %           RCT. Creates BIAS RCT using bicas.tools.create_RCT() as part of
    %           the test.
    %
    % PROPOSAL: Use TestMethodSetup and TestMethodTeardown.
    %   CON: Directory and files can only be shared among tests for find_RCT_regexp().
    %       CON-PROPOSAL: Move find_RCT_regexp() test code to separate file.
    % PROPOSAL: Use TestClassSetup and TestClassTeardown for logger object.
    %   PRO: Can share logger object.
    %   



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        function test_find_RCT_regexp_empty(testCase)
            [tempDir, L] = bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(testCase, {});
            
            testCase.verifyError(...
                @() bicas.proc.L1L2.cal.rct.findread.find_RCT_regexp(...
                tempDir, '20[0-9][0-9]\.cdf', L), ...
                'BICAS:CannotFindRegexMatchingRCT')            
            
        end
        
        
        
        function test_find_RCT_regexp_no_match(testCase)
            [tempDir, L] = bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
                testCase, {'20201.cdf', '2020.CDF'});
            
            testCase.verifyError(...
                @() bicas.proc.L1L2.cal.rct.findread.find_RCT_regexp(...
                tempDir, '20[0-9][0-9]\.cdf', L), ...
                'BICAS:CannotFindRegexMatchingRCT')
        end
        
        
        
        function test_find_RCT_regexp_1_match(testCase)
            [tempDir, L] = bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
                testCase, {'2020.cdf', 'asdsf'});
            
            path = bicas.proc.L1L2.cal.rct.findread.find_RCT_regexp(...
                tempDir, '20[0-9][0-9]\.cdf', L);
            
            testCase.verifyEqual(...
                path, ...
                fullfile(tempDir, '2020.cdf'))            
        end
        
        
        
        function test_find_RCT_regexp_2_match(testCase)
            [tempDir, L] = bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
                testCase, {'2020.cdf', '2021.cdf'});
            
            path = bicas.proc.L1L2.cal.rct.findread.find_RCT_regexp(...
                tempDir, '20[0-9][0-9]\.cdf', L);
            
            testCase.verifyEqual(...
                path, ...
                fullfile(tempDir, '2021.cdf'))            
        end
        
        
        
        function test_find_RCT_regexp_realistic(testCase)
            FN_1 = 'SOLO_CAL_RPW-BIAS_V202111191204.cdf';
            FN_2 = 'SOLO_CAL_RPW-BIAS_V202011191204.cdf';
            
            [tempDir, L] = bicas.proc.L1L2.cal.rct.findread___UTEST.setup_files(...
                testCase, {FN_1, FN_2});
            
            path = bicas.proc.L1L2.cal.rct.findread.find_RCT_regexp(...
                tempDir, 'SOLO_CAL_RPW-BIAS_V20[0-9]{10,10}.cdf', L);
            
            testCase.verifyEqual(...
                path, ...
                fullfile(tempDir, FN_1))
        end
        
        
        
    end    % methods(Test)
        
        
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Create (temporary) directory with specified empty files.
        %
        function [tempDir, L] = setup_files(testCase, filenamesCa)
            L = bicas.Logger('none', false);
            
            DirFixture = testCase.applyFixture(...
                matlab.unittest.fixtures.TemporaryFolderFixture(...
                    'WithSuffix', ['.', mfilename()]));
            
            tempDir = DirFixture.Folder;
            %fprintf('tempDir = %s\n', tempDir)

            for fileCa = filenamesCa(:)'
                % Create empty file.
                irf.fs.write_file(...
                    fullfile(tempDir, fileCa{1}), ...
                    zeros(1,0,'uint8'))
            end

        end
   
        
        
    end    % methods(Static, Access=private)

    
    
end
