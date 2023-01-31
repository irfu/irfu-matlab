%
% matlab.unittest automatic test code for
% solo.shk.read_RSH_file(), AND
% solo.shk.read_RSH_file_many().
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-09-07
%
classdef read_RSH_file___UTEST < matlab.unittest.TestCase



    properties(Constant)
        EMPTY_D_EXP = struct(...
            'Name',             {cell(0,1)}, ...
            'EngineeringValue', {cell(0,1)}, ...
            'TimeStampAsciiA',  {cell(0,1)});


        FILE_0_TXT = {...
'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
'<ns2:ResponsePart xmlns:ns2="http://edds.egos.esa/model">'
'  <Response>'
'    <ParamResponse>'
'      <ParamSampleList>'
'      </ParamSampleList>'
'    </ParamResponse>'
'  </Response>'
'</ns2:ResponsePart>'
        };


        FILE_1_TXT = {...
'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
'<ns2:ResponsePart xmlns:ns2="http://edds.egos.esa/model">'
'  <Response>'
'    <ParamResponse>'
'      <ParamSampleList>'
'        <ParamSampleListElement>'
'          <Name>NCAT11X0</Name>'
'          <TimeStampAsciiA>2021-09-04T00:00:02.018306</TimeStampAsciiA>'
'          <Unit>none</Unit>'
'          <Description>THR 1A Cumulative OnTime</Description>'
'          <EngineeringValue>11.5746627642289</EngineeringValue>'
'          <RawValue>11.57466276422888</RawValue>'
'        </ParamSampleListElement>'
'      </ParamSampleList>'
'    </ParamResponse>'
'  </Response>'
'</ns2:ResponsePart>'
}


    FILE_2_TXT = {...
'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
'<ns2:ResponsePart xmlns:ns2="http://edds.egos.esa/model">'
'  <Response>'
'    <ParamResponse>'
'      <ParamSampleList>'
'        <ParamSampleListElement>'
'          <Name>NCAT11X0</Name>'
'          <TimeStampAsciiA>2021-09-04T00:00:32.018343</TimeStampAsciiA>'
'          <Unit>none</Unit>'
'          <Description>THR 1A Cumulative OnTime</Description>'
'          <EngineeringValue>11.5746627642289</EngineeringValue>'
'          <RawValue>11.57466276422888</RawValue>'
'        </ParamSampleListElement>'
'        <ParamSampleListElement>'
'          <Name>NRUD2169</Name>'
'          <TimeStampAsciiA>2021-09-04T00:13:05.405659</TimeStampAsciiA>'
'          <Unit>none</Unit>'
'          <Description>RSA_3_3_7 RPW 2_MLI B</Description>'
'          <EngineeringValue>Deployed</EngineeringValue>'
'          <RawValue>0</RawValue>'
'        </ParamSampleListElement>'
'      </ParamSampleList>'
'    </ParamResponse>'
'  </Response>'
'</ns2:ResponsePart>'
}

    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)

            txtFile = solo.shk.read_RSH_file___UTEST.create_text_file(...
                testCase, testCase.FILE_0_TXT);

            DAct = solo.shk.read_RSH_file(txtFile, []);

            testCase.verifyEqual(DAct, testCase.EMPTY_D_EXP)

            DAct = solo.shk.read_RSH_file(txtFile, {'NCFT55P0'});
            testCase.verifyEqual(DAct, testCase.EMPTY_D_EXP)
        end



        function test1(testCase)

            txtFile = solo.shk.read_RSH_file___UTEST.create_text_file(...
                testCase, testCase.FILE_1_TXT);

            % Cf namesCa = {}.
            DAct = solo.shk.read_RSH_file(txtFile, []);
            testCase.verifyEqual(DAct, ...
                struct(...
                    'Name',             {{'NCAT11X0'}}, ...
                    'EngineeringValue', {{'11.5746627642289'}}, ...
                    'TimeStampAsciiA',  {{'2021-09-04T00:00:02.018306'}}))

            DAct = solo.shk.read_RSH_file(txtFile, {'NCAT11X0'});
            testCase.verifyEqual(DAct, ...
                struct(...
                    'Name',             {{'NCAT11X0'}}, ...
                    'EngineeringValue', {{'11.5746627642289'}}, ...
                    'TimeStampAsciiA',  {{'2021-09-04T00:00:02.018306'}}))

            % Cf namesCa = [].
            DAct = solo.shk.read_RSH_file(txtFile, {});
            testCase.verifyEqual(DAct, testCase.EMPTY_D_EXP)

            DAct = solo.shk.read_RSH_file(txtFile, ...
                {'Nonexistent1', 'Nonexistent2'});
            testCase.verifyEqual(DAct, testCase.EMPTY_D_EXP)
        end



        function test2(testCase)

            txtFile = solo.shk.read_RSH_file___UTEST.create_text_file(...
                testCase, testCase.FILE_2_TXT);

            DAct = solo.shk.read_RSH_file(txtFile, []);
            testCase.verifyEqual(DAct, ...
                struct(...
                    'Name',             {{'NCAT11X0'; 'NRUD2169'}}, ...
                    'EngineeringValue', {{'11.5746627642289'; 'Deployed'}}, ...
                    'TimeStampAsciiA',  {{'2021-09-04T00:00:32.018343'; '2021-09-04T00:13:05.405659'}}))

            DAct = solo.shk.read_RSH_file(txtFile, {'NRUD2169'});
            testCase.verifyEqual(DAct, ...
                struct(...
                    'Name',             {{'NRUD2169'}}, ...
                    'EngineeringValue', {{'Deployed'}}, ...
                    'TimeStampAsciiA',  {{'2021-09-04T00:13:05.405659'}}))
        end



        % solo.shk.read_RSH_file_many()
        % Read 0 files.
        function test10(testCase)

            DAct = solo.shk.read_RSH_file_many(...
                {}, []);

            testCase.verifyEqual(DAct, testCase.EMPTY_D_EXP)
        end



        % solo.shk.read_RSH_file_many()
        % Read 1 file.
        function test11(testCase)
            txtFile1 = solo.shk.read_RSH_file___UTEST.create_text_file(...
                testCase, testCase.FILE_1_TXT);

            DAct = solo.shk.read_RSH_file_many(...
                {txtFile1}, []);

            testCase.verifyEqual(DAct, ...
                struct(...
                    'Name',             {{'NCAT11X0'}}, ...
                    'EngineeringValue', {{'11.5746627642289'}}, ...
                    'TimeStampAsciiA',  {{'2021-09-04T00:00:02.018306'}}))
        end



        % solo.shk.read_RSH_file_many()
        % Read multiple files.
        function test12(testCase)
            txtFile1 = solo.shk.read_RSH_file___UTEST.create_text_file(...
                testCase, testCase.FILE_1_TXT);
            txtFile2 = solo.shk.read_RSH_file___UTEST.create_text_file(...
                testCase, testCase.FILE_2_TXT);

            DAct = solo.shk.read_RSH_file_many(...
                {txtFile1, txtFile2}, []);

            testCase.verifyEqual(DAct, ...
                struct(...
                    'Name',             {{'NCAT11X0'; 'NCAT11X0'; 'NRUD2169'}}, ...
                    'EngineeringValue', {{'11.5746627642289'; '11.5746627642289'; 'Deployed'}}, ...
                    'TimeStampAsciiA',  {{'2021-09-04T00:00:02.018306'; '2021-09-04T00:00:32.018343'; '2021-09-04T00:13:05.405659'}}))
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        function txtFile = create_text_file(testCase, rowsCa)

            DirFixture = testCase.applyFixture(...
                matlab.unittest.fixtures.TemporaryFolderFixture(...
                    'WithSuffix', ['.', mfilename()]));

            txtFile = fullfile(DirFixture.Folder, 'ROC_SOLO_HK.xml');

            irf.fs.write_file(...
                txtFile, ...
                uint8(strjoin(rowsCa, newline)))

        end



    end    % methods(Static, Access=private)



end
