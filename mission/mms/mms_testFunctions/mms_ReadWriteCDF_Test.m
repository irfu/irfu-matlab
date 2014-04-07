% MMS_READWRITECDF_TEST is a unit testing framework for testing various MMS
% Matlab script functions.
%       results = MMS_READWRITECDF_TEST creates a unit testing framework of
%       several test. Each designed to test various parts of the MMS
%       processing. These functions require Matlab R2013b or later.
%
%       Example:
%               results = MMS_READWRITECDF_TEST
%               results.run
%
%       See also MATLAB.UNITTEST.


function tests = mms_ReadWriteCDF_Test
    ve=version;
    if( str2double(ve(1))<8 )
        error('Require at least R2013b to run this test. Please upgrade.');
    elseif( (str2double(ve(3))<2) && (str2double(ve(1))==8) )
        error('Require at least R2013b to run this test. Please upgrade.');
    end
    tests = functiontests(localfunctions);
end

function testReadCDF(testCase)
    % Read one of the predefined MMS SDP CDF file. This one is used as
    % source file for processing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    dataObj = dataobj( [DATA_PATH_ROOT, ...
        '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf'],...
        'KeepTT2000');
    actSolution = dataObj.data.mms2_sdp_dce_sensor.nrec;
    expSolution = 445919;
    verifyEqual(testCase,actSolution,expSolution);
end

function testUscProcessAndReadCDF(testCase)
    % Test to write one Usc CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('usc', ...
       [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dcv_20150410_v0.1.3.cdf'], ...
       [DATA_PATH_ROOT, '/tmp_HK/mms2_fields_hk_l1b_101_20150410_v0.0.1.cdf'], ...
       [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf']);

% Or to force an error, simply do not provide all required inputs.
%  mms_sdc_sdp_proc('usc', ...
%        [DATA_PATH_ROOT, '/tmp_HK/mms2_fields_hk_l1b_101_20150410_v0.0.1.cdf'],...
%        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf'],'');
%    
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf'], 'KeepTT2000');
    actSolution = dataObjIn.data.mms2_sdp_escp_dcv.nrec;
    expSolution = 445919;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf
end


function testSITLprocessAndReadCDF(testCase)
    % Test to write one SITL CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.end
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('sitl', ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf'], ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dcv_20150410_v0.1.3.cdf'], ...
        [DATA_PATH_ROOT, '/tmp_HK/mms2_fields_hk_l1b_101_20150410_v0.0.1.cdf']);
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_sdp_fast_sitl_dce2d_20150410000000_v0.0.0.cdf'], 'KeepTT2000');
    actSolution = dataObjIn.data.mms2_sdp_dce_xyz_dsl.nrec;
    expSolution = 445919;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT//mms2_sdp_fast_sitl_dce2d_20150410000000_v0.0.0.cdf
end


function testQuickLookProcessAndReadCDF(testCase)
    % Test to write one QuickLook CDFend file to $DROPBOX_ROOT. The output 
    % file is removed afterwards to ensure it does not interfer with future
    % writing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('ql', ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf'], ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dcv_20150410_v0.1.3.cdf'], ...
        [DATA_PATH_ROOT, '/tmp_HK/mms2_fields_hk_l1b_101_20150410_v0.0.1.cdf']);
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf'], 'KeepTT2000');
    actSolution = dataObjIn.data.mms2_sdp_dce_xyz_dsl.nrec;
    expSolution = 445919;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf
end