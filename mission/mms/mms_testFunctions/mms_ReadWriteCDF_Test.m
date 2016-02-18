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
    % Verify Matlab R2013b or later.
    if( verLessThan('matlab', '8.3') )
        error('Require R2013b or later to run this test. Please upgrade.');
    end
    tests = functiontests(localfunctions);
end

function testReadCDF(testCase)
    % Read one of the predefined MMS SDP CDF file. This one is used as
    % source file for processing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    dataObj = dataobj( [DATA_PATH_ROOT, ...
        '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dce_20160101_v2.0.1.cdf']);
    actSolution = dataObj.data.mms2_sdp_dce_sensor.nrec;
    expSolution = 1605920;
    verifyEqual(testCase,actSolution,expSolution);
end

function testSCpotProcessAndReadCDF(testCase)
    % Test to write one Usc CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('scpot', ...
       [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dcv_20160101_v2.0.1.cdf'], ...
       [DATA_PATH_ROOT, '/hk/mms2/fields/2016/01/mms2_fields_hk_l1b_101_20160101_v0.1.0.cdf'], ...
       [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dce_20160101_v2.0.1.cdf'],...
       [DATA_PATH_ROOT, '/MRT9/mms2/mms2_fields_hk_l1b_10e_20160101_v0.1.0.cdf']);

% Or to force an error, simply do not provide all required inputs.
%  mms_sdc_sdp_proc('usc', ...
%        [DATA_PATH_ROOT, '/hk/mms2/fields/2015/04/mms2_fields_hk_l1b_101_20150410_v0.1.0.cdf'],...
%        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v1.0.1.cdf'],'');
%
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_edp_fast_l2_scpot_20160101000000_v1.0.0.cdf']);
    actSolution = dataObjIn.data.mms2_edp_scpot.nrec;
    expSolution = 1605920;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_edp_fast_l2_scpot_20160101000000_v1.0.0.cdf
end


function testSITLprocessAndReadCDF(testCase)
    % Test to write one SITL CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.end
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('sitl', ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dce_20160101_v2.0.1.cdf'], ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dcv_20160101_v2.0.1.cdf'], ...
        [DATA_PATH_ROOT, '/hk/mms2/fields/2016/01/mms2_fields_hk_l1b_101_20160101_v0.1.0.cdf']);
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_edp_fast_sitl_dce2d_20160101000000_v1.0.0.cdf']);
    actSolution = dataObjIn.data.mms2_edp_dce_xyz_dsl.nrec;
    expSolution = 1605920;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_edp_fast_sitl_dce2d_20160101000000_v1.0.0.cdf
end


function testQuickLookProcessAndReadCDF(testCase)
    % Test to write one QuickLook CDFend file to $DROPBOX_ROOT. The output 
    % file is removed afterwards to ensure it does not interfer with future
    % writing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('ql', ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dce_20160101_v2.0.1.cdf'], ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dcv_20160101_v2.0.1.cdf'], ...
        [DATA_PATH_ROOT, '/hk/mms2/fields/2016/01/mms2_fields_hk_l1b_101_20160101_v0.1.0.cdf']);
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_edp_fast_ql_dce2d_20160101000000_v1.0.0.cdf']);
    actSolution = dataObjIn.data.mms2_edp_dce_xyz_dsl.nrec;
    expSolution = 1605920;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_edp_fast_ql_dce2d_20160101000000_v1.0.0.cdf
end

function testL2PRErocessAndReadCDF(testCase)
    % Test to write one QuickLook CDFend file to $DROPBOX_ROOT. The output
    % file is removed afterwards to ensure it does not interfer with future
    % writing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    mms_sdc_sdp_proc('l2pre', ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dce_20160101_v2.0.1.cdf'], ...
        [DATA_PATH_ROOT, '/science/mms2/sdp/fast/l1b/2016/01/01/mms2_sdp_fast_l1b_dcv_20160101_v2.0.1.cdf'], ...
        [DATA_PATH_ROOT, '/hk/mms2/fields/2016/01/mms2_fields_hk_l1b_101_20160101_v0.1.0.cdf']);
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT, ...
        '/mms2_edp_fast_l2_dce2d_20160101000000_v1.0.0.cdf']);
    actSolution = dataObjIn.data.mms2_edp_dce_xyz_dsl.nrec;
    expSolution = 1605920;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_edp_fast_l2_dce2d_20160101000000_v1.0.0.cdf
end
