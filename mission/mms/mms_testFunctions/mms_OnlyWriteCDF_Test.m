% MMS_ONLYWRITECDF_TEST is a unit testing framework for testing writing MMS
% CDF functions only.
%       results = MMS_ONLYWRITECDF_TEST creates a unit testing framework
%       for testing writing to MMS specific CDF files. And make sure that
%       they are readable as expected afterwards. 
%       These functions require Matlab R2013b or later.
%
%       Example:
%               results = MMS_ONLYWRITECDF_TEST
%               results.run
%
%       See also MATLAB.UNITTEST.


function tests = mms_OnlyWriteCDF_Test
    ve=version;
    if(str2double(ve(1))<8)
        error('Require at least R2013b to run this test. Please upgrade.');
    elseif(str2double(ve(3))<2 && (str2double(ve(1))==8))
        error('Require at least R2013b to run this test. Please upgrade.');
    end
    global ENVIR MMS_CONST;
    if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
    ENVIR = mms_sdc_sdp_init; % dummy scNumber used for this test..
    tests = functiontests(localfunctions);
end


function testUscWriteCDF(testCase)
    % Test to write one Usc CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.
    global ENVIR;
    % Fill data with at least 10 points as we have experienced issues with
    % 3x3 data in output file.
    epochTmp = [int64(481907669295761138); int64(481907669327011502); ...
        int64(481907669358261866); int64(481907669389512230); ...
        int64(481907669420762594); int64(481907669452012958); ...
        int64(481907669483263323); int64(481907669514513687); ...
        int64(481907669545764051); int64(481907669577014415)];
    % Fill data with easily identifiable numerical values.
    data1Tmp = [single(0.1); single(2); single(3); single(4); ...
        single(5); single(6); single(7); single(8); single(9); single(10)];
    data2Tmp = [single(1); single(0.12); single(13); single(14); ...
        single(15); single(16); single(17); single(18); single(19); ...
        single(20)];
    data3Tmp = [single(1); single(22); single(0.23); single(24); ...
        single(25); single(26); single(27); single(28); single(29); ...
        single(30)];
    
    psp_p=[data1Tmp, data2Tmp, data3Tmp, data1Tmp, data2Tmp, data3Tmp];
    bitmask= [uint16(1); uint16(2); uint16(3); uint16(4); uint16(5); ...
        uint16(6); uint16(7); uint16(8); uint16(9); uint16(10)];
    oldDir=pwd;
    cd(ENVIR.DROPBOX_ROOT);
    
    mms_sdc_sdp_cdfwrite( ...
        'mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf', int8(2), ...
        'usc', epochTmp, data1Tmp, data2Tmp, data3Tmp, psp_p, bitmask);
    
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj( [ENVIR.DROPBOX_ROOT, ...
        '/mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf']);
    
    % Do some checks that the written output and subsequent reading was
    % as expected.
    actSolution = dataObjIn.data.mms2_sdp_escp_dcv.nrec;
    expSolution = 10;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Verify exactly same output as input.
    
    % EstimatedSpaceCraftPotential, ESCP
    actSolution = dataObjIn.data.mms2_sdp_escp_dcv.data;
    expSolution = data1Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % ProbetoSpacecraftPotential, PSP
    actSolution = dataObjIn.data.mms2_sdp_psp_dcv.data;
    expSolution = data2Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Delta 
    actSolution = dataObjIn.data.mms2_sdp_delta_dcv.data;
    expSolution = data3Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Probe to Spacecraft Potential individual Probe, PSP_P
    actSolution = dataObjIn.data.mms2_sdp_psp_probes_dcv.data;
    expSolution = psp_p;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Bitmask
    actSolution = dataObjIn.data.mms2_sdp_bitmask_dcv.data;
    expSolution = bitmask;
    verifyEqual(testCase,actSolution,expSolution);
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf
    cd(oldDir);
end


function testQuickLookWriteCDF(testCase)
    % Test to write one SITL CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.
    global ENVIR;
    % Fill data with at least 10 points as we have experienced issues with
    % 3x3 data in output file.
    epochTmp = [int64(481907669295761138); int64(481907669327011502); ...
        int64(481907669358261866); int64(481907669389512230); ...
        int64(481907669420762594); int64(481907669452012958); ...
        int64(481907669483263323); int64(481907669514513687); ...
        int64(481907669545764051); int64(481907669577014415)];
    % Fill data with easily identifiable numerical values.
    data1Tmp = [single(0.1); single(2); single(3); single(4); ...
        single(5); single(6); single(7); single(8); single(9); single(10)];
    data2Tmp = [single(1); single(0.12); single(13); single(14); ...
        single(15); single(16); single(17); single(18); single(19); ...
        single(20)];
    data3Tmp = [single(1); single(22); single(0.23); single(24); ...
        single(25); single(26); single(27); single(28); single(29); ...
        single(30)];
    
    data4Tmp = [data1Tmp, data2Tmp, data3Tmp];
    data5Tmp = [data3Tmp, data2Tmp, data1Tmp];
    
    bitmask= [uint16(1); uint16(2); uint16(3); uint16(4); uint16(5); ...
        uint16(6); uint16(7); uint16(8); uint16(9); uint16(10)];
    qualityMark= [uint16(10); uint16(9); uint16(8); uint16(7); ...
        uint16(6); uint16(5); uint16(4); uint16(3); uint16(2); uint16(10)];
    
    oldDir = pwd;
    cd(ENVIR.DROPBOX_ROOT);
    
    mms_sdc_sdp_cdfwrite( ...
        'mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf', int8(2), ...
        'ql', epochTmp,  data4Tmp, data5Tmp, bitmask, qualityMark);
    
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj( [ENVIR.DROPBOX_ROOT, ...
        '/mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf']);
    
    % Do some checks that the written output and subsequent reading was
    % as expected.
    actSolution = dataObjIn.data.mms2_sdp_epoch_dce.nrec;
    expSolution = 10;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Verify exactly same output as input.
    
    % Epoch times in TT2000
    actSolution = dataObjIn.data.mms2_sdp_epoch_dce.data;
    expSolution = epochTmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % DCE xyz in PGSE ref. frame
    actSolution = dataObjIn.data.mms2_sdp_dce_xyz_pgse.data;
    expSolution = data4Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % DCE xyz in DSL ref. frame
    actSolution = dataObjIn.data.mms2_sdp_dce_xyz_dsl.data;
    expSolution = data5Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Bitmask
    actSolution = dataObjIn.data.mms2_sdp_dce_bitmask.data;
    expSolution = bitmask;
    verifyEqual(testCase,actSolution,expSolution);
    
    % QualityMark
    actSolution = dataObjIn.data.mms2_sdp_dce_quality.data;
    expSolution = qualityMark;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf
    cd(oldDir);
end


function testSITLwriteCDF(testCase)
     % Test to write one SITL CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.
    global ENVIR;
    % Fill data with at least 10 points as we have experienced issues with
    % 3x3 data in output file.
    epochTmp = [int64(481907669295761138); int64(481907669327011502); ...
        int64(481907669358261866); int64(481907669389512230); ...
        int64(481907669420762594); int64(481907669452012958); ...
        int64(481907669483263323); int64(481907669514513687); ...
        int64(481907669545764051); int64(481907669577014415)];
    % Fill data with easily identifiable numerical values.
    data1Tmp = [single(0.1); single(2); single(3); single(4); ...
        single(5); single(6); single(7); single(8); single(9); single(10)];
    data2Tmp = [single(1); single(0.12); single(13); single(14); ...
        single(15); single(16); single(17); single(18); single(19); ...
        single(20)];
    data3Tmp = [single(1); single(22); single(0.23); single(24); ...
        single(25); single(26); single(27); single(28); single(29); ...
        single(30)];
    
    data4Tmp = [data1Tmp, data2Tmp, data3Tmp];
    data5Tmp = [data3Tmp, data2Tmp, data1Tmp];
    
    bitmask= [uint16(1); uint16(2); uint16(3); uint16(4); uint16(5); ...
        uint16(6); uint16(7); uint16(8); uint16(9); uint16(10)];
    
    oldDir = pwd;
    cd(ENVIR.DROPBOX_ROOT);
    
    mms_sdc_sdp_cdfwrite( ...
        'mms2_sdp_sitl_l1b_dce2d_20150410000000_v0.0.0.cdf', int8(2), ...
        'sitl', epochTmp,  data4Tmp, data5Tmp, bitmask);
    
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj( [ENVIR.DROPBOX_ROOT, ...
        '/mms2_sdp_sitl_l1b_dce2d_20150410000000_v0.0.0.cdf']);
    
    % Do some checks that the written output and subsequent reading was
    % as expected.
    actSolution = dataObjIn.data.mms2_sdp_epoch_dce.nrec;
    expSolution = 10;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Verify exactly same output as input.
    
    % Epoch times in TT2000
    actSolution = dataObjIn.data.mms2_sdp_epoch_dce.data;
    expSolution = epochTmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % DCE xyz in PGSE ref. frame
    actSolution = dataObjIn.data.mms2_sdp_dce_xyz_pgse.data;
    expSolution = data4Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % DCE xyz in DSL ref. frame
    actSolution = dataObjIn.data.mms2_sdp_dce_xyz_dsl.data;
    expSolution = data5Tmp;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Bitmask
    actSolution = dataObjIn.data.mms2_sdp_dce_bitmask.data;
    expSolution = bitmask;
    verifyEqual(testCase,actSolution,expSolution);
    
    % Delete the output file created, or next run will automatically have
    % errors when trying to write to the same file.
    !rm $DROPBOX_ROOT/mms2_sdp_sitl_l1b_dce2d_20150410000000_v0.0.0.cdf
    cd(oldDir);
end


function testPhaseFromSunpulse(testCase)
    % Verify that the phase calculation from sunpulse works as expected.

    % Fill data with at least 10 points.
    epochData = [...
        int64(481907549000000000); int64(481907550000000000); ...
        int64(481907551000000000); int64(481907552000000000); ...
        int64(481907553000000000); int64(481907554000000000); ...
        int64(481907555000000000); int64(481907556000000000); ...
        int64(481907557000000000); int64(481907558000000000)];
    
    % Fill sunpulse data with at least 10 points, 20 seconds between each
    % (corresponding to 3 rpm).
    epochSunpulse = [...
        int64(481907509000000000); int64(481907529000000000); ...
        int64(481907549000000000); int64(481907569000000000); ...
        int64(481907589000000000); int64(481907609000000000); ...
        int64(481907629000000000); int64(481907649000000000); ...
        int64(481907669000000000); int64(481907689000000000)];
    
    phase = mms_sdc_sdp_phase(epochData, epochSunpulse);
    
    % Compare results with expected results.
    actSolution = phase;
    
    expSolution = [...
        162
        180
        198
        216
        234
        252
        270
        288
        306
        324];
    
    verifyEqual(testCase, actSolution, expSolution);

end
