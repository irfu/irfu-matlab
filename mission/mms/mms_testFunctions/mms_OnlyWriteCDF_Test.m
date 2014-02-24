% MMS_ONLYWRITECDF_TEST is a unit testing framework for testing writing MMS
% CDF functions only.
%       results = MMS_ONLYWRITECDF_TEST creates a unit testing framework
%       for testing writing to MMS specific CDF files. And make sure that
%       they are readable as expected afterwards..
%
%       Example:
%               results = MMS_ONLYWRITECDF_TEST
%               results.run
%
%       See also MATLAB.UNITTEST.


function tests = mms_ReadWriteCDF_Test
    tests = functiontests(localfunctions);
end


function testUscProcessAndReadCDF(testCase)
    % Test to write one Usc CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
    % writing.
    DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
    DROPBOX_ROOT = getenv('DROPBOX_ROOT');
    % Fill data with at least 10 points as we have experienced issues with
    % 3x3 data in output file.
    epochTmp = [int64(481907669295761138); int64(481907669327011502); int64(481907669358261866); int64(481907669389512230); int64(481907669420762594); int64(481907669452012958); int64(481907669483263323); int64(481907669514513687); int64(481907669545764051); int64(481907669577014415)];
    % Fill data with easily identifiable numerical values.
    data1Tmp = [single(0.1); single(2); single(3); single(4); single(5); single(6); single(7); single(8); single(9); single(10)];
    data2Tmp = [single(1); single(0.12); single(13); single(14); single(15); single(16); single(17); single(18); single(19); single(20)];
    data3Tmp = [single(1); single(22); single(0.23); single(24); single(25); single(26); single(27); single(28); single(29); single(30)];
    
    psp_p=[data1Tmp, data2Tmp, data3Tmp, data1Tmp, data2Tmp, data3Tmp];
    bitmask= [uint16(1); uint16(2); uint16(3); uint16(4); uint16(5); uint16(6); uint16(7); uint16(8); uint16(9); uint16(10)];
    cd(DROPBOX_ROOT);
    try
            mms_sdc_cdfwrite('mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf', int8(2), 'usc', epochTmp, data1Tmp, data2Tmp, data3Tmp, psp_p', bitmask);
        catch err
            % An error occured.
            % Give more information for mismatch.
            if (strcmp(err.identifier,'MATLAB:mms_sdc_cdfwrite:filename_output:exists'))
                % If our cdfwrite code resulted in error write proper log message.
                error('MATLAB:SDCcode', '183.. Matlab:MMS_sdc_cdfwrite:filename_output:exist a file with the requested filename already exist in the output dir. Not going to overwrite it!');
            else
                % Display any other errors as usual.
                rethrow(err);
            end
    end % End of try
    
    % If no error was return for full processing try reading the output
    % file created and verify number of record is correct.
    dataObjIn = dataobj([DROPBOX_ROOT,'/mms2_sdp_fast_l2_uscdcv_20150410000000_v0.0.0.cdf'],'tint',0,'true');
    
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
end


function testSITLprocessAndReadCDF(testCase)
    % Test to write one SITL CDF file to $DROPBOX_ROOT. The output file is 
    % removed afterwards to ensure it does not interfer with future
%     % writing.end
%     DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
%     DROPBOX_ROOT = getenv('DROPBOX_ROOT');
%     mms_sitl_dce([DATA_PATH_ROOT,'/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf'],[DATA_PATH_ROOT,'/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dcv_20150410_v0.1.3.cdf']);
%     % If no error was return for full processing try reading the output
%     % file created and verify number of record is correct.
%     dataObjIn = mms_cdf_in_process([DROPBOX_ROOT,'/mms2_sdp_sitl_l1b_dce2d_20150410000000_v0.0.0.cdf'],'sci');
%     actSolution = dataObjIn.data.mms2_sdp_dce_xyz_dsl.nrec;
%     expSolution = 445919;
%     verifyEqual(testCase,actSolution,expSolution);
%     % Delete the output file created, or next run will automatically have
%     % errors when trying to write to the same file.
%     !rm $DROPBOX_ROOT/mms2_sdp_sitl_l1b_dce2d_20150410000000_v0.0.0.cdf
end


function testQuickLookProcessAndReadCDF(testCase)
    % Test to write one QuickLook CDFend file to $DROPBOX_ROOT. The output 
    % file is removed afterwards to ensure it does not interfer with future
    % writing.
%     DATA_PATH_ROOT = getenv('DATA_PATH_ROOT');
%     DROPBOX_ROOT = getenv('DROPBOX_ROOT');
%     mms_ql_dce([DATA_PATH_ROOT,'/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dce_20150410_v0.1.3.cdf'], [DATA_PATH_ROOT,'/science/mms2/sdp/fast/l1b/2015/04/10/mms2_sdp_fast_l1b_dcv_20150410_v0.1.3.cdf']);
%     % If no error was return for full processing try reading the output
%     % file created and verify number of record is correct.
%     dataObjIn = mms_cdf_in_process([DROPBOX_ROOT,'/mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf'],'sci');
%     actSolution = dataObjIn.data.mms2_sdp_dce_xyz_dsl.nrec;
%     expSolution = 445919;
%     verifyEqual(testCase,actSolution,expSolution);
%     % Delete the output file created, or next run will automatically have
%     % errors when trying to write to the same file.
%     !rm $DROPBOX_ROOT/mms2_sdp_fast_ql_dce2d_20150410000000_v0.0.0.cdf
end