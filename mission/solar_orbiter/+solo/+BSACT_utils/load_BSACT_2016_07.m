%
% Load the July 2016 BSACT files.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-12-12
%
function [TcDcc, TcDcv, TcTf, TcIc] = load_BSACT_2016_07(bsactRootPath)

    %=================================
    % Register DCC calibration tables
    %=================================
    TcDcc = solo.BSACT_utils.reader_DCC_DCV_TF_IC();
    add_DCC_dir(-25, '4_3_DC_CURRENT_TEST', 'SO_BIAS_DC_CURRENT_ID%03i_Ver_00_FM1_-25_4.3.txt', 'testlogbook_2016-07-23 _19-41-06__FS.txt');
    add_DCC_dir(0,   '4_3_DC_CURRENT_TEST', 'SO_BIAS_DC_CURRENT_ID%03i_Ver_00_FM1_0_4.3.txt',   'testlogbook_2016-07-24 _20-08-56__FM1_0_4_3.txt');
    add_DCC_dir(25,  '4_3_AC_CURRENT_TEST', 'SO_BIAS_DC_CURRENT_ID%03i_Ver_00_FM1_+25C.txt',    'testlogbook_2016-07-23 _13-11-53__FM1.txt');
    add_DCC_dir(50,  '4_3_AC_CURRENT_TEST', 'SO_BIAS_DC_CURRENT_ID%03i_Ver_00_FM1_+50c.txt',    'testlogbook_2016-07-24 _18-08-22__FM1_+50C.txt');
    add_DCC_dir(70,  '4_3_AC_CURRENT_TEST', 'SO_BIAS_DC_CURRENT_ID%03i_Ver_00_FM1.txt',         'testlogbook_2016-07-24 _13-38-12__FM1.txt');

    %=================================
    % Register DCV calibration tables
    %=================================
    TcDcv = solo.BSACT_utils.reader_DCC_DCV_TF_IC();
    add_DCV_dir(-25, 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_-25_4.4.txt', 'testlogbook_2016-07-23 _20-01-05__VER_FS.txt');
    add_DCV_dir(0,   'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_0_4.4.txt',   'testlogbook_2016-07-24 _20-28-36__VER_FM1_0_4.4.txt');
    add_DCV_dir(25,  'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_+25C.txt',    'testlogbook_2016-07-23 _12-00-09__VER_FM1.txt');
    add_DCV_dir(50,  'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_+50c.txt',    'testlogbook_2016-07-24 _16-26-18__VER_FM1_+50C.txt');
    add_DCV_dir(70,  'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1.txt',         'testlogbook_2016-07-24 _12-26-06__VER_FM1.txt');

    %================================
    % Register TF calibration tables
    %================================
    TcTf = solo.BSACT_utils.reader_DCC_DCV_TF_IC();
    add_TF_dir(-25, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_-25_4.5.txt', 'testlogbook_2016-07-23 _21-13-50__VER_FS.txt');
    add_TF_dir(  0, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_0_4.5.txt',   'testlogbook_2016-07-24 _18-52-41__VER_4_5_FM1_0C.txt');
    add_TF_dir( 25, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_+25C.txt',    'testlogbook_2016-07-23 _10-46-54__VER_FM1.txt');
    add_TF_dir( 50, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_+50c.txt',    'testlogbook_2016-07-24 _15-14-42__VER_Fm1.txt');
    add_TF_dir( 70, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1.txt',         'testlogbook_2016-07-24 _11-13-26__VER_FM1.txt');

    %================================
    % Register IC calibration tables
    %================================
    TcIc = solo.BSACT_utils.reader_DCC_DCV_TF_IC();
    add_IC_dir(-25, 'SO_BIAS_INT_CAL_ID%03i_Ver_00_FM1_-25_4.8.txt', 'testlogbook_2016-07-23 _23-01-49__VER_FS.txt');
    add_IC_dir( 25, 'SO_BIAS_INT_CAL_ID%03i_Ver_00_FM1_+25C.txt',    'testlogbook_2016-07-23 _17-10-23__VER_FM1.txt');
    add_IC_dir( 50, 'SO_BIAS_INT_CAL_ID%03i_Ver_00_FM1_+50c.txt',    'testlogbook_2016-07-24 _17-41-02__VER_FM1_+50C.txt');
    add_IC_dir( 70, 'SO_BIAS_INT_CAL_ID%03i_Ver_00_FM1.txt',         'testlogbook_2016-07-24 _13-56-06__VER_FM1.txt');



    % Function to compress the code. (Use bsactRootPath, TcDcc as "global variables".)
    % Needs argument "testDirName" since it can be both 4_3_DC_CURRENT_TEST and 4_3_AC_CURRENT_TEST (a misspelling:
    % DC-->AC).
    function add_DCC_dir(mebTempCelsius, testDirName, cTableFilesPattern, testlogbookFilename)
        subdir = sprintf('TEMP%dC', mebTempCelsius);
        TcDcc.add_test_directory(...
            'DCC', ...
            fullfile(bsactRootPath, subdir, testDirName, cTableFilesPattern), ...
            fullfile(bsactRootPath, subdir, testDirName, testlogbookFilename), ...
            mebTempCelsius ...
            );
    end



    % Function to compress the code. (Use bsactRootPath, TcDcv as "global variables".)
    function add_DCV_dir(mebTempCelsius, cTableFilesPattern, testlogbookFilename)
        subdir = sprintf('TEMP%dC', mebTempCelsius);
        TcDcv.add_test_directory(...
            'DCV', ...
            fullfile(bsactRootPath, subdir, '4_4_DC_VOLTAGE_TEST', cTableFilesPattern), ...
            fullfile(bsactRootPath, subdir, '4_4_DC_VOLTAGE_TEST', testlogbookFilename), ...
            mebTempCelsius ...
            );
    end



    % Function to compress the code. (Use bsactRootPath, TcTf as "global variables".)
    function add_TF_dir(mebTempCelsius, cTableFilesPattern, testlogbookFilename)
        subdir = sprintf('TEMP%dC', mebTempCelsius);
        TcTf.add_test_directory(...
            'TF', ...
            fullfile(bsactRootPath, subdir, '4_5_TRANSFER_FUNCTION', cTableFilesPattern), ...
            fullfile(bsactRootPath, subdir, '4_5_TRANSFER_FUNCTION', testlogbookFilename), ...
            mebTempCelsius ...
            );
    end

    % Function to compress the code. (Use bsactRootPath, TcTf as "global variables".)
    function add_IC_dir(mebTempCelsius, cTableFilesPattern, testlogbookFilename)
        subdir = sprintf('TEMP%dC', mebTempCelsius);
        TcIc.add_test_directory(...
            'IC', ...
            fullfile(bsactRootPath, subdir, '4_8_INTERNAL_CALIBRATION', cTableFilesPattern), ...
            fullfile(bsactRootPath, subdir, '4_8_INTERNAL_CALIBRATION', testlogbookFilename), ...
            mebTempCelsius ...
            );
    end
end
