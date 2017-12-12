% Load the July 2016 BSACT files.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden.
% First created 2017-12-12
function [TcDcv, TcTf] = load_BSACT_2016_07(bsactRootPath)

%=================================
% Register DCV calibration tables
%=================================
TcDcv = bicas.BSACT_utils.reader_DCV_TF();
add_DCV_dir(-25, 'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_-25_4.4.txt', 'testlogbook_2016-07-23 _20-01-05__VER_FS.txt');
add_DCV_dir(0,   'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_0_4.4.txt',   'testlogbook_2016-07-24 _20-28-36__VER_FM1_0_4.4.txt');
add_DCV_dir(25,  'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_+25C.txt',    'testlogbook_2016-07-23 _12-00-09__VER_FM1.txt');
add_DCV_dir(50,  'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1_+50c.txt',    'testlogbook_2016-07-24 _16-26-18__VER_FM1_+50C.txt');
add_DCV_dir(70,  'SO_BIAS_DC_VOLTAGE_ID%i_Ver_00_FM1.txt',         'testlogbook_2016-07-24 _12-26-06__VER_FM1.txt');

    function add_DCV_dir(mebTempCelsius, cTableFilesPattern, testlogbookFilename)
        tempSubdir = sprintf('TEMP%iC', mebTempCelsius);
        TcDcv.add_test_directory(...
            fullfile(bsactRootPath, tempSubdir, '4_4_DC_VOLTAGE_TEST', cTableFilesPattern), ...
            fullfile(bsactRootPath, tempSubdir, '4_4_DC_VOLTAGE_TEST', testlogbookFilename), ...
            mebTempCelsius ...
        );
    end

%================================
% Register TF calibration tables
%================================
TcTf = bicas.BSACT_utils.reader_DCV_TF();
add_TF_dir(-25, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_-25_4.5.txt', 'testlogbook_2016-07-23 _21-13-50__VER_FS.txt');
add_TF_dir(  0, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_0_4.5.txt',   'testlogbook_2016-07-24 _18-52-41__VER_4_5_FM1_0C.txt');
add_TF_dir( 25, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_+25C.txt',    'testlogbook_2016-07-23 _10-46-54__VER_FM1.txt');
add_TF_dir( 50, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1_+50c.txt',    'testlogbook_2016-07-24 _15-14-42__VER_Fm1.txt');
add_TF_dir( 70, 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FM1.txt',         'testlogbook_2016-07-24 _11-13-26__VER_FM1.txt');

    function add_TF_dir(mebTempCelsius, cTableFilesPattern, testlogbookFilename)
        tempSubdir = sprintf('TEMP%iC', mebTempCelsius);
        TcTf.add_test_directory(...
            fullfile(bsactRootPath, tempSubdir, '4_5_TRANSFER_FUNCTION', cTableFilesPattern), ...
            fullfile(bsactRootPath, tempSubdir, '4_5_TRANSFER_FUNCTION', testlogbookFilename), ...
            mebTempCelsius ...
        );
    end

end
