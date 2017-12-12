% Load the June 2016 BSACT files.
%
% NOTE: The MEB temperature used for this calibration test is not yet known (2017-12-12).
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden.
% First created 2017-12-12
function [TcDcv, TcTf] = load_BSACT_2016_06(bsactRootPath)

MEB_TEMPERATURE_CELSIUS = 22;  % Guessing ~room temperature. /Erik P G Johansson 2017-12-12.

%=================================
% Register DCV calibration tables
%=================================
TcDcv = bicas.BSACT_utils.reader_DCV_TF();
TcDcv.add_test_directory(...
    fullfile(bsactRootPath, '4-4_BIAS_DC_VOLTAGE/', 'SO_BIAS_DC_VOLTAGE_ID%03i_Ver_00_FS0_PAFM.txt'), ...
    fullfile(bsactRootPath, '4-4_BIAS_DC_VOLTAGE/', 'Log', 'testlogbook_2016-06-21 _16-06-21__VER_FS.txt'), ...
    MEB_TEMPERATURE_CELSIUS);

%================================
% Register TF calibration tables
%================================
TcTf  = bicas.BSACT_utils.reader_DCV_TF();
TcTf.add_test_directory(...
    fullfile(bsactRootPath, '4-5_TRANSFER_FUNCTION', 'SO_BIAS_AC_VOLTAGE_ID%02i_Ver_00_FS0_PAFM.txt'), ...
    fullfile(bsactRootPath, '4-5_TRANSFER_FUNCTION', 'Log', 'testlogbook_2016-06-22 _09-19-38__VER_FS.txt'), ...
    MEB_TEMPERATURE_CELSIUS);

end
