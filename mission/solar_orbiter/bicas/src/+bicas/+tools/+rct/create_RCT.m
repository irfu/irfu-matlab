%
% Script/utility for creating the ROC-SGSE and RODP BIAS RCTs. The code reads a
% master file, and using it, creates a version with calibration data added to
% it.
%
% NOTE: The actual calibration data is hard-coded in this file.
% NOTE: Will overwrite old RCT file.
% NOTE: Generates using the default CDF format version. Might not be the one
%       that is expected for archiving purposes.
%
%
% ARGUMENTS
% =========
% rctMasterCdfFile
%       Path to empty BIAS master RCT CDF file (existing) to read.
% destDir
%       Path to directory of the RCT file to be created (compliant filename will
%       be generated by the function). NOTE: Any pre-existing destination file
%       will be overwritten.
%
%
% RETURN VALUE
% ============
% rctPath : Path to the file created. Useful for printing log messages.
%
%
% RATIONALE
% =========
% An alternative to using this function is to modify the master RCT manually
% with cdfedit (from NASA SPDF's CDF software). Using this function is however
% better because
% (1) There are two RCTs with identical calibration data content (ROC-SGSE +
%     RODP),
% (2) it saves work when updating the master RCT .cdf (generated from master
%     Excel file), including making unofficial intermediate versions in the
%     master (e.g. fixing typos)
% (3) it is easier to edit (compared to using cdfedit), e.g. (a) editing
%     existing information as well as (b) adding records
% (4) it is easier to edit if the RCT is extended with more data by changing the
%     format (in particular more transfer functions, by some multiple).
% (5) it could be modified to read data from external files, e.g. text files.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-03-09
%
function rctPath = create_RCT(rctMasterCdfFile, destDir, beginDt, endDt, versionNbr)
% PROPOSAL: Change function name: Something which implies using a master file and "filling it".
% PROPOSAL: Somehow separate the code with the hard-coded data into a separate file.
%
% PROPOSAL: Merge definition of bicas.ga.mods.Database() with that of BICAS
%           proper (bicas.const.GA_MODS_DB).
%   CON: Must then use BICAS versions.
%   CON: Obscures the definition, since RCTs are a bit separate from BICAS
%        proper. BICAS must be compatible with the RCT format, but the method of
%        deriving the calibration values (e.g. curve fitting) can be independent
%        of BICAS.

GMDB = bicas.ga.mods.Database({bicas.const.RCT_DSI});
% NOTE: Using BICAS version in MODS. Not obvious that one should.
GMDB.add_GMVE({bicas.const.RCT_DSI}, ...
  bicas.ga.mods.VersionEntry('2024-07-12', '8.1.0', ...
  {'Updated with compliant metadata and filename.'}))
GA_MODS = GMDB.get_MODS_strings_CA(bicas.const.RCT_DSI);



rctFilename = bicas.tools.rct.create_RCT_filename(beginDt, endDt, versionNbr);
rctPath     = fullfile(destDir, rctFilename);

[RctZvL, RctZvH] = set_RCT_content();
create_RCT_file(rctMasterCdfFile, rctPath, RctZvL, RctZvH, versionNbr, GA_MODS);
end



function [RctZvL, RctZvH] = set_RCT_content()

ADD_DEBUG_RECORD_L = 0;
ADD_DEBUG_RECORD_H = 0;
C.N_ZV_COEFF       = 8;

%===================================================================
% Create EMPTY (not zero-valued) variables representing zVariables.
%===================================================================
RctZvL.Epoch_L                  = int64( zeros(0,1) );
RctZvL.BIAS_CURRENT_OFFSET      = zeros(0, 3);
RctZvL.BIAS_CURRENT_GAIN        = zeros(0, 3);
RctZvL.TRANSFER_FUNCTION_COEFFS = zeros(0, 2, C.N_ZV_COEFF, 4);

RctZvH.Epoch_H  = int64( zeros(0,1) );
RctZvH.E_OFFSET = zeros(0,3);
RctZvH.V_OFFSET = zeros(0,3);



%===================================================================================================================
% --------------------
% DC single
%            -1.041e10 s^2 + 8.148e14 s - 5.009e20
%   ---------------------------------------------------------
%   s^4 + 8.238e05 s^3 + 2.042e12 s^2 + 2.578e17 s + 8.556e21
%
% DC diff
%             2.664e11 s^2 - 1.009e18 s - 2.311e23
%   ---------------------------------------------------------
%   s^4 + 3.868e06 s^3 + 7.344e12 s^2 + 4.411e18 s + 2.329e23
%
% AC, diff, low-gain (gamma=5)
%            -1.946e12 s^2 - 1.365e18 s - 2.287e18
%   -------------------------------------------------------
%   s^4 + 3.85e06 s^3 + 4.828e12 s^2 + 2.68e17 s + 1.348e19
%
% AC, diff, gamma=100
%             1.611e24 s^4 - 2.524e30 s^3 - 1.258e35 s^2 - 4.705e39 s + 2.149e40
%   ---------------------------------------------------------------------------------------
%   s^6 + 7.211e17 s^5 + 6.418e23 s^4 + 6.497e28 s^3 + 2.755e33 s^2 + 4.817e37 s + 2.114e39
%
% Based on BIAS standalone calibrations 2016-06-21/22, 100 kOhm stimuli,
% (there is only one temperature for these tests), TEST ID=0-3 Fits have
% been made using MATLAB function invfreqs with weights = 1 for freqHz <=
% 199e3.
%-------------------------------------------------------------------------------------------------------------------
% NOTE: Above fits for DC single/diff, AC low gain (NOT AC high gain) can be
% re-created using
%   * Files 20160621_FS0_EG_Test_Flight_Harness_Preamps/4-5_TRANSFER_FUNCTION/SO_BIAS_AC_VOLTAGE_ID{00..02}*.txt
%   * N_ZEROS = 2;
%     N_POLES = 4;
%   * N_ITERATIONS = 30;
%   * weights = double( (Data.freqHz <= 199e3) );
%   * [b, a] = invfreqs(Data.z, Data.freqRps, N_ZEROS, N_POLES, weights, N_ITERATIONS);
% NOTE: Unclear how to re-create the fit for AC high gain, but it should be
% similar but using file
%   20160621_FS0_EG_Test_Flight_Harness_Preamps/4-5_TRANSFER_FUNCTION/SO_BIAS_AC_VOLTAGE_ID03*.txt
% NOTE: All ABOVE TFs (from quoted e-mail) except AC diff high-gain, invert
% the sign at 0 Hz. This sign change appears to be wrong. The source files
% (four) all have ~sign inversion at 10 Hz (the lowest tabulated frequency);
% -141-142 degrees phase for both AC TFs.
% --
% NOTE: Above transfer functions were e-mailed Erik P G Johansson to Thomas
%       Chust 2017-12-13 on request.
% NOTE: Above transfer functions were used in officially delivered RCTs
%       below, but with the OPPOSITE SIGN (numerator):
%   SOLO_CAL_RPW_BIAS_V202004062127.cdf
%   SOLO_CAL_RPW_BIAS_V202003101607.cdf
%===================================================================================================================

%===================================================================================================================
% SC = Sign Changed
TF_DC_single_2017_12_13_SC = {-[-5.009e20,  8.148e14, -1.041e10],                      [8.556e21, 2.578e17, 2.042e12, 8.238e05, 1]};
TF_DC_diff_2017_12_13_SC   = {-[-2.311e23, -1.009e18,  2.664e11],                      [2.329e23, 4.411e18, 7.344e12, 3.868e06, 1]};
TF_AC_LG_2017_12_13_SC     = {-[-2.287e18, -1.365e18, -1.946e12],                      [1.348e19, 2.68e17 , 4.828e12, 3.85e06,  1]};
TF_AC_HG_2017_12_13_SC     = {-[ 2.149e40, -4.705e39, -1.258e35, -2.524e30, 1.611e24], [2.114e39, 4.817e37, 2.755e33, 6.497e28, 6.418e23,  7.211e17, 1]};
%===================================================================================================================
if 0
  % Version used up until 2020-11-19, ~internally.
  RctZvL = add_RCT_ZVs_L(RctZvL, int64(0), ...
    [2.60316e-09, -4.74234e-08, -4.78828e-08]', [1.98008e-09, 1.97997e-09, 1.98021e-09]', ...
    create_tfc_ZV_record(C, ...
    'DC_single', TF_DC_single_2017_12_13_SC, ...
    'DC_diff',   TF_DC_diff_2017_12_13_SC, ...
    'AC_lg',     TF_AC_LG_2017_12_13_SC, ...
    'AC_hg',     TF_AC_HG_2017_12_13_SC));
end
%===================================================================================================================
if 1
  % Version used from 2020-11-19, ~internally.
  % SOLO_CAL_RPW-BIAS_V202011191147.cdf
  TF_AC_LG_2020_11_18 = {...
    fliplr([2.66744e+19, -1.52549e+27, 1.8496e+34, 1.33551e+40, 0]), ...
    fliplr([1, 5.94792e+14, 6.98964e+21, 3.59175e+28, 4.71518e+34, 2.62533e+39, 1.27572e+41])};
  TF_AC_HG_2020_11_18 = {...
    fliplr([3.89496e+32, 2.13364e+37, 6.91633e+41, 0]), ...
    fliplr([1, 1.71685e+20, 9.97355e+25, 1.05546e+31, 4.33219e+35, 7.05355e+39, 3.4202e+41])};

  RctZvL = add_RCT_ZVs_L(RctZvL, int64(0), ...
    [2.60316e-09, -4.74234e-08, -4.78828e-08]', [1.98008e-09, 1.97997e-09, 1.98021e-09]', ...
    create_tfc_ZV_record(C, ...
    'DC_single', TF_DC_single_2017_12_13_SC, ...
    'DC_diff',   TF_DC_diff_2017_12_13_SC, ...
    'AC_lg',     TF_AC_LG_2020_11_18, ...
    'AC_hg',     TF_AC_HG_2020_11_18));
end
%===================================================================================================================



if ADD_DEBUG_RECORD_L
  % TEST: Add another record for Epoch_L.
  RctZvL = add_RCT_ZVs_L(RctZvL, ...
    RctZvL.Epoch_L(end) + 1e9, ...
    RctZvL.BIAS_CURRENT_OFFSET(end)', ...
    RctZvL.BIAS_CURRENT_GAIN(end)', ...
    RctZvL.TRANSFER_FUNCTION_COEFFS(end, :,:,:) ...
    );
  warning('Creating RCT with added test data.')
end



%========================================================================================
% Values from 20160621_FS0_EG_Test_Flight_Harness_Preamps.
% V_OFFSET values from mheader.reg6 for tests with stimuli=1e5 Ohm.
% E_OFFSET values from mheader.reg6 for tests with stimuli=1e5 Ohm, non-inverted inputs.
% Sign is uncertain
% Uncertain whether it is correct to use value for stimuli=1e5 Ohm instead 1e6 Ohm.
% Uncertain whether it is correct to use the reg6 value instead of own fit.
%========================================================================================
RctZvH = add_RCT_ZVs_H(RctZvH, int64(0), -[0.001307, 0.0016914, 0.0030156]', -[0.015384, 0.01582, 0.017215]');
if ADD_DEBUG_RECORD_H
  RECORD_DELAY_NS = 2e9;
  % TEST: Add another record for Epoch_H.
  RctZvH = add_RCT_ZVs_H(RctZvH, ...
    RctZvH.Epoch_H(end) + RECORD_DELAY_NS, ...
    RctZvH.V_OFFSET(end, :)', ...
    RctZvH.V_OFFSET(end, :)');
  warning('Creating RCT with added test data.')
end

end    % function set_RCT_content()



% Create array corresponding to one CDF record of zVar TRANSFER_FUNCTION_COEFFS.
%
% TFC = (zVar) TRANSFER_FUNCTION_COEFFS
%
% ARGUMENTS
% =========
% varargin :
%   (For n = 1,2,3,4:)
%   varargin{2n-1} : String constant representing channel: "DC_single" etc.
%   varargin{2n}   : Length 2 cell array.
%                   {1} = Numerator coefficients
%                   {2} = Denominator coefficients
%                   NOTE: Unusual "syntax" for argument list. String constants
%                   must be the same every time. Only there for safety.
%
function zvRecord = create_tfc_ZV_record(C, varargin)
% PROPOSAL: Convert varargin to struct directly.

zvRecord = double(zeros(1, 2, C.N_ZV_COEFF, 4));
INDEX_LABEL_LIST = {'DC_single', 'DC_diff', 'AC_lg', 'AC_hg'};
for i=1:4
  if strcmp(varargin{2*i-1}, INDEX_LABEL_LIST{i})
    zvRecord(1,:,:,i) = ca2na(varargin{2*i});
  else
    error('Illegal argument string constant.')
  end
end

% ASSERTIONS
assert(all(isfinite(zvRecord), 'all'), ...
  'create_RCT:Assertion', 'zvRecord contains non-finite values.')
irf.assert.sizes(zvRecord, [1, 2, C.N_ZV_COEFF, 4])

%###########################################################################

% Convert TF format
% Convert 1x2 cell array of 1D arrays --> "2D" array, size 1x2xN
  function na = ca2na(ca)
    % ASSERTIONS
    assert(iscell(ca))
    assert(numel(ca) == 2)
    irf.assert.sizes(...
      ca{1}, [1, NaN], ...
      ca{2}, [1, NaN])
    assert(ca{2}(end) == 1, ...
      ['Submitted TF does not have a coefficient 1 for the highest', ...
      ' order denominator term. This assertion should be satisified', ...
      ' if the TF comes from a MATLAB fit, but is otherwise not required.'])

    % Convert cell array to numeric array.
    na = [...
      padarray(ca{1}, [0, C.N_ZV_COEFF-numel(ca{1})], 'post'); ...
      padarray(ca{2}, [0, C.N_ZV_COEFF-numel(ca{2})], 'post') ...
      ];    % 2 x N_ZV_COEFF
    na = permute(na, [3,1,2]);    % 1 x 2 x N_ZV_COEFF

    % ASSERTIONS
    irf.assert.sizes(na, [1, 2, C.N_ZV_COEFF])
  end
end



% Add 1 record to zVariables associated with Epoch_L.
function RctZvL = add_RCT_ZVs_L(...
  RctZvL, ...
  Epoch_L, ...
  BIAS_CURRENT_OFFSET, ...
  BIAS_CURRENT_GAIN, ...
  TRANSFER_FUNCTION_COEFFS)

RctZvL.Epoch_L                 (end+1, 1)       = Epoch_L;
RctZvL.BIAS_CURRENT_OFFSET     (end+1, :)       = BIAS_CURRENT_OFFSET;
RctZvL.BIAS_CURRENT_GAIN       (end+1, :)       = BIAS_CURRENT_GAIN;
RctZvL.TRANSFER_FUNCTION_COEFFS(end+1, :, :, :) = TRANSFER_FUNCTION_COEFFS;
end



% Add 1 record to zVariables associated with Epoch_H.
function RctZvH = add_RCT_ZVs_H(RctZvH, Epoch_H, V_OFFSET, E_OFFSET)
RctZvH.Epoch_H (end+1, 1) = Epoch_H;
RctZvH.V_OFFSET(end+1, :) = V_OFFSET;
RctZvH.E_OFFSET(end+1, :) = E_OFFSET;
end



% RCT GLOBAL ATTRIBUTES
% =====================
% RCTs should contain multiple GAs. Most can be set in skeletons, but below GAs
% must be set by code:
%
% CALIBRATION_TABLE
%   "Filename of the calibration table, without the “.cdf” extension."
%   "Same definition than for RCS output CDF metadata"
% CALIBRATION_VERSION
%   "Version of the calibration table."
%   "Same definition than for RCS output CDF metadata"
% /Quotes from ROC-PRO-PIP-ICD-00037-LES, RCS ICD, version 1/7, Section 4.2.6
%
%
% Xavier Bonnin also implies that "Datetime" should be set as in datasets:
% """"""""
% As a reminder the link between gAttr and filename’s fields is:
%
% <Source_name>_<LEVEL>_<Descriptor>_<Datetime>_V<Data_version>_<Free_field>.cdf
%
% Where,
%
% <Source_name> 	—> Lower case prefix of the Source_name gAttr (here is always "solo")
% <LEVEL>			—> Upper case prefix of the LEVEL gAttr (here is "CAL")
% Descriptor 		—> Lower case prefix of the Descriptor gAttr (here is "rpw-scm")
% Datetime 			—> Datetime gAttr value in upper case (here "20200210-20990101")
% Data_version		—> Data_version gAttr value
% Free_field		—> Lower case value of Free_field gAttr (Free_field gAttr and field are optional)
% """"""""
% /Xavier Bonnin, e-mail 2024-07-09
%
% Datetime
%   "Datetime field as provided in the CDF file name"
%   "e.g., “20200301” (for a daily file)"
% /Quotes from ROC-PRO-PIP-ICD-00037-LES, RCS ICD, version 1/7, Section 4.1.2
% RCS output CDF metadata setting convention (i.e. when describing datasets, not
% RCTs directly).
%
%
% ARGUMENTS
% =========
% RctL
%     Struct with ZVs associated with Epoch_L.
% RctH
%     Struct with ZVs associated with Epoch_H.
%
%
% NOTES ON RCT
% ============
% """"""""
% Helen has confirmed that only the filename matters on SOAR side.
% So CAL CDF files do not need to be ISTP compliant.
%
% I would just recommend that there enough information inside the CAL CDF to
% identify the file, namely providing the following global attributes:
% Source_name, LEVEL, Descriptor, Data_version, Datetime* and Free_field (if
% any).
% """"""""   /Xavier Bonnin, e-mail, 2024-07-08
% NOTE: This is not the same list as in the RCS ICS 01/07
%
%
function create_RCT_file(rctMasterCdfFile, rctFilePath, RctL, RctH, rctVersionNbr, ga_MODS)
% PROPOSAL: Assertion for matching number of records Epoch_L+data, Epoch_H+data.
%   PROPOSAL: Read from master file which should match.
% TODO-DEC: Require correct MATLAB classes (via write_CDF_dataobj)?
% PROPOSAL: Use BICAS code for writing datasets.
%   CON: BICAS can not handle CALIBRATION_TABLE, CALIBRATION_VERSION it is intended
%        here.

assert(isnumeric(rctVersionNbr))

DT_FORMAT_STR = 'yyyy-MM-dd''T''HH:mm:ss''Z''';



DataObj = dataobj(rctMasterCdfFile);

%=========
% Set GAs
%=========
% NOTE: Overwriting previous value in skeleton (SOLO_CAL_RPW-BIAS_V02.cdf). The
%       value should be empty in the skeleton to be less deceiving.
rctFilename = irf.fs.get_name(rctFilePath);
% R = solo.adm.dsfn.parse_dataset_filename(rctFilename);
assert(strcmp(rctFilename(end-3:end), '.cdf'))
ga_CALIBRATION_TABLE = rctFilename(1:end-4);
DataObj.GlobalAttributes.CALIBRATION_TABLE   = {ga_CALIBRATION_TABLE};

DataObj.GlobalAttributes.CALIBRATION_VERSION = {sprintf('%02i', rctVersionNbr)};

% Uncertain if GAs "Data_version" and "CALIBRATION_VERSION" should be identical in
% RCTs, but it does make some sense. They are different science datasets.
DataObj.GlobalAttributes.Data_version = DataObj.GlobalAttributes.CALIBRATION_VERSION;
DataObj.GlobalAttributes.MODS         = ga_MODS;
% GA Datetime is not required by RCS ICD (ROC-PRO-PIP-ICD-00037-LES, 01/07), but
% Xavier Bonnin requested it in e-mail 2024-07-08.
DataObj.GlobalAttributes.Datetime     = R.timeIntervalStr;

% TIME_MIN, TIME_MAX are not required for CAL, but they are in the skeleton, so
% one can just as well set them in the same way they are set for datasets (or
% remove them from the skeleton)
% beginStr = char(datetime(beginDt, 'Format', DT_FORMAT_STR));
% endStr   = char(datetime(endDt,   'Format', DT_FORMAT_STR));
% DataObj.GlobalAttributes.TIME_MIN    = beginStr;
% DataObj.GlobalAttributes.TIME_MAX    = endStr;



%================
% Set zVariables
%================
DataObj.data.Epoch_L.data = RctL.Epoch_L;
DataObj.data.Epoch_H.data = RctH.Epoch_H;

DataObj.data.BIAS_CURRENT_OFFSET.data      = RctL.BIAS_CURRENT_OFFSET;         % Epoch_L
DataObj.data.BIAS_CURRENT_GAIN.data        = RctL.BIAS_CURRENT_GAIN;           % Epoch_L
DataObj.data.TRANSFER_FUNCTION_COEFFS.data = RctL.TRANSFER_FUNCTION_COEFFS;    % Epoch_L
DataObj.data.E_OFFSET.data                 = RctH.E_OFFSET;                    % Epoch_H
DataObj.data.V_OFFSET.data                 = RctH.V_OFFSET;                    % Epoch_H



%============
% Create CDF
%============
irf.cdf.write_dataobj(...
  rctFilePath, ...
  DataObj.GlobalAttributes, ...
  DataObj.data, ...
  DataObj.VariableAttributes, ...
  DataObj.Variables, ...
  'calculateMd5Checksum', 1);
end
