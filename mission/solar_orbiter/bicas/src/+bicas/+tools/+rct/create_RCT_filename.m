% Create RCT filename (time-stamped) following official filenaming convention.
%
%
% RCT FILENAME CONVENTION
% =======================
% NOTE: The official RPW RCT filenaming convention has changed multiple times and can be confused.
% --
% See implementation for comments.
% See comments for settings PROCESSING.RCT_REGEXP.* (all RCTs), in
% bicas.create_default_BSO().
%
%
% """"""""
% 4.2.3 File naming
% The RPW CAL file shall comply the following file naming convention [AD1]:
% solo_[LEVEL]_[Descriptor]_[Datetime]_V[CALIBRATION_VERSION].cdf
% Where:
% - [Descriptor] is the prefix of the Descriptor CDF global attribute value (before “>”).
% It shall be of the form rpw-[equipment]-[*], where is the name of the RPW
% equipment, e.g., “tds”, “lfr”, “bia”, “scm”, ... and [*] is an optional field that can be
% used to specify the content of the file (e.g., “bias-f0”)
% - [LEVEL] is the RCT data processing level as defined in [AD1]. It shall be “CAL”.
% - [Datetime] gives the date when the RCT CDF file has been released. The format is
% “yyyymmdd”.
% - [CALIBRATION_VERSION] is a 2-digit number of the version of the calibration table
% file (e.g., “01”).
% IMPORTANT:
% - Each field in the file name shall be separated by an underscore “_” (use hyphens “-“
% inside each field to distinguish sub-strings).
% - The RPW CAL file shall use lower case convention, except for the level and the version.
%
% 4.2.4 Data versioning
% The version of a RPW CAL CDF file is a 2-digit number, which shall be iterated each time a
% new version of the file is released for a given descriptor and date.
% Initial version number shall be “01”.
% In the RPW CAL filename, the version number shall appear with the “V” prefix in capital letter
% (e.g., “V02”).
% The value of the CALIBRATION_VERSION global attribute (see next section) shall be set
% with the file version number.
% """"""""
% /ROC-PRO-PIP-ICD-00037-LES, RCS ICD, version 1/7
% NOTE: The "Datetime" description appears to be wrong and contradicts
% SOL-SGS-TN-0009, "Metadata Definition for Solar Orbiter Science Data", version
% 2/6, Section 2.1.3.5
%
%
% """"""""
% 2.1.3.5 Calibration and Ancillary Files
% Calibration and Ancillary Data will follow the same standard. Calibration files will have CAL
% included as the level; further subcategories can be added as appropriate in descriptor, separated by
% hyphens. E.g.
%              solo_CAL_mag-ibs-offset_20191001-20201001_V01.cdf
% """"""""
% /SOL-SGS-TN-0009, "Metadata Definition for Solar Orbiter Science Data",
% version 2/6
%
%
% RETURN VALUES
% =============
% destFilename
% gaCALIBRATION_VERSION
%   GA "CALIBRATION_VERSION"
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function destFilename = create_RCT_filename(DtBegin, DtEnd, versionNbr)

irf.dt.assert_UTC_midnight(DtBegin)
irf.dt.assert_UTC_midnight(DtEnd)

% NOTE: Should not contain seconds.
% gaCALIBRATION_VERSION = char(datetime("now","Format","uuuuMMddHHmm"));

% IMPLEMENTATION NOTE: The official filenaming convention is not followed
% here!! Not sure how to comply with it either (which receiver should the
% BIAS RCT specify?).
S.datasetId          = 'SOLO_CAL_RPW-BIAS';
S.isCdag             = false;
S.Dt1                = DtBegin;
S.Dt2                = DtEnd;
S.timeIntervalFormat = 'DAY_TO_DAY';
S.versionNbr = versionNbr;
S.lesTestStr = [];
S.cneTestStr = [];
destFilename = solo.adm.dsfn.DatasetFilename(S).filename;

end
