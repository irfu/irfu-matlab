%
% Create RCT and RCT JSON files to be delivered together.
%
% NOTE: Should probably ideally generate RCT and RCT JSON files with the
% official CDF version, i.e. using the BICAS delivery git repo, not irfu-matlab.
%
%
% ARGUMENTS
% =========
% rctVersionNbr
%       RCT version number.
%       NOTE: Shall presumably be incremented for new RCT versions which have
%       the same time interval of validity, just like datasets.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-02.
%
function create_RCT_and_RCT_JSON(rctMasterCdfFile, destDir, rctVersionNbr)
% bicas.tools.rct.create_RCT_and_RCT_JSON('/nonhome_data/work_files/SOLAR_ORBITER/DataPool/SOLO/RPW/CDF/Master/SOLO_CAL_RPW-BIAS_V03.cdf', '~/temp/temp', 1)

% Time interval for which the RCT is valid. Should cover the entire mission.
% --------------------------------------------------------------------------
% NOTE: Solar Orbiter was launched "10 February 2020, 04:03 UTC" (Wikipedia),
%       but 2022-02-09 local Florida time.
% NOTE: Cf. solo_CAL_rpw-scm_20200210-20990101_V10.cdf which uses the SolO
%       launch date and presumably tries to cover the entire mission too.
% NOTE: The end of the time interval in END_DT is represented differently in
%       (1) the RCT filename and (2) the JSON file content.
%       since the RCT filename (a) only specifies dates (not time) and
%       (b) the RCT filename end date is inclusive.
%       Ex:     2100-01-01T00:00:00Z in JSON file content.
%           ==> 2099-12-31 in the RCT filename.
% IMPORTANT NOTE: One should presumably not change the time interval
%                 unnecessarily, since it is represented in the filename and
%                 since the version number is (presumably) tied to the exact
%                 chosen time interval.
BEGIN_DT = datetime('2020-02-10T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
END_DT   = datetime('2100-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');



% IMPLEMENTATION NOTE: Log messages run AFTER the completion of sub-tasks, since
% one wants to log the paths to the created files.

rctPath     = bicas.tools.rct.create_RCT(rctMasterCdfFile, destDir, BEGIN_DT, END_DT, rctVersionNbr);
fprintf(1, 'Created RCT file      "%s"\n', rctPath);

rctJsonPath = bicas.tools.rct.create_RCT_JSON(destDir, irf.fs.get_name(rctPath), BEGIN_DT, END_DT);
fprintf(1, 'Created RCT JSON file "%s"\n', rctJsonPath);
end
