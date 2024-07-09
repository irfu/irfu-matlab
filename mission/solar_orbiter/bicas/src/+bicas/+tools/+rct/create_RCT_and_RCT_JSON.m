%
% Create RCT and RCT JSON files to be delivered together.
%
% NOTE: Should probably ideally generate RCT and RCT JSON files with the
% official CDF version, i.e. using the BICAS delivery git repo, not irfu-matlab.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-02.
%
function create_RCT_and_RCT_JSON(rctMasterCdfFile, destDir)
% NOTE <2020-11, EJ: Should NOT(!!!) use RCT master CDF at
%                    /nonhome_data/work_files/SOLAR_ORBITER/skeletons_BIAS_RCT/.
%   2020-11-19, EJ: Suspect that faulty RCT master CDF has been removed.
%
% MANUAL CALL
% bicas.tools.rct.create_RCT_and_RCT_JSON('/nonhome_data/work_files/SOLAR_ORBITER/DataPool/SOLO/RPW/CDF/Master/SOLO_CAL_RPW-BIAS_V02.cdf', '/home/erjo/temp/temp/')

rctPath     = bicas.tools.rct.create_RCT(rctMasterCdfFile, destDir);
fprintf(1, 'Created RCT file      "%s"\n', rctPath);

rctJsonPath = bicas.tools.rct.create_RCT_JSON(destDir, irf.fs.get_name(rctPath));
fprintf(1, 'Created RCT JSON file "%s"\n', rctJsonPath);
end
