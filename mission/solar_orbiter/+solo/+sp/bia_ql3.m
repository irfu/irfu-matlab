%
% Main function for software "BIA_QL3" for creating BIAS summary plots. The
% software is intended to be delivered to, and used by, LESIA/ROC.
%
% See irfu-matlab/mission/solar_orbiter/bia_ql3/README.TXT for more
% documentation of "BIA_QL3".
%
% It therefore:
% ** has a simpler interface specified by ROC
% ** is meant to be called from bash wrapper script
%    irfu-matlab/mission/solar_orbiter/bia_ql3/roc/bia_ql3
%
%
% INTERFACE REQUIREMENTS FOR THIS FUNCTION
% ========================================
% LESIA/ROC requries a certain OS SHELL interface. See
%   https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/42
%   > shellExecutable YYYYMMDD input_bia_path input_lfr_wf_e_path output_dir log_dir
% This function's interface mirrors the bash/OS script interface to simplify the
% bash wrapper script. The function is therefore not very versatile.
%
%
% NOTES
% =====
% * Makes assumptions on filenames due to using ~globbing patterns.
% * Assumes there is only one matching dataset in the respective directories
%   (non-recursive).
%   NOTE: Can handle the absence of matching datasets/directories, but fails on
%   multiple matching datasets (assertion).
% * Works for both CDAG/non-CDAG filenames.
% * Code does not use other "library" code, in particular not solo.adm for
%   parsing filenames.
%
%
% ARGUMENTS
% =========
% yyyyMmDdStr
%       String yyyymmdd. yyyy=year (4 digits), mm=Month (2 digits), dd=day of
%       month (2 digits).
% hkBiaDir
%       Path to (day) directory with BIAS HK datasets directly underneath it.
% lfrWfDirPath
%       Path to (day) directory with L2 LFR CWF + SWF datasets directly
%       underneath it.
% outputDir
%       Path to output directory.
% --
% NOTE: The code fails gracefully and continues if the HK or CWF/SWF datasets or
% directories do not exist. This is useful since datasets as well as directories
% might legitimately not exist for days without data.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-29.
%
function bia_ql3(yyyyMmDdStr, hkBiaDir, lfrWfDirPath, outputDir)
% PROPOSAL: Eliminate use of glob.attr. Do not read file an extra time.
%   PRO: DATASET_ID is known when calling plot_save_SP_pattern().
%   CON: Dataset data version has to be extracted from filename.

% IMPLEMENTATION NOTE: Needed to get paths correct when code is being called
% from bash script.
% IMPLEMENTATION NOTE: Not using irf(...) commands that produce file
%   ~/.matlab_datastore_<hostname>
% since this causes some kind of problems for ROC (Quynh Nhu NGUYEN
% when running BICAS in parallel).
% See https://gitlab.obspm.fr/ROC/RCS/BIA_QL3/-/issues/3/.
irf('check_path');

irf.log('notice')    % Set log level.



% ASSERTIONS
assert(ischar(yyyyMmDdStr))
assert(numel(yyyyMmDdStr) == 8, ...
  'yyyyMmDdStr="%s" does not have exactly 8 characters.', yyyyMmDdStr)
% NOTE: Does not necessarily give a perfect error message, but does print
% the flawed string.
irf.assert.castring_regexp(yyyyMmDdStr, '20[1-9][0-9][0-1][0-9][0-3][0-9]')

plot_save_SP_pattern(outputDir, hkBiaDir,     'solo_HK_rpw-bia_%s_V*.cdf',             yyyyMmDdStr)
plot_save_SP_pattern(outputDir, lfrWfDirPath, 'solo_L2_rpw-lfr-surv-cwf-e*_%s_V*.cdf', yyyyMmDdStr)
plot_save_SP_pattern(outputDir, lfrWfDirPath, 'solo_L2_rpw-lfr-surv-swf-e*_%s_V*.cdf', yyyyMmDdStr)
end



function plot_save_SP_pattern(outputDir, datasetDir, filenamePattern, yyyyMmDdStr)

datasetPathPattern = fullfile(...
  datasetDir, ...
  sprintf(filenamePattern, yyyyMmDdStr));

% NOTE: dir() fails gracefully/continues if the pattern has no matches, even
% if the directory itself does not exist. This is useful since datasets as
% well as directories might legitimately not exist for days without data.
Fi = dir(datasetPathPattern);
if isscalar(Fi)
  datasetPath = fullfile(Fi.folder, Fi.name);
elseif isempty(Fi)
  % Do nothing.
  irf.log('w', sprintf('Can not find any dataset matching "%s".\n', datasetPathPattern))
  return    % EXIT
else
  error('Found multiple files matching pattern "%s".', datasetPathPattern)
end

irf.log('n', sprintf('Creating summary plot for dataset "%s".', datasetPath))



%===============================================================================
% Read CDF to obtain DSI and dataset version from global attributes
% -----------------------------------------------------------------
% NOTE: Not using data for plotting. Only reading this to avoid having to parse
% the filename (and rely on filenaming conventions).
%===============================================================================
Do = dataobj(datasetPath);
dsi           = Do.GlobalAttributes.Dataset_ID{1};
dataVersionGa = Do.GlobalAttributes.Data_version{1};

% Normalize: dataVersionNbr := normalized dataVersionGa
if isnumeric(dataVersionGa)
  dataVersionNbr = dataVersionGa;
elseif ischar(dataVersionGa)
  dataVersionNbr = str2double(dataVersionGa);
else
  error(...
    ['Can not determine dataset version from global', ...
    ' attribute Data_version in file "%s"'], datasetPath)
end



%=====================================
% dateVec3 := reformatted yyyyMmDdStr
%=====================================
dateVec3 = [...
  str2double(yyyyMmDdStr(1:4)), ...
  str2double(yyyyMmDdStr(5:6)), ...
  str2double(yyyyMmDdStr(7:8))];



plot_save_SP_file(outputDir, datasetPath, dsi, dateVec3, dataVersionNbr)
end



function plot_save_SP_file(outputDir, datasetPath, dsi, dateVec3, dataVersionNbr)

switch(dsi)

  case 'SOLO_HK_RPW-BIA'
    plotFunc = @solo.sp.plot_HK;

  case {'SOLO_L2_RPW-LFR-SBM1-CWF-E', ...
      'SOLO_L2_RPW-LFR-SBM2-CWF-E', ...
      'SOLO_L2_RPW-LFR-SURV-CWF-E'}
    plotFunc = @solo.sp.plot_LFR_CWF;

  case 'SOLO_L2_RPW-LFR-SURV-SWF-E'
    plotFunc = @solo.sp.plot_LFR_SWF;

  otherwise
    error('Can not plot dataset "%s" using DATASET_ID="%s".', ...
      datasetPath, dsi)
end



spFilename = solo.sp.create_SP_filename(dsi, dateVec3, dataVersionNbr);
spFilePath = fullfile(outputDir, spFilename);



% IMPLEMENTATION NOTE: The solo.sp.plot_* functions create their own
% figures. Can therefore not create the figure beforehand and set
% visibility=off.
irf.log('n', sprintf('Creating in-memory summary plot for "%s"', spFilePath))
t = tic();
plotFunc(datasetPath);
wallTimeSec = toc(t);
irf.log('n', sprintf('Wall time to create in-memory summary plot: %.2f [s]', wallTimeSec))



hFig = gcf();
irf.log('n', sprintf('Saving in-memory plot to file "%s"', spFilePath))
solo.sp.save_SP_figure(spFilePath, hFig)
close(hFig)
end
