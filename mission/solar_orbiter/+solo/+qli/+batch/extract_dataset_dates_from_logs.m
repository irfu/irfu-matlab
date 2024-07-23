%
% Given a file path pattern which matches on or multiple log files, select the
% file with the "last" filename if filenames sorted alphabetically. Search that
% one log file for day-long dataset filenames with specified dataset IDs. Return
% the set of unique (starting) timestamps (midnight) for those datasets.
%
% Is primarily intended to be used for selecting dates for which to generate
% quicklooks (QLIs).
%
%
% ARGUMENTS
% =========
% logFileDirPattern
%       String. String pattern for dir() command describing one or multiple log
%       files.
%       NOTE: Must match at least one file (~failsafe).
% dsiCa
%       Cell array of dataset IDs for datasets which should be searched for.
%       NOTE: This excludes any "-cdag" suffix. The function will match both
%       CDAG and non-CDAG files.
%
%
% RETURN VALUES
% =============
% DatasetsDtArray
%       Column array of datetime. The unique dates for the matching filenames
%       found in the selected log file. Beginning of time interval covered by
%       dataset according to filename (midnight for day-long datasets).
%       TimeZone = 'UTCLeapSeconds'.
% logFilePath
%       Path to selected log file.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [DatasetsDtArray, logFilePath] = extract_dataset_dates_from_logs(...
  logFileDirPattern, dsiCa)

% PROPOSAL: Specify dataset filename patterns (not DSIs).
% PROPOSAL: Require at least one DSI.
% PROPOSAL: For every log, collect the log FMD. Only update QLIs (implicated by
%           the datasets mentioned in the log) whose QLI FMD is older than the log FMD.
%   PRO: Rerunning a algorithm will not update the same QLIs again.
%   PROBLEM: Has no way of obtaining the QLI FMDs for selected dates. Must
%            retreive QLI FMDs for all dates (which is somewhat slow).

logFilePath = select_log_file(logFileDirPattern);
s = fileread(logFilePath);

datasetFileNameCa = cell(0, 1);
for i = 1:numel(dsiCa)
  dsi = dsiCa{i};

  assert(strcmp(dsi, upper(dsi)), 'dsi="%s" is not uppercase (convention).', dsi)

  % IMPLEMENTATION NOTE: Important to prevent maximal munch from making matches
  % covering multiple datasets/filenames.
  % Ex: Over multiple rows. ==> Exclude line feed in filename.
  % Ex: On the same row.    ==> Exclude period in filename (except before file suffix).
  %
  % NOTE: solo.adm.parse_dataset_filename()'s support for the unofficial
  % basename extension can cause problems if there are multiple dataset
  % filenames on the same row and one does not exclude e.g. period.
  % NOTE: Must permit filenames with and without "-cdag".
  pattern = sprintf('%s(|-cdag)_[^\\n.]*\\.cdf', dsi);

  % IMPLEMENTATION NOTE: Using case-insensitive reg. exp. matching to handle
  % that (1) dataset filenames contain DSI in lowercase (mostly), and (2) that
  % uppercase dataset IDs have not been historically required in this code.
  % Note: Some dataset filenames actually have mixed case in the dataset ID part
  % (I think) but they are not relevant here (yet).
  matchCa = regexpi(s, pattern, 'match');

  datasetFileNameCa = [datasetFileNameCa; matchCa(:)];
end
% assert(iscolumn(datasetFileNameCa))

DsmdArray = solo.adm.paths_to_DSMD_array(datasetFileNameCa(:));

% IMPLEMENTATION NOTE: Special case for zero matches/DSMDs since
% [DsmdArray.dt1] behaves differently for empty DSMD array (it returns 0x0
% double).
if isempty(DsmdArray)
  Dt1Array = solo.qli.const.EMPTY_DT_ARRAY;
else
  Dt1Array = [DsmdArray.dt1];
end
Dt1Array = unique(Dt1Array);
Dt1Array = sort(Dt1Array);

DatasetsDtArray = Dt1Array(:);

end







% Select log file: File with last filename, if filenames are sorted
% alphabetically.
%
%
% NOTES
% -----
% PROBLEM: It is not obvious how to select log file when there are multiple
% simultaneous log files being built on simultaneously (most relevant when
% running multiple batch BICAS processing runs). One could use (1) the last file
% name (last in alphabetic order) or (2) last file modification timestamp, but
% neither truly solves the problem alone.
%
% IMPLEMENTATION NOTE: One could use file modification date, but this is bad
% when running multiple processes simultaneously since that means multiple files
% are modified simultaneously.
%
% For example, if this function is called after one BICAS batch processing has
% just finished, but another BICAS batch processing is still underway, the
% latter's log file is still continuously updated and may have a later file
% modification date.
function logFilePath = select_log_file(logFileDirPattern)

% FSOI = File System Object Info
FsoiArray = dir(logFileDirPattern);
FsoiArray = FsoiArray(~[FsoiArray.isdir]);

% Require non-zero matching log files
% -----------------------------------
% IMPLEMENTATION NOTE: This is a failsafe against setting the wrong path
% pattern.
if isempty(FsoiArray)
  error('No files match logFileDirPattern=%s.', logFileDirPattern)
end



%[~, iSort] = sort([FsioArray.datenum], 'ascend');   % Sort by log file modification date.
[~, iSort]  = sort({FsoiArray.name});   % Sort by log filename (not entire path).
Fsoi        = FsoiArray(iSort(end));
logFilePath = fullfile(Fsoi.folder, Fsoi.name);

end