%
% Given a file path pattern which matches on or multiple log files, select the
% file with the "last" filename if filenames sorted alphabetically. Search that
% one log file for day-long dataset filenames with specified dataset IDs. Return
% the set of unique (starting) timestamps for those datasets.
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
%       NOTE: Must match at least one file.
% datasetIdCa
%       Cell array of dataset IDs for datasets which should be searched for.
%       NOTE: This excludes any "-cdag" suffix. The code will match both CDAG
%       and non-CDAG files.
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
function [DatasetsDtArray, logFilePath] = extract_dataset_dates_from_logs(logFileDirPattern, datasetIdCa)
  % PROPOSAL: Require correct case for dataset ID?
  % PROPOSAL: Specify filename patterns (not dataset IDs).
  % PROPOSAL: Require at least one dataset ID.
  %
  % NOTE: Must handle zero matches.

  % FSOI File System Object Info
  FsoiArray = dir(logFileDirPattern);
  FsoiArray = FsoiArray(~[FsoiArray.isdir]);
  assert(~isempty(FsoiArray), 'No files match logFileDirPattern=%s.')

  % Select log file: File with last filename, if filenames are sorted
  % -----------------------------------------------------------------
  % IMPLEMENTATION NOTE: Could use file modification date, but this is bad
  % when running multiple processes simultaneously since that means multiple
  % files are modified simultaneously.

  % For example, if this function is called after one BICAS batch processing has
  % just finished, but another BICAS batch processing is still underway, the
  % latter's log file is still continuously updated and may have a later file
  % modification date.
  %[~, iSort] = sort([FsioArray.datenum], 'ascend');   % Sort by file modification date.
  [~, iSort]  = sort({FsoiArray.name});   % Sort by filename (not entire path).
  Fsoi        = FsoiArray(iSort(end));
  logFilePath = fullfile(Fsoi.folder, Fsoi.name);

  s = fileread(logFilePath);

  datasetFileNameCa = cell(0, 1);
  for i = 1:numel(datasetIdCa)
    datasetId = datasetIdCa{i};
    % IMPLEMENTATION NOTE: Important to prevent maximal munch from making
    % matches covering multiple datasets/filenames.
    % Ex: Over multiple rows. ==> Exclude line feed in filename.
    % Ex: On the same row.    ==> Exclude period in filename (except before file suffix).
    % NOTE: solo.adm.parse_dataset_filename()'s "unofficial" basename extension
    % can cause problems if there are multiple dataset filenames on the same row
    % and one does not exclude e.g. period.
    % NOTE: Must permit filenames with and without "-cdag".
    pattern = sprintf('%s(|-cdag)_[^\\n.]*\\.cdf', datasetId);

    matchCa = regexpi(s, pattern, 'match');   % NOTE: Case-insensitive.

    datasetFileNameCa = [datasetFileNameCa; matchCa(:)];
  end
  % assert(iscolumn(datasetFileNameCa))

  DsmdArray = solo.adm.paths_to_DSMD_array(datasetFileNameCa(:));

  % IMPLEMENTATION NOTE: Special case for zero matches/DSMDs since
  % [DsmdArray.dt1] behaves differently for empty DSMD array (it returns 0x0
  % double).
  if isempty(DsmdArray)
    Dt1Array = datetime.empty(0, 1);
  else
    Dt1Array = [DsmdArray.dt1];
  end
  Dt1Array = unique(Dt1Array);
  Dt1Array = sort(Dt1Array);

  DatasetsDtArray = Dt1Array(:);
end
