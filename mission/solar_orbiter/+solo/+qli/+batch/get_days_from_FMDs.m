%
% Generate array of dates derived from file modification dates, which should
% indicate datasets which are newer than the corresponding quicklooks.
%
%
% ARGUMENTS
% =========
% Settings
%       Struct with fields for relevant settings.
%
%
% RETURN VALUES
% =============
% DaysDtArray
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function DaysDtArray = get_days_from_FMDs(Settings, varargin)

assert(numel(varargin) == 0, 'Illegal number of additional arguments.')

dsiCa = [solo.qli.batch.const.SOURCE_DSI_DICT.values{:}]';

DaysDtArray = solo.qli.batch.fmd.get_days_from_FMDs(...
  Settings.datasetDirsCa, ...
  Settings.qliDir, ...
  dsiCa);
end
