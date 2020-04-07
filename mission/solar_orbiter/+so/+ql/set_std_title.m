%
% Set a standardized plot title (call "title).
%
% ARGUMENTS
% =========
% plotTypeStr : Top string, e.g. "LFR CWF".
% filePath    : Path to file. Will only use the filename but permits submitting entire path for convenience.
%
function set_std_title(plotTypeStr, filePath, hTopAxes)
    assert(isscalar(hTopAxes))
    
    [~, basename, suffix] = fileparts(filePath);
    filename = [basename, suffix];
    
    labelTimestamp = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
    title(hTopAxes, {plotTypeStr, so.ql.escape_str(sprintf('Plot time: %s, %s', labelTimestamp, filename))})
end