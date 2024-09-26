%
% Given a ROC SolO HK (RSH) filename, return a parsed date and version number.
%
%
% ARGUMENTS
% =========
% filename
%       Ex: solo_HK_platform_20210901_V02.xml
%
%
% RETURN VALUES
% =============
% isRsh
%       Logical. Whether filename is recognized as an RSH file at all. Can be
%       used to detect whether filename is a ROC SolO HK file or not.
% Dt
%       datetime object for day.
%       NaT if filename is not parsable.
% versionNbr
%       Non-negative integer number.
%       NaN if filename is not parsable.
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-09-08.
%
function [isRsh, Dt, versionNbr] = parse_RSH_filename(filename)
% NOTE: Permits version strings of length >=2 (not just exacly 2).
[subStrCa, ~, isPerfectMatch] = irf.str.regexp_str_parts(...
  filename, ...
  {'solo_HK_platform_', '[0-9]{4}', '[0-1][0-9]', '[0-3][0-9]', ...
  '_V', '[0-9]{2,}', '\.xml'}, ...
  'PERMIT_NON_MATCH');

if ~isPerfectMatch
  isRsh      = false;
  Dt         = NaT;
  versionNbr = NaN;
else
  yearStr  = subStrCa{2};
  monthStr = subStrCa{3};
  dayStr   = subStrCa{4};
  verStr   = subStrCa{6};

  isRsh = true;
  Dt         = datetime([yearStr, '-', monthStr, '-', dayStr]);
  versionNbr = str2double(verStr);
end

end
