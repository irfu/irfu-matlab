function mat_date = epoch2date(isdat_epoch)
%EPOCH2DATE converts ISDAT epoch to MATLAB datenum
%   ISDAT epoch is the number of seconds since 1-Jan-1970 and
%   MATLAB datenum is the number of days since 0-Jan-0000.
%
%   See also DATENUM, DATESTR, TOEPOCH, FROMEPOCH

% Yuri Khotyaintsev, 2003

% 719529 is the number of days from 0-Jan-0000 to 1-Jan-1970
mat_date = double(719529 + double(double(isdat_epoch(:))/double(24 * 3600)));