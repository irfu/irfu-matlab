function isdat_epoch = date2epoch(mat_date)
%DATE2EPOCH converts MATLAB datenum to ISDAT epoch
%   ISDAT epoch is the number of seconds since 1-Jan-1970 and
%   MATLAB datenum is the number of days since 0-Jan-0000.
%
%   See also EPOCH2DATE, DATENUM, DATESTR, TOEPOCH, FROMEPOCH

% Yuri Khotyaintsev, 2004

% 719529 is the number of days from 0-Jan-0000 to 1-Jan-1970
isdat_epoch = double(mat_date(:) - 719529)*double(24 * 3600);
