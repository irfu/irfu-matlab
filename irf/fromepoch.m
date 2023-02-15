function full_time = fromepoch(second)
%FROMEPOCH - Convert seconds since 1970 to [YYYY MM DD hh mm dd] format
%
% timevec = fromepoch(epochsec)
%   Converts the time epochsec, which is given as seconds since
%   the epoch 1 Jan 1970, to a numerical vector in the form
%   timevec=[year mon day hour min sec].
%   The seconds since epoch time format is the time specification
%   used by the ISDAT system. To convert into Matlabs standard time
%   format use DATENUM.
%
%   Examples:
%     tv = fromepoch(0)          returns tv = [1970 1 1 0 0 0],
%     mt = datenum(fromepoch(0)) returns mt = 719529,
%     ts = datestr(datenum(fromepoch(0))) returns ts='01-Jan-1970'.
%
%   See also: TOEPOCH, EPOCH2ISO, IRF_TIME, DATENUM, DATESTR

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if(verLessThan('matlab','8.4'))
  % Offset by 0.1 seconds to avoid rounding problems in
  % Matlab versions prior to R2009a (7.8)
  % See http://www.mathworks.com/support/bugreports/485847
  t = datevec(epoch2date(fix(double(second(:)))+0.1));
  t(:,6)=fix(t(:,6));
  % correct fractions of second. This preserves double-precision accuracy
  % using isdat epoch (1970-01-01) rather than matlab epoch (year 0).
  % Resulting accuracy ~1e-6 sec for year 2004.
  t(:,6) = t(:,6) + double(second(:)) - fix(double(second(:)));
  full_time = t;
else
  % Use built in Matlab function (as of Matlab 8.4, R2014b)
  t = datetime(second(:),'ConvertFrom','posixtime');
  full_time=[t.Year, t.Month, t.Day, t.Hour, t.Minute, t.Second];
end

end
