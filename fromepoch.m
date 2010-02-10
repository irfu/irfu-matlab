function full_time = fromepoch(second)
%fromepoch - Convert seconds since 1970 to [YYYY MM DD hh mm dd] time format.
%
% timevec = fromepoch(epochsec) converts the time epochsec, which
%   is given as seconds since the epoch 1 Jan 1970, to a numerical vector
%   in the form timevec=[year mon day hour min sec].
%   The seconds since epoch time format is the time specification
%   used by the ISDAT system. To convert into Matlabs standard time
%   format use DATENUM.
%
%   Examples:
%     tv = fromepoch(0)          returns tv = [1970 1 1 0 0 0],
%     mt = datenum(fromepoch(0)) returns mt = 719529,
%     ts = datestr(datenum(fromepoch(0))) returns ts='01-Jan-1970'.
%
%   See also: TOEPOCH, ISGETDATALITE, DATENUM, DATESTR.
%
% $Id$
%

% Yuri Khotyaintsev, 2003
t = datevec(epoch2date(fix(double(second(:)))));

% correct fractions of second. This actually preserves 
% accuracy ~1e-6 sec for year 2004.
t(:,6) = t(:,6) + double(second(:)) - fix(double(second(:))); 

full_time = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

