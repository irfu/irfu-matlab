function out = epoch2iso(t)
%EPOCH2ISO converts ISDAT epoch to ISO time string
% ISDAT epoch is the number of seconds since 1-Jan-1970 and
% ISO time string has a format yyyy-mm-ddTHH:MM:ss.wwwwwwZ
% as described in the CEF data file syntax
% http://www.space-plasma.qmul.ac.uk/csds/welcome.html
%
%   See also epoch2date, datenum, datestr, toepoch, fromepoch
%
% $Id$

s1=datestr(epoch2date(t),'yyyy-mm-ddTHH:MM:SS');
s2=zeros(length(t),8);
% take care about fraction of second
for j=1:length(t)
	ss = sprintf('%6f',t(j)-fix(t(j)));
	s2(j,:) = [ss(1,2:8) 'Z'];
end
out = [s1 s2];
