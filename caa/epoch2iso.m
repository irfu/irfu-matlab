function out = epoch2iso(t)
%EPOCH2ISO converts ISDAT epoch to ISO time string
% ISDAT epoch is the number of seconds since 1-Jan-1970 and
% ISO time string has a format yyyy-mm-ddTHH:MM:ss.wwwwwwZ
% as described in the CEF data file syntax
% http://www.space-plasma.qmul.ac.uk/csds/welcome.html
%
% t_iso_s = iso2epoch(t_epoch)
%
%   See also EPOCH2DATE, ISO2EPOCH, TOEPOCH, FROMEPOCH
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev

% We need to do all this because DATESTR rounds seconds
d = fromepoch(t);

for j=2:5, s1(j-1) = {add_zero(d(:,j),num2str(d(:,j)))}; end

% Take care about seconds separately
s2 = add_zero(d(:,6),num2str(d(:,6),'%6f'));

sZ = s2(:,1); sZ(:) = 'Z';
sT = sZ; sT(:) = 'T';
sdash = sZ; sdash(:) = '-';
scol = sZ; scol(:) = ':';

out = [num2str(d(:,1)) sdash s1{1} sdash s1{2} sT s1{3} scol s1{4} scol s2 sZ];

% Help function to insert zeros for numbers containing only one digit 
function out = add_zero(d,s)
out = s;
ii = find(d<10);
if ~isempty(ii)
	ss = s(ii,1);
	ss(:) = '0';
	if length(ii)==length(d)
		% Add to all lines
		out = [ss s];
	else
		out(ii,:) = [ss s(ii,2:end)];
	end
end
