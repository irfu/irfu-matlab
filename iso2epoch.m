function out = iso2epoch(s)
%ISO2EPOCH convert ISO time string to ISDAT epoch
%
% t_epoch = iso2epoch(t_iso_s)
%
% Input:
%	s - string (may be a matrix) in format:
%		yyyy-mm-ddTHH:MM:ss.wwwZ	// ISO Time
%		yyyy-mm-ddTHH:MM:ss.wwwwwwZ	// Extended ISO Time
%
%    See also EPOCH2ISO, TOEPOCH, FROMEPOCH
%
% $Id$

% Copyright 2004,2007 Yuri Khotyaintsev

mask = '%4d-%2d-%2dT%2d:%2d:%fZ';

% If we have multiple rows, we need to turn the matrix
if min(size(s))>1
	if size(s,2)==27 || size(s,2)==24 || size(s,2)==20, s=s'; end
	n_column = size(s,2);
else n_column = 1;
end

a = sscanf(s,mask);

N = length(a)/6;
if N~=fix(N) || N~=n_column, disp('something is wrong with input'), end
a = reshape(a,6,fix(N));
a = a';
out = toepoch([a(:,1) a(:,2) a(:,3) a(:,4) a(:,5) a(:,6)]);
