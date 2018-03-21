function out = iso2epoch(s)
%ISO2EPOCH convert ISO time string to ISDAT epoch
%
% t_epoch = iso2epoch(t_iso_s)
%
% Input:
%	s - string (may be a matrix) in format:
%		yyyy-mm-ddTHH:MM:ss.mmmZ	// ISO Time
%		yyyy-mm-ddTHH:MM:ss.mmmuuuZ	// Extended ISO Time
%   yyyy-mm-ddTHH:MM:ss.mmmuuunnnZ // VERY Extended ISO Time
%
%    See also EPOCH2ISO, TOEPOCH, FROMEPOCH

% Copyright 2004,2007 Yuri Khotyaintsev

mask = '%4d-%2d-%2dT%2d:%2d:%fZ';

% If we have multiple rows, we need to turn the matrix
if min(size(s))>1
	if any(size(s,2) == [30 27 24 20]), s=s'; end
	n_column = size(s,2);
else, n_column = 1;
end

a = sscanf(s,mask);

N = length(a)/6;
if N~=fix(N) || N~=n_column, irf.log('warning','something is wrong with input'), end
a = reshape(a,6,fix(N));
a = a';
out = toepoch([a(:,1) a(:,2) a(:,3) a(:,4) a(:,5) a(:,6)]);
