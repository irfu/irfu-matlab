function out = epoch2iso(t,fmt)
%EPOCH2ISO   Convert ISDAT epoch to ISO time string
%
% ISDAT epoch is the number of seconds since 1-Jan-1970 and
% ISO time string has a format yyyy-mm-ddTHH:MM:ss.wwwwwwZ
% as described in the CEF data file syntax
% http://www.space-plasma.qmul.ac.uk/csds/welcome.html
%
% t_iso_s = epoch2iso(t_epoch,[fmt])
%   fmt - 0: long (default), 1: short
%
%   See also EPOCH2DATE, ISO2EPOCH, TOEPOCH, FROMEPOCH
%
% $Id$

% Copyright 2004-2007 Yuri Khotyaintsev

if nargin<2, fmt = 0; end

if length(t)<5
	% We need to do all this because DATESTR rounds seconds
	d = fromepoch(t);
	
	for j=2:5, s1(j-1) = {add_zero(d(:,j),num2str(d(:,j),'%d'))}; end
	
	% Take care about seconds separately
	s2 = add_zero(d(:,6),num2str(d(:,6),'%6f'));
	if fmt, s2 = s2(:,1:6); end
	
	sZ = s2(:,1); sZ(:) = 'Z';
	sT = sZ; sT(:) = 'T';
	sdash = sZ; sdash(:) = '-';
	scol = sZ; scol(:) = ':';
	
	out = [num2str(d(:,1)) sdash s1{1} sdash s1{2} sT s1{3} scol s1{4} scol s2 sZ];
else
	% This approach is faster for data with many samples per minute, as we run 
	% from epoch only once per minute
	out = char(zeros(length(t),27-fmt*3));
	out(:,[5 8]) = '-';
	out(:,11) = 'T';
	out(:,[14 17]) = ':';
	out(:,end) = 'Z';
	
	tss = fromepoch(t(1));
	tee = fromepoch(t(end));
	
	mins = toepoch([tss(1:5) 0]):60:toepoch([tee(1:5) 0]);
	d = fromepoch(mins);
	
	s1 = {'', '', '', '',''}; sl=[0 4; 5 2; 8 2; 11 2; 14 2];
	j_start = 0;
	
	for j=1:5
		if d(1,j)==d(end,j)
			ss = add_zero(d(1,j),num2str(d(1,j),'%d'));
			for jj=1:sl(j,2), out(:,sl(j,1)+jj) = ss(jj); end
		else
			j_start = j; 
			for jj=j:5, s1(jj) = {add_zero(d(:,jj),num2str(d(:,jj),'%d'))}; end
			break
		end
	end
	
	for j=1:length(mins)
		if j==length(mins), ii = find(t>=mins(j));
        else ii = find(t>=mins(j) & t<mins(j+1));
		end;
		if isempty(ii), continue, end
		if j_start
			for kk=j_start:5
				for jj=1:sl(kk,2), out(ii,sl(kk,1)+jj) = s1{kk}(j,jj); end
			end
		end
		s2 = add_zero(t(ii)-mins(j),num2str(t(ii)-mins(j),'%6f'));
		if fmt, s2 = s2(:,1:6); end
		out(ii,18:end-1) = s2;
	end
end

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
		out(ii,:) = [ss s(ii,1:end-1)];
	end
end

