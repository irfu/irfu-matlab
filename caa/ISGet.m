function [t,d] = ISGet(db_s,st,dt,cli,ins,sig,sen,cha,par)
% ISGet help wrapper around isGetDataLite
% [t,d] = ISGet(db_s,st,dt,cli,ins,sig,sen,cha,par)
%
% Last 3 input arguments can be skipped
%
% $Id$
%
% see also isGetDataLite


if nargin < 9, par = ' '; end
if nargin < 8, cha = ' '; end
if nargin < 7, sen = ' '; end

p = tokenize(db_s,'|');

for i=1:length(p)
	% reset errors, otherwise try/catch fails
	lasterr('')
	try
		dbase = Mat_DbOpen(p{i});
		disp(lasterr)

		[t, d] = isGetDataLite(dbase, st, dt, ...
		'Cluster',num2str(cli),ins,sig,sen,cha,par);

		if exist('dbase'), Mat_DbClose(dbase), clear dbase, end

		if ~isempty(d), return, end

	catch
		warning('ISDAT:getData',...
		'Error getting data from ISDAT server %s: %s', p{i},lasterr)
		t = []; d = [];
	end
end
