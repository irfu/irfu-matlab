function [t,d] = caa_is_get(db_s,st,dt,cli,ins,sig,sen,cha,par)
% caa_is_get help wrapper around isGetDataLite
% [t,d] = caa_is_get(db_s,st,dt,cli,ins,sig,sen,cha,par)
%
% INPUT
%  db_s - string defining dbh data base, if several separate by '|', e.g. 'disco:10|sanna:1'
%       - if db_s is not string assume that it is handler of dbh connection, see Mat_DbOpen
%  st   - start time vector, e.g. [2002 2 6 0 0 0]
%  dt   - time interval in seconds
%  cli,.. - client, instrument, signal, ...
%
% Last 3 input arguments can be skipped
%
% Eample: [t,mlt]=caa_is_get('disco:10',[2002 2 6 0 0 0],100,4,'ephemeris', 'mlt');
%
% $Id$
% see also isGetDataLite


if nargin < 9, par = ' '; end
if nargin < 8, cha = ' '; end
if nargin < 7, sen = ' '; end

if ischar(db_s) % db_s is dbh database list
  p = tokenize(db_s,'|');
  db_s_is_database_list=1;
else            % assume db_s is dbh handler
  p= 'null';
  db_s_is_database_list=0;
end

for i=1:length(p)
	% reset errors, otherwise try/catch fails
	lasterr('')
	try
    if db_s_is_database_list,
  		dbase = Mat_DbOpen(p{i});
    else % assume db_s is databse handler
      dbase=db_s;
    end
	if ~isempty(lasterr), irf_log('dsrc',lasterr), end

		[t, d] = isGetDataLite(dbase, st, dt, ...
		'Cluster',num2str(cli),ins,sig,sen,cha,par);

		if db_s_is_database_list,
      if exist('dbase'), Mat_DbClose(dbase), clear dbase, end
    end

		if ~isempty(d), return, end

	catch
		irf_log('dsrc',...
		'Error getting data from ISDAT server %s: %s', p{i},lasterr)
		t = []; d = [];
	end
end
