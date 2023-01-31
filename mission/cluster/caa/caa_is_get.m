function [t,d] = caa_is_get(db_s,st,dt,cli,ins,sig,sen,cha,par,units)
% CAA_IS_GET  get Cluster data from ISDAT database
%
% [time,data] = caa_is_get(db_s,st,dt,cl_id,ins,sig,[sen,cha,par,units])
%  get data from ISDAT database
%
% OUTPUT
%  data  - data matrix in AV format,first column is time (epoch)
% INPUT
%  db_s  - string defining dbh data base
%  st    - start time epoch, e.g. toepoch([2002 2 6 0 0 0])
%  dt    - time interval in seconds
%  cl_id - Cluster ID 1..4
%  ins,sig,sen,cha,par - data request parameters
%  units - 'tm' or 'phys' (default)
%
% In case ISDAT server returned DbBAD_INTERNAL (may be caused by the server
% crash) program will sleep for some time and retry the request.
%
% Example:
%  [t,d] = caa_is_get('disco:10',toepoch([2002 2 6 0 0 0]),100,...
%          4,'ephemeris','mlt');
%
% See also isGetDataLite
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% TODO: we may try to open persistent connection when ISDAT will get stable
% TODO: can be use Mat_DbErrorString

MAXOPENATTEMPTS = 10;  % Number of times to retry in case open fails
MAXNRET = 2;  % Number of times to retry in case of internal server error
SLEEPINT = 5; % Number of seconds to sleep when waiting for server restart

narginchk(6,10)
if nargin < 10, units = ''; end
if nargin < 9, par = ' '; end
if nargin < 8, cha = ' '; end
if nargin < 7, sen = ' '; end

% Check imput
if ~ischar(db_s), error('DB_S must be a string'), end
if ~isnumeric(st) || length(st) ~= 1
  error('ST must be a double number')
end
if ~isnumeric(st), error('DT must be a number'), end
if ~isnumeric(cli) || cli<0 || cli>4, error('CL_ID must be a number 1..4'), end
if ~ischar(sig), error('SIG must be a string'), end
if ~ischar(sen), error('SEN must be a string'), end
if ~ischar(cha), error('CHA must be a string'), end
if ~ischar(par), error('PAR must be a string'), end

nret = 0;
while nret<MAXNRET
  nret = nret + 1;
  
  openAttempts=1;
  while (openAttempts ~= 0)
    try
      dbase = Mat_DbOpen(db_s);
      openAttempts=0;
    catch
      if (openAttempts>MAXOPENATTEMPTS)
        error('ISDAT:Mat_DbOpen',['ISDAT : error opening database ' db_s])
      end
      openAttempts=openAttempts+1;
      pause(0.01*openAttempts);
    end
  end
  
  if ~isempty(units)
    [t, d, iserr] = isGetDataLite(dbase, st, dt, ...
      'Cluster',num2str(cli),ins,sig,sen,cha,par,units);
  else
    [t, d, iserr] = isGetDataLite(dbase, st, dt, ...
      'Cluster',num2str(cli),ins,sig,sen,cha,par);
  end
  
  Mat_DbClose(dbase);
  d = double(d);
  
  % If we got DbBAD_INTERNAL we suppose it was cause by the server crash,
  % so we sleep for SLEEPINT sec and try again
  if iserr ~= 15, break, end
  pause(SLEEPINT)
end

if iserr==0 && ~isempty(d), return, end

if is_nodata_err(iserr)
  irf_log('dsrc',['ISDAT : ' is_err2str(iserr)])
elseif isempty(d)
  irf_log('dsrc','ISDAT : empty return')
else
  error('ISDAT:Mat_DbGetDataLite',['ISDAT : ' is_err2str(iserr)])
end


function res = is_err2str(err_id)
% Help function to get error message from ISDAT error code
% Isdat_root/lib/Db/StringFunc.c : DbErrorString

switch err_id
  case 0  % DbSUCCESS
    res = 'no error';
  case 1  % DbBAD_DROP
    res = 'data drop';
  case 2  % DbBAD_TIME
    res = 'no data';
  case 3  % DbBAD_PROJECT
    res = 'bad project';
  case 4  % DbBAD_MEMBER
    res = 'bad member';
  case 5  % DbBAD_EOF
    res = 'end of file';
  case 6  % DbBAD_INSTRUMENT
    res = 'bad instrument';
  case 7  % DbBAD_SENSOR
    res = 'bad sensor';
  case 8  % DbBAD_SIGNAL
    res = 'bad signal';
  case 9  % DbBAD_CHANNEL
    res = 'bad channel';
  case 10 % DbBAD_PARAMETER
    res = 'bad parameter';
  case 11 % DbBAD_UNITS
    res = 'bad units';
  case 12 % DbBAD_REDUCTION
    res = 'bad reduction';
  case 13 % DbBAD_GAPFILL
    res = 'bad gapfill';
  case 14 % DbBAD_ALLOC
    res = 'out of memory';
  case 15 % DbBAD_INTERNAL
    res = 'internal error, possible server crash.';
  case 16 % DbBAD_GAP
    res = 'data gap';
  case 17 % DbBAD_ZONE
    res = 'bad zone (between two samples)';
  case 21 % DbBAD_REQUEST
    res = 'bad protocol request';
  case 22 % DbBAD_VALUE
    res = 'bad parameter value';
  case 23 % DbBAD_ACCESS
    res = 'data access denied';
  case 24 % DbBAD_LENGTH
    res = 'protocol length error';
  case 25 % DbNOT_IMPLEMENTED
    res = 'not implemented';
  case 26 % DbBAD_INTERVAL
    res = 'bad interval';
  case 28 % DbBAD_MEMLIMIT
    res = 'memory limit exceeded';
  case 29 % DbBAD_FILE
    res = 'file system error, (missing file or wrong permissions).';
  otherwise
    res = 'unknown error';
end

function res = is_nodata_err(err_id)
res = 0;
% DbBAD_DROP, DbBAD_TIME and DbBAD_EOF are 'no data' errors
if err_id==1 || err_id==2 || err_id==5, res = 1; end
