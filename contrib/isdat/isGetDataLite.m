function [timedat,reqdat,iserr]=isGetDataLite(db,starttime,duration,proj,mem,inst,sig,sen,chan,param,units)
%isGetDataLite  Get science data from ISDAT server
%
% [tim,dat]=isGetDataLite(DB,st,dur,pro,mem,ins,sig,sen,cha,par,[units])
% [tim,dat,iserr]=isGetDataLite(DB,st,dur,pro,mem,ins,sig,sen,cha,par,[units])
%
% Retrieves data from an isdat server.
%
% This assumes that one has previously opened a filehandle
% DB to the database connection using the Mat_DbOpen function.
%
% The server attempts to retrieve data according to the requested
% start time st = [year mon day hour min sec.msec] (vector of doubles) or
% st = isdat_epoch (double, see below),
% duration dur in units of seconds,
% the dataset: the project pro, the member mem, the instrument ins,
% the signal sig, the sensor sen, the channel cha, and the parameter par.
%
% units - 'phys' (default) or 'tm'
%
% This specification list need not be given in full if the
% dataset is specified with fewer levels.
% If data according to the request is available, it is returned
% in doubles array dat for which the sample index is in matlab dim 1.
%
% The timing of the individual samples of data is returned in tim as a double
% representing seconds since 1-Jan-1970. To convert epoch time to year, minute,
% second format use function FROMEPOCH.
%
%  Example ( assuming dbh running locally as "dbh :0 &" )
%
%  >> DB = Mat_DbOpen('unix:0');
%  >> [tim, dat] = isGetDataLite(DB, [1993 08 01 15 10 00], 2, ...
%       'Test', '*', 'wave', 'sine', 'wf');
%
% should return a sinusoid from the "Test" dataset (which always exist in the
% database for the times given here).
%
% See also: Mat_DbOpen, FROMEPOCH.
%

narginchk(4,11)

if nargin < 11, units = 'phys';
else
  if ~( strcmpi(units,'tm') || strcmpi(units,'phys') )
    error('Invalid value for UNITS')
  end
end

request = zeros(17,1);
if length(starttime) > 1, request(1) = toepoch(starttime);
else, request(1) = starttime; end
request(2) = 0;
request(3) = duration;
request(4) = 0;

if nargin<10
  param=' ';
  if nargin < 9
    chan = ' ';
    if nargin < 8
      sen = ' ';
      if nargin < 7
        sig = ' ';
        if nargin < 6
          inst = ' ';
          if nargin < 5
            mem = ' ';
          end
        end
      end
    end
  end
end
[dataSetId, err]=Mat_DbName2Spec(db,proj,mem,inst,sig,sen,chan,param);

request(5)  = dataSetId(1);
request(6)  = dataSetId(2);
request(7)  = dataSetId(3);
request(8)  = dataSetId(4);
request(9)  = dataSetId(5);
request(10) = dataSetId(6);
request(11) = dataSetId(7);

if strcmpi(units,'tm')
  request(12) = 1;	%data_request_units       = DbUN_TM
else
  request(12) = 3;	%data_request_units       = DbUN_PHYS
end
request(13) = 1;	%data_request_reduction   = DbRED_NONE
request(14) = 0;	%data_request_samples	  = DbUNDEF
request(15) = 1;	%data_request_gapFill	  = DbGAP_NAN
request(16) = 1;	%data_request_pack        = DbPACK_TIMETAG
request(17) = -1;	%data_request_dataVersion = -1

lasterr('')
try
  [reqdat, result, info, timedat] = Mat_DbGetDataLite(db,request);
catch
  error('ISDAT','Error getting data from ISDAT')
end

if nargout>2, iserr = result; end
