function pos = onera_desp_lib_sgp4_tle(varargin)
%***************************************************************************************************
% Copyright 2007, T.P. O'Brien
%
% This file is part of IRBEM-LIB.
%
%    IRBEM-LIB is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Lesser General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    IRBEM-LIB is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License
%    along with IRBEM-LIB.  If not, see <http://www.gnu.org/licenses/>.
%
%***************************************************************************************************
%
% pos = onera_desp_lib_sgp4_tle(...)
% computes the position [alt(km),lat(deg),lon(deg)] of satellite orbiting Earth
% into struct: pos.date, pos.alt, pos.lat, pos.lon
% pos = onera_desp_lib_sgp4_tle(InFile)
%   defines start, stop, and delta time to propagate each TLE according to input file
% onera_desp_lib_sgp4_tle(InFile,OutFile)
%   leaves data in the OutFile rather than loading into memory
% pos = onera_desp_lib_sgp4_tle(startsfe,stopsfe,deltasec,InFile)
%   propagates each TLE according to user start and stop time
% startsfe: start time in seconds from date provided in each TLE. This  number can be negative.
% stopsfe: stop time in seconds from date provided in each TLE. This number can be negative.
% deltasec: step time in seconds to propagate TLE.
% InFile: provides the path and name to locate the input TLE to be propagated.
%

if nargin == 1
    % (InFile)
    runtype = 0;
    startsfe = 0;
    stopsfe = 0;
    deltasec = 0;
    InFile = varargin{1};
elseif nargin == 4
    %    startsfe,stopsfe,deltasec,InFile
    runtype = 1;
    startsfe = varargin{1};
    stopsfe = varargin{2};
    deltasec = varargin{3};
    InFile = varargin{4};
elseif nargin == 2
    % (InFile,OutFile)
    runtype = 0;
    startsfe = 0;
    stopsfe = 0;
    deltasec = 0;
    InFile = varargin{1};
    OutFile = varargin{2};
else
    error('%s: Incorrect number of input arguments',mfilename);
end

if ~exist(InFile,'file')
    error('%s: InFile "%s" not found',mfilename,InFile);
end

onera_desp_lib_load;
strlenIn = length(InFile);

if ~exist('OutFile','var')
    delOutFile = 1; % delete the output file
    OutFile = 'onera_desp_lib_sgp4_tle.tmp';
    while exist(OutFile,'file')
        OutFile = sprintf('onera_desp_lib_sgp4_tle.tmp.%.0f',10000*rand(1));
    end
else
    delOutFile = 0; % save the output file
end
strlenOut = length(OutFile);

% void sgp4_tle1_ ( long int * runtype , double * startsfe , double * stopsfe , double * deltasec , char * InFileByte , long int * strlenIn , char * OutFileByte , long int * strlenOut );
calllib('onera_desp_lib','sgp4_tle1_',runtype,startsfe,stopsfe,deltasec,InFile,strlenIn,OutFile,strlenOut);

if delOutFile
    % read the OutFile and then delete it

    % 1: date (dd/mm/yyyy)
    % 2: time (hh:mm:ss)
    % 3: decimal year
    % 4: altitude (km)
    % 5: latitude (deg)
    % 6: longitude (deg)
    % 26/01/2007   0: 4:59.999982 2007.06850263      20096.506220        -0.002629      -124.124210
    % 26/01/2007   0:10: 0.000004 2007.06851218      20096.506204        -0.002508      -122.858292
    [dd,mm,yyyy,hh,minute,ss,decyear,alt,lat,lon] = textread(OutFile,'%d/%d/%d %d:%d:%f %f %f %f %f');

    pos.date = datenum(yyyy,mm,dd,hh,minute,ss);
    pos.alt = alt;
    pos.lat = lat;
    pos.lon = lon;
    pos.sysaxes = 'gdz';

    delete(OutFile); % clean up
else
    pos = OutFile;
end
