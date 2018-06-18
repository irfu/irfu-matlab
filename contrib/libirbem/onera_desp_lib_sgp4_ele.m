function pos = onera_desp_lib_sgp4_ele(elements,startdate,enddate,deltasec,sysAxesOUT)
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
% function pos = onera_desp_lib_sgp4_ele(elements,startdate,enddate,deltasec,sysAxesOUT)
% computes the position of satellite orbiting Earth
% into struct: pos.date (Nx1), pos.X (Nx3)
% also returns pos.elements and pos.sysaxes for record keeping
% elements is a structure containing orbital elements
% startdate: start time in matlab date format
% enddate: end time in matlab date format
% deltasec: step time in ***seconds*** to produce orbit outputs
% sysAxesOUT: coordinate system for output pos.X (uses convention of onera_desp_lib_sysaxes)
% all elemety types rely on elements.epoch (defaults to 0)
% The following element types are supported:
% elements.type = 'onera' or 'o', requires i, A_p, A_a, Omega, omega, M0
% elements.type = 'classical' or 'c', requires a, e, i, Omega, omega, T
%   (omega may be replaced by Pi)
%   (T may be replaced by nu0, u0, l0, M0)
% elements.type = 'mean' or 'm', requires n, e, i, Omega, omega, M0
% elements.type = 'rv', requires r,v (km and km/s, ECI)
% elements.type = 'solar' or 's', requires i, A_p, A_a, H_a, H_i, T
% elements.type = 'geo' or 'g', requires geo_lon
%
% definitions:
% epoch - time for which elements are given (default startdate)
% a - semi-major axis (Re)
% e - eccentricity
% i - inclination (deg)
% Omega - longitude of ascending node (deg)
% omega - argument of perigee (deg)
% Pi - longitude of perigee (deg)
% T - time of perigee passage (matlab date number)
% n - mean motion (rev/day)
% r - ECI position at epoch, km
% v - ECI velocity at epoch, km/s
% A_p - altitude of perigee (km, geocentric)
% A_a - altitude of apogee (km, geocentric)
% H_a - local time of apogee (hrs)
% H_i - local time of maximum inclination (hrs)
% nu0 - true anomaly at epoch (deg)
% u0 - argument of latitude at epoch (deg)
% l0 - true longitude at epoch (deg)
% M0 - mean anomaly at epoch (deg)
% geo_lon - longitude in GEO coordinates (deg)

if nargin < 5
    sysAxesOUT = 'gdz';
end

if ~isfield(elements,'epoch')
    elements.epoch = startdate;
end

pos.elements = elements;
pos.sysaxes = sysAxesOUT;

sysAxesOUT = onera_desp_lib_sysaxes(sysAxesOUT);

ele_opts = zeros(5,1);
rv = 0;
switch(lower(elements.type(1)))
    case 'g' % geosynchronous orbit, do manually
        geo_alt = 35786; % wikipedia value
        sys_in = onera_desp_lib_sysaxes('gdz');
        pos.date = (startdate:(deltasec/24/60/60):enddate)';
        pos.X = onera_desp_lib_coord_trans(repmat([geo_alt 0 elements.geo_lon],length(pos.date),1),[sys_in sysAxesOUT],pos.date);
        return
    case 'o'
        ele_opts(1) = 1;
        e_names = {'i','A_p','A_a','Omega','omega','M0'};
    case 'c'
        ele_opts(1) = 2;
        e_names = {'a','e','i','Omega','omega','T'};
    case 'r'
        ele_opts(1) = 3;
        rv = 1;
    case 's'
        ele_opts(1) = 4;
        e_names = {'i','A_p','A_a','H_a','H_i','T'};
    case 'm'
        ele_opts(1) = 5;
        e_names = {'n','e','i','Omega','omega','M0'};
end

eles = nan(6,1);

if rv
    if ~isfield(elements,'r')
        error('%s: element "r" not found', mfilename);
    elseif ~isfield(elements,'v')
        error('%s: element "v" not found', mfilename);
    else
        eles(1:3) = elements.r;
        eles(4:6) = elements.v;
    end
else   
    % handle alternatives for e5
    if isfield(elements,'omega')
        e_names{5} = 'omega';
        ele_opts(2) = 1;
    elseif isfield(elements,'Pi')
        e_names{5} = 'Pi';
        ele_opts(2) = 2;
    end

    % handle alternatives for e6
    if isfield(elements,'T')
        elements.Tsfe = (elements.T-elements.epoch)*24*60*60; 
        e_names{6} = 'Tsfe';
        ele_opts(3) = 1;
    elseif isfield(elements,'nu0')
        e_names{6} = 'nu0';
        ele_opts(3) = 2;
    elseif isfield(elements,'u0')
        e_names{6} = 'u0';
        ele_opts(3) = 3;
    elseif isfield(elements,'l0')
        e_names{6} = 'l0';
        ele_opts(3) = 4;
    elseif isfield(elements,'M0')
        e_names{6} = 'M0';
        ele_opts(3) = 5;
    end

    for i = 1:6
        if isfield(elements,e_names{i})
            eles(i) = elements.(e_names{i});
        else
            error('%s: element "%s" not found', mfilename, e_names{i});
        end
    end

end

onera_desp_lib_load;
[Yr,Mon,Day,Hr,Minute,Sec] = datevec(elements.epoch);

startsfe = (startdate-elements.epoch)*24*60*60; % seconds from epoch
stopsfe = (enddate-elements.epoch)*24*60*60; % seconds from epoch
nmax = onera_desp_lib_ntime_max;
ntimes = floor((stopsfe-startsfe)/deltasec);

pos.date = nan(ntimes,1);
pos.X = nan(ntimes,3);

i0 = 1;
while i0 <= ntimes
    ii = i0:(min(ntimes,i0+nmax-1));
    ii0 = ii-i0+1;
    vstartsfe = startsfe + (ii(1)-1)*deltasec;
    vstopsfe = startsfe + (ii(end))*deltasec;
    iYr = nan(nmax,1);
    iYrPtr = libpointer('int32Ptr',iYr);
    iDoy = nan(nmax,1);
    iDoyPtr = libpointer('int32Ptr',iDoy);
    UT = nan(nmax,1);
    UTptr = libpointer('doublePtr',UT);
    x1 = nan(nmax,1);
    x1Ptr = libpointer('doublePtr',x1);
    x2 = nan(nmax,1);
    x2Ptr = libpointer('doublePtr',x2);
    x3 = nan(nmax,1);
    x3Ptr = libpointer('doublePtr',x3);
    calllib('onera_desp_lib','sgp4_ele1_',sysAxesOUT,Yr,Mon,Day,Hr,Minute,Sec,...
        eles(1),eles(2),eles(3),eles(4),eles(5),eles(6),ele_opts,...
        vstartsfe,vstopsfe,deltasec,iYrPtr,iDoyPtr,UTptr,x1Ptr,x2Ptr,x3Ptr);
    iYear = double(get(iYrPtr,'value'));
    ikeep = iYear(1:length(ii))>0; % trim out entries not reached by sgp4_ele1
    ii = ii(ikeep); 
    if isempty(ii)
        break;
    end
    ii0 = ii0(ikeep);
    date = datenum(iYear,1,double(get(iDoyPtr,'value')),0,0,double(get(UTptr,'value')));
    ii0 = ii0(ikeep);
    pos.date(ii) = date(ii0);
    x1 = get(x1Ptr,'value');
    pos.X(ii,1) = x1(ii0);
    x2 = get(x2Ptr,'value');
    pos.X(ii,2) = x2(ii0);
    x3 = get(x3Ptr,'value');
    pos.X(ii,3) = x3(ii0);
    i0 = ii(end)+1;
end

ikeep = isfinite(pos.date); % trim out entries not reached by sgp4_ele1
pos.date = pos.date(ikeep);
pos.X = pos.X(ikeep,:);

