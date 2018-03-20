function elements = onera_desp_lib_rv2coe(r,v)
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
% function elements = onera_desp_lib_rv2coe(r,v)
% returns classical orbit elements for vehicle 
%   at position r, w/ velocity v
% r - 3x1 ECI position at epoch, km (*** NOTE UNITS km not Re ***)
% v - 3x1 ECI velocity at epoch, km/s 
% elements is a structure with the following elements:
% h - semi-latus rectum (Re)
% a - semi-major axis (Re)
% e - eccentricity
% i - inclination (deg)
% Omega - longitude of ascending node (deg)
% these elements may or may not be defined, depending on the orbit type
%   (e.g., circular, equatorial, etc)
% "5th" element possibilities
% omega - argument of perigee (deg)
% Pi - longitude of perigee (deg)
% "6th" element possibilities, assumes r,v provided at "epoch"
% M0 - mean anomaly at epoch (deg)
% nu0 - true anomaly at epoch (deg)
% u0 - argument of latitude at epoch (deg)
% l0 - true longitude at epoch (deg)

onera_desp_lib_load;
Re =  6371.2; % mean Earth Radius used in GDZ calculation in IRBEM
elements.r = r;
elements.v = v;

h = nan;
hPtr = libpointer('doublePtr',h);

a = nan;
aPtr = libpointer('doublePtr',a);

e = nan;
ePtr = libpointer('doublePtr',e);

i = nan;
iPtr = libpointer('doublePtr',i);

Omega = nan;
OmegaPtr = libpointer('doublePtr',Omega);

omega = nan;
omegaPtr = libpointer('doublePtr',omega);

M0 = nan;
M0Ptr = libpointer('doublePtr',M0);

nu0 = nan;
nu0Ptr = libpointer('doublePtr',nu0);

u0 = nan;
u0Ptr = libpointer('doublePtr',u0);

l0 = nan;
l0Ptr = libpointer('doublePtr',l0);

Pi = nan;
PiPtr = libpointer('doublePtr',Pi);

%       SUBROUTINE rv2coe      ( R, V, P, A, Ecc, Incl, Omega, Argp, Nu,M, ArgLat, TrueLon, LonPer )
calllib('onera_desp_lib','rv2coe_',r,v,hPtr,aPtr,ePtr,iPtr,OmegaPtr,omegaPtr,nu0Ptr,M0Ptr,u0Ptr,l0Ptr,PiPtr);
undef = 999999.1;
elements.h = get(hPtr,'value');
elements.a = get(aPtr,'value');
elements.e = get(ePtr,'value');
elements.i = get(iPtr,'value');
elements.Omega = get(OmegaPtr,'value');
elements.omega = get(omegaPtr,'value');
elements.nu0 = get(nu0Ptr,'value');
elements.M0 = get(M0Ptr,'value');
elements.u0 = get(u0Ptr,'value');
elements.l0 = get(l0Ptr,'value');
elements.Pi = get(PiPtr,'value');
fnames = fieldnames(elements);
for i = 1:length(fnames)
    var = fnames{i};
    if abs(elements.(var)-undef)<abs(undef)/100
        elements = rmfield(elements,var);
    end
end

% convert distances to Re
elements.h = elements.h/Re;
elements.a = elements.a/Re;

% turn the remaining angles into degrees
angles = {'i','Omega','omega','nu0','M0','u0','l0','Pi'};
for i = 1:length(angles)
    var = angles{i};
    if isfield(elements,var)
        elements.(var) = elements.(var)*180/pi;
    end
end
elements.type = 'c';
elements.Re = Re;
