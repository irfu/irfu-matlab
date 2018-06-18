function sysaxes = onera_desp_lib_sysaxes(sysaxes)
%***************************************************************************************************
% Copyright 2006, T.P. O'Brien
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
% function sysaxes = onera_desp_lib_sysaxes(sysaxes)
% converts sysaxes string to number, leaves number unchanged
% sysaxes: long integer to define which coordinate system is provided in 
% 
% 
% 0: GDZ (alti, lati, longi - km,deg.,deg) 
% 1: GEO (cartesian) - Re 
% 2: GSM (cartesian) - Re 
% 3: GSE (cartesian) - Re 
% 4: SM (cartesian) - Re 
% 5: GEI (cartesian) - Re (ECI)
% 6: MAG (cartesian) - Re 
% 7: SPH (geo in spherical) - (radial distance, lati, longi - Re, deg., deg.)
% 8: RLL  (radial distance, lati, longi - Re, deg., deg. - prefered to 7) 
% 9: HEE
% 10: HAE
% 11: HEEQ
if isempty(sysaxes)
    sysaxes = 0;
    return
end
if isnumeric(sysaxes)
    return
end

switch(upper(sysaxes))
    case {'GDZ'},  sysaxes = 0;
    case {'GEO'},  sysaxes = 1;
    case {'GSM'},  sysaxes = 2;
    case {'GSE'},  sysaxes = 3;
    case {'SM'},   sysaxes = 4;
    case {'GEI','ECI'},  sysaxes = 5;
    case {'MAG'},  sysaxes = 6;
    case {'SPH'},  sysaxes = 7;
    case {'RLL'},  sysaxes = 8;
    case {'HEE'},  sysaxes = 9;
    case {'HAE'},  sysaxes = 10;
    case {'HEEQ'},  sysaxes = 11;
    otherwise
        error('Unknown value for sysaxes "%s" in %s',sysaxes,mfilename);
end
