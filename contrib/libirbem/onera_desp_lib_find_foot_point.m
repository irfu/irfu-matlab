function [Xfoot,Bfoot,BfootMag] = onera_desp_lib_find_foot_point(kext,options,sysaxes,matlabd,x1,x2,x3,stop_alt,hemi_flag,maginput)
%***************************************************************************************************
% Copyright 2005, T.P. O'Brien
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
% function [Xfoot,Bfoot,BfootMag] = onera_desp_lib_find_foot_point(kext,options,sysaxes,matlabd,x1,x2,x3,stop_alt,hemi_flag,maginput)
% finds the foot point for a field line
% where the foot point is defined as the point at stop_alt km in gdz
% hemi_flag has the following meaning:
% '','same',0 - same Hemisphere as start point
% 'N','North', +1 - Northern Hemisphere
% 'S','South', -1 - Southern Hemisphere
% 'O','Opposite', 2 - opposite Hemisphere as start point
%
% outputs:
%   Xfoot (gdz), Bfoot (GEO, nT), and BfootMag (nT)
% kext - specified the external field model
% For the kext argument, see helps for onera_desp_lib_kext
% options - controls the field tracing
% For the options argument, see helps for onera_desp_lib_options
% sysaxes - sets the coordinate system for the input points
% For the sysaxes argument, see helps for onera_desp_lib_sysaxes
% x1, x2, and x3 are the points of interest in the system specified by sysaxes
% alpha is the local pitch angle in degrees
% maginput - [length(x1) x 25] provides inputs to dynamic external field models
% (if maginput is omitted or empty, then a matrix of zeros is assumed)
% maginput(1st element,*) =Kp: value of Kp as in OMNI2 files but has to be double instead of integer type
% maginput(2nd element,*) =Dst: Dst index (nT)
% maginput(3rd element,*) =dens: Solar Wind density (cm-3)
% maginput(4th element,*) =velo: Solar Wind velocity (km/s)
% maginput(5th element,*) =Pdyn: Solar Wind dynamic pressure (nPa)
% maginput(6th element,*) =ByIMF: GSM y component of IMF mag. field (nT)
% maginput(7th element,*) =BzIMF: GSM z component of IMF mag. field (nT)
% maginput(8th element,*) =G1:  G1=< Vsw*(Bperp/40)^2/(1+Bperp/40)*sin^3(theta/2) > where the <> mean an average over the previous 1 hour, Vsw is the solar wind speed, Bperp is the transverse IMF component (GSM) and theta it's clock angle.
% maginput(9th element,*) =G2: G2=< a*Vsw*Bs > where the <> mean an average over the previous 1 hour, Vsw is the solar wind speed, Bs=|IMF Bz| when IMF Bz < 0 and Bs=0 when IMF Bz > 0, a=0.005.
% maginput(10th element,*) =G3:  G3=< Vsw*Dsw*Bs /2000.>
% where the <> mean an average over the previous 1 hour, Vsw is the solar wind speed, Dsw is the solar wind density, Bs=|IMF Bz| when IMF Bz < 0 and Bs=0 when IMF Bz > 0.
% maginput(11th element,*) =W1 see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
% maginput(12th element,*) =W2 see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
% maginput(13th element,*) =W3 see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
% maginput(14th element,*) =W4 see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
% maginput(15th element,*) =W5 see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
% maginput(16th element,*) =W6 see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
% maginput(17th element,*) =AL the auroral index
%
% maginput(18th element,*) to maginput(25th element,*): for future use
%
% IMPORTANT: all inputs must be present. For those which are not used a dummy value can be provided.
%

if nargin < 10
    maginput = [];
end

if ischar(hemi_flag)
    switch(upper(hemi_flag))
        case {'','SAME'}, hemi_flag = 0;
        case {'N','NORTH'}, hemi_flag = +1;
        case {'S','SOUTH'}, hemi_flag = -1;
        case {'O','OPPOSITE'}, hemi_flag = 2;
        otherwise
            error('Unknown hemi_flag "%s"',hemi_flag);
    end
end


matlabd = datenum(matlabd);

onera_desp_lib_load;

ntime = length(x1);
kext = onera_desp_lib_kext(kext);
options = onera_desp_lib_options(options);
sysaxes = onera_desp_lib_sysaxes(sysaxes);
if isempty(maginput)
    maginput = nan(ntime,25);
end
if (size(maginput,1)==25) && (size(maginput,2)~=25) % 25xN
    maginput = maginput'; % Nx25
end
if size(maginput,1) ~= ntime
    maginput = repmat(maginput,ntime,1);
end
if length(matlabd)==1
    matlabd = repmat(matlabd,ntime,1);
end

maginput = onera_desp_lib_maginputs(maginput); % NaN to baddata


[iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);
Xfoot = nan(ntime,3);
Bfoot = nan(ntime,3);
BfootMag = nan(ntime,1);
XfootPtr = libpointer('doublePtr',nan(3,1));
BfootPtr = libpointer('doublePtr',nan(3,1));
BfootMagPtr = libpointer('doublePtr',nan);
for i = 1:ntime
    calllib('onera_desp_lib','find_foot_point1_',kext,options,sysaxes,iyear(i),idoy(i),UT(i),x1(i),x2(i),x3(i),stop_alt,hemi_flag,maginput(i,:),...
        XfootPtr,BfootPtr,BfootMagPtr);
    % have to do this next bit because Ptr's aren't really pointers
    Xfoot(i,:) = get(XfootPtr,'value');
    Bfoot(i,:) = get(BfootPtr,'value');
    BfootMag(i) = get(BfootMagPtr,'value');
end

% the flag value is actually -1d31
% BfootMag<=0 also indicates error
Xfoot(Xfoot<-1e30) = nan;
Bfoot(Bfoot<-1e30) = nan;
Xfoot(BfootMag<=0,:) = nan; % discard debug info
Bfoot(BfootMag<=0,:) = nan; % discard debug info
BfootMag(BfootMag<=0) = nan;
