function [Bmin,xGEO] = onera_desp_lib_find_magequator(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)
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
% function [Bmin,xGEO] = onera_desp_lib_find_magequator(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)
% finds the magnetic equator for a particle
% xGEO has size length(x1) x 3, Bmin has size length(x1) x 1
% kext - specified the external field model
% For the kext argument, see helps for onera_desp_lib_kext
% options - controls the field tracing
% For the options argument, see helps for onera_desp_lib_options
% sysaxes - sets the coordinate system for the input points
% For the sysaxes argument, see helps for onera_desp_lib_sysaxes
% x1, x2, and x3 are the points of interest in the system specified by sysaxes
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

if nargin < 8
    maginput = [];
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
Bmin = nan(ntime,1);
xGEO = nan([ntime,3]);
bmin = nan;
xgeo = nan(1,3);
bminPtr = libpointer('doublePtr',bmin);
xgeoPtr = libpointer('doublePtr',xgeo);
for i = 1:ntime
    calllib('onera_desp_lib','find_magequator1_',kext,options,sysaxes,iyear(i),idoy(i),UT(i),x1(i),x2(i),x3(i),maginput(i,:),...
        bminPtr,xgeoPtr);
    % have to do this next bit because Ptr's aren't really pointers
    Bmin(i) = get(bminPtr,'value');
    xGEO(i,:) = get(xgeoPtr,'value');
end

% the flag value is actually -1d31
Bmin(Bmin<-1e30) = nan;
xGEO(xGEO<-1e30) = nan;
