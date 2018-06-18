function [Lm,Lstar,Blocal,Bmin,Bmir,J,POSIT,hmin,hmin_lon] = onera_desp_lib_drift_bounce_orbit(kext,options,sysaxes,matlabd,x1,x2,x3,alpha,maginput,R0)
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
%***************************************************************************************************
%
% [Lm,Lstar,Blocal,Bmin,Bmir,J,POSIT,hmin,hmin_lon] = onera_desp_lib_drift_bounce_orbit(kext,options,sysaxes,matlabd,x1,x2,x3,alpha,maginput,R0)
% returns spatial coordinates of drift-bounce orbit starting from a single point
% as described in the ONERA/DESP
% library documentation
% Lm, Lstar, Bmin, J are for the starting point
% Bmir is the mirror field strength
% Blocal is a 25 cell array of Nbnc x 1 local B values
% POSIT is 25 cell arrays of Nbnc x 3
% and the cells contain the field line at each azimuth (up to 1000 points)
% where the 2nd dimension spans x, y, z GEO
% hmin - minimum altitude along drift orbit (not necessarily a point in
% POSIT)
% hmin_lon - longitude of minimum altitude
% kext - specified the external field model
% For the kext argument, see helps for onera_desp_lib_kext
% options - controls the field tracing
% For the options argument, see helps for onera_desp_lib_options
% sysaxes - sets the coordinate system for the input points
% For the sysaxes argument, see helps for onera_desp_lib_sysaxes
% x1, x2, and x3 is the starting point in the system specified by sysaxes
% alpha local pitch angle at starting point, degrees
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
% R0 minmum radial distance allowed for field line traces & bounce orbit
% (R0<1 allowed)
%
% IMPORTANT: all inputs must be present. For those which are not used a dummy value can be provided.
%

if nargin < 9
    maginput = [];
end

if nargin < 10
    R0=1;
end

matlabd = datenum(matlabd);

onera_desp_lib_load;

kext = onera_desp_lib_kext(kext);
in_options = options;
options = onera_desp_lib_options(options);
if isempty(in_options)
    options(1) = 1; % by default, force full drift shell to be traced, otherwise only get one field line
end
sysaxes = onera_desp_lib_sysaxes(sysaxes);
if isempty(maginput)
    maginput = nan(1,25);
end
if (size(maginput,1)==25) && (size(maginput,2)~=25) % 25xN
    maginput = maginput'; % Nx25
end
maginput = onera_desp_lib_maginputs(maginput); % NaN to baddata

Nbounce = 1000; % maximum size of bounce array
Naz = 25; % number of azimuths
Lm = nan;
Lstar = Lm;
Blocalar = nan(Nbounce,Naz);
Bmin = Lm;
Bmir = Lm;
J = Lm;
[iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);
LmPtr = libpointer('doublePtr',Lm);
LstarPtr = libpointer('doublePtr',Lstar);
BlocalPtr = libpointer('doublePtr',Blocalar);
BminPtr = libpointer('doublePtr',Bmin);
BmirPtr = libpointer('doublePtr',Bmir);
JPtr = libpointer('doublePtr',J);
POSITarray = nan([3 Nbounce Naz]);
POSITPtr = libpointer('doublePtr',POSITarray);
NPOSIT = zeros(Naz,1);
NPOSITPtr = libpointer('int32Ptr',NPOSIT);
hminPtr = libpointer('doublePtr',nan);
hmin_lonPtr = libpointer('doublePtr',nan);
maginput = maginput';
calllib('onera_desp_lib','drift_bounce_orbit2_1_',kext,options,sysaxes,iyear,idoy,UT,x1,x2,x3,alpha,maginput,R0,...
    LmPtr,LstarPtr,BlocalPtr,BminPtr,BmirPtr,JPtr,POSITPtr,NPOSITPtr,hminPtr,hmin_lonPtr);
% have to do this next bit because Ptr's aren't really pointers
Lm = get(LmPtr,'value');
Lstar = get(LstarPtr,'value');
Blocalar = get(BlocalPtr,'value');
Bmin = get(BminPtr,'value');
Bmir = get(BmirPtr,'value');
J = get(JPtr,'value');
POSITarray = get(POSITPtr,'value');
POSITarray = reshape(POSITarray,3,Nbounce,Naz);
Blocalar = reshape(Blocalar,Nbounce,Naz);
POSITarray(POSITarray<-1e30) = nan; % flags
Blocalar(Blocalar<-1e30) = nan; % flags
NPOSIT = get(NPOSITPtr,'value');
POSIT = cell(Naz,1);
Blocal = cell(Naz,1);
for i = 1:Naz
    n = NPOSIT(i);
    POSIT{i} = squeeze(POSITarray(:,1:n,i))';
    Blocal{i} = Blocalar(1:n,i);
end
hmin = hminPtr.value;
hmin_lon = hmin_lonPtr.value;

% the flag value is actually -1d31
Lm(Lm<-1e30) = nan;
Lstar(Lstar<-1e30) = nan;
Bmin(Bmin<-1e30) = nan;
Bmir(Bmir<-1e30) = nan;
J(J<-1e30) = nan;
hmin(hmin<-1e30) = nan;
hmin_lon(hmin_lon<-1e30) = nan;

