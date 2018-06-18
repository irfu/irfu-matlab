function [Lm,Blocal,Bmin,J,POSIT] = onera_desp_lib_trace_field_line(kext,options,sysaxes,matlabd,x1,x2,x3,maginput,R0)
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
% function [Lm,Blocal,Bmin,J,POSIT] = onera_desp_lib_trace_field_line(kext,options,sysaxes,matlabd,x1,x2,x3,maginput,R0)
% returns spatial coordinates of field line starting from a single point
% as described in the ONERA/DESP
% library documentation
% Lm, Blocal, Bmin, J are for the starting point
% POSIT is Nbnc x 3
% and the field line x, y, z GEO
% kext - specified the external field model
% For the kext argument, see helps for onera_desp_lib_kext
% options - controls the field tracing
% For the options argument, see helps for onera_desp_lib_options
% sysaxes - sets the coordinate system for the input points
% For the sysaxes argument, see helps for onera_desp_lib_sysaxes
% x1, x2, and x3 is the starting point in the system specified by sysaxes
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
% For default maginput, use []
% if maginput is omitted, [] is assumed
% if R0 is omitted, R0=1 is assumed

if nargin < 8
    maginput = [];
end

if nargin < 9
    R0 = 1;
end

matlabd = datenum(matlabd);

onera_desp_lib_load;

kext = onera_desp_lib_kext(kext);
options = onera_desp_lib_options(options);
sysaxes = onera_desp_lib_sysaxes(sysaxes);
if isempty(maginput)
    maginput = nan(1,25);
end
if size(maginput,2) == 1 % make column vector into row vector
    maginput = maginput';
end

maginput = onera_desp_lib_maginputs(maginput); % NaN to baddata


Nbounce = 20*150; % maximum size of bounce array
Lm = nan;
Blocal = nan(Nbounce,1);
Bmin = Lm;
J = Lm;
[iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);
LmPtr = libpointer('doublePtr',Lm);
BlocalPtr = libpointer('doublePtr',Blocal);
BminPtr = libpointer('doublePtr',Bmin);
JPtr = libpointer('doublePtr',J);
POSIT = nan([3 Nbounce]);
POSITPtr = libpointer('doublePtr',POSIT);
NPOSIT = 0;
NPOSITPtr = libpointer('int32Ptr',NPOSIT);
maginput = maginput';
calllib('onera_desp_lib','trace_field_line2_1_',kext,options,sysaxes,iyear,idoy,UT,x1,x2,x3,maginput,R0,...
    LmPtr,BlocalPtr,BminPtr,JPtr,POSITPtr,NPOSITPtr);
% have to do this next bit because Ptr's aren't really pointers
Lm = get(LmPtr,'value');
Blocal = get(BlocalPtr,'value');
Bmin = get(BminPtr,'value');
J = get(JPtr,'value');
POSIT = get(POSITPtr,'value');
POSIT = reshape(POSIT,3,Nbounce)';
NPOSIT = get(NPOSITPtr,'value');
POSIT = POSIT(1:NPOSIT,:);
Blocal = Blocal(1:NPOSIT);

% the flag value is actually -1d31
Lm(Lm<-1e30) = nan;
Blocal(Blocal<-1e30) = nan;
Bmin(Bmin<-1e30) = nan;
J(J<-1e30) = nan;
POSIT(POSIT<-1e30) = nan;

