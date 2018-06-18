function [Lm,Lstar,Bmirror,Bmin,J,MLT] = onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd,x1,x2,x3,varargin)
%***************************************************************************************************
% Copyright 2006-2009, T.P. O'Brien
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
% [Lm,Lstar,Bmirror,Bmin,J,MLT] = onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd,x1,x2,x3,maginput)
% wrapper for variants of make_lstar for locally mirroring particles
% [Lm,Lstar,Bmirror,Bmin,J,MLT] = onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd,x1,x2,x3,alpha,maginput)
% wrapper for variants of make_lstar_shell_splitting for particles w/ local
% pitch angle alpha, degrees
% func_name is a string that identifies which DLL function to call.
% other inputs/outputs identical to onera_desp_lib_make_lstar

switch(lower(func_name))
    case 'onera_desp_lib_make_lstar'
        libfunc_name = 'make_lstar1_';
        splitting = false;
    case 'onera_desp_lib_landi2lstar'
        libfunc_name = 'landi2lstar1_';
        splitting = false;
    case 'onera_desp_lib_make_lstar_shell_splitting'
        libfunc_name = 'make_lstar_shell_splitting1_';
        splitting = true;
    case 'onera_desp_lib_landi2lstar_shell_splitting'
        libfunc_name = 'landi2lstar_shell_splitting1_';
        splitting = true;
    otherwise
        error('Unknown func_name %s',func_name);
end

if splitting
    alpha = varargin{1};
    imaginput = 2;
    Nmaxpa = 25; % maximum number of pitch angles for splitting functions
else
    alpha = 90;
    imaginput = 1;
    Nmaxpa = 1; % maximum number of pitch angles is 1 for non-splitting case
end

if length(varargin)>=imaginput
    maginput = varargin{imaginput};
else
    maginput = [];
end

matlabd = datenum(matlabd);

onera_desp_lib_load;

ntime = length(x1);
nipa = length(alpha);
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

Lm = nan(ntime,nipa);
Lstar = nan(ntime,nipa);
Bmirror = nan(ntime,nipa);
Bmin = nan(ntime,1);
J = nan(ntime,nipa);
MLT = nan(ntime,1);
Nmax = onera_desp_lib_ntime_max; % maximum array size in fortran library
if ntime>Nmax
    % break up the calculation into chunks the libarary can handle
    for i = 1:Nmax:ntime
        ii = i:min(i+Nmax-1,ntime);
        if splitting
            [Lm(ii,:),Lstar(ii,:),Bmirror(ii,:),Bmin(ii),J(ii,:),MLT(ii)] = ...
                onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd(ii),x1(ii),x2(ii),x3(ii),alpha,maginput(ii,:));
        else
            [Lm(ii,:),Lstar(ii,:),Bmirror(ii,:),Bmin(ii),J(ii,:),MLT(ii)] = ...
                onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd(ii),x1(ii),x2(ii),x3(ii),maginput(ii,:));
        end
    end
elseif nipa>Nmaxpa
    % break up the calculation into chunks the libarary can handle
    for i = 1:Nmaxpa:nipa
        ii = i:min(i+Nmaxpa-1,nipa);
        if splitting
            [Lm(:,ii),Lstar(:,ii),Bmirror(:,ii),Bmin_tmp,J(:,ii),MLT_tmp] = ...
                onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd,x1,x2,x3,alpha(ii),maginput);
        else
            [Lm(:,ii),Lstar(:,ii),Bmirror(:,ii),Bmin_tmp,J(:,ii),MLT_tmp] = ...
                onera_desp_lib_make_lstar_core(func_name,kext,options,sysaxes,matlabd,x1,x2,x3,maginput);
        end
        nanMLT = ~isfinite(MLT);
        MLT(nanMLT) = MLT_tmp(nanMLT);
        nanBmin = ~isfinite(Bmin);
        Bmin(nanBmin) = Bmin_tmp(nanBmin);
    end
else
    % reinitialize Lm, etc for full size
    Lm = nan(Nmax,Nmaxpa);
    Lstar = nan(Nmax,Nmaxpa);
    Bmirror = nan(Nmax,Nmaxpa);
    Bmin = nan(Nmax,1);
    J = nan(Nmax,Nmaxpa);
    MLT = nan(Nmax,1);
    
    [iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);
    LmPtr = libpointer('doublePtr',Lm);
    LstarPtr = libpointer('doublePtr',Lstar);
    BmirrorPtr = libpointer('doublePtr',Bmirror);
    BminPtr = libpointer('doublePtr',Bmin);
    JPtr = libpointer('doublePtr',J);
    MLTPtr = libpointer('doublePtr',MLT);
    if nipa<Nmaxpa
        alpha = [alpha(:)',nan(1,Nmaxpa-nipa)]; % pad alpha
    end
    maginput = maginput';
    % expand arrays
    iyear = [iyear(:)', nan(1,Nmax-ntime)];
    idoy = [idoy(:)', nan(1,Nmax-ntime)];
    UT = [UT(:)', nan(1,Nmax-ntime)];
    x1 = [x1(:)', nan(1,Nmax-ntime)];
    x2 = [x2(:)', nan(1,Nmax-ntime)];
    x3 = [x3(:)', nan(1,Nmax-ntime)];
    maginput = [maginput, nan(25,Nmax-ntime)];
    
    if splitting
        calllib('onera_desp_lib',libfunc_name,ntime,nipa,kext,options,sysaxes,iyear,idoy,UT,x1,x2,x3,alpha,maginput,...
            LmPtr,LstarPtr,BmirrorPtr,BminPtr,JPtr,MLTPtr);
    else
        calllib('onera_desp_lib',libfunc_name,ntime,kext,options,sysaxes,iyear,idoy,UT,x1,x2,x3,maginput,...
            LmPtr,LstarPtr,BmirrorPtr,BminPtr,JPtr,MLTPtr);
    end
    
    % have to do this next bit because Ptr's aren't really pointers
    Lm = get(LmPtr,'value');
    Lstar = get(LstarPtr,'value');
    Bmirror = get(BmirrorPtr,'value');
    Bmin = get(BminPtr,'value');
    J = get(JPtr,'value');
    MLT = get(MLTPtr,'value');
    
    % shrinkwrap Lm, etc
    Lm = Lm(1:ntime,1:nipa);
    Lstar = Lstar(1:ntime,1:nipa);
    Bmirror = Bmirror(1:ntime,1:nipa);
    Bmin = Bmin(1:ntime);
    J = J(1:ntime,1:nipa);
    MLT = MLT(1:ntime);
end

% flags to nan
% the flag value is actually -1d31
Lm(Lm<-1e30) = nan;
Lstar(Lstar<-1e30) = nan;
Bmirror(Bmirror<-1e30) = nan;
Bmin(Bmin<-1e30) = nan;
J(J<-1e30) = nan;
MLT(MLT<-1e30) = nan;
