function Lstar = onera_desp_lib_empiricallstar(kext,options,matlabd,maginput,Lm,J)
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
% Lstar = onera_desp_lib_empiricallstar(kext,options,matlabd,maginput,Lm,J)
% compute fast Lstar
% inputs/outputs have identical meaning to onera_desp_lib_make_lstar
% Works for Olson-Pfitzer Quiet only (ignores kext and field model options)
% ignores noLstar and makePhi option, always gives Lstar
% set maginput to [] if not needed by kext
% currently only works for Olson-Pfitzer Quiet

matlabd = datenum(matlabd);

onera_desp_lib_load;

ntime = length(matlabd);
kext = onera_desp_lib_kext(kext);
options = onera_desp_lib_options(options);
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


siz_in = size(matlabd);

Nmax = onera_desp_lib_ntime_max; % maximum array size in fortran library
Lstar = nan(Nmax,1);
if ntime>Nmax
    % break up the calculation into chunks the libarary can handle
    for i = 1:Nmax:ntime
        ii = i:min(i+Nmax-1,ntime);
        Lstar(ii)= ...
            onera_desp_lib_empiricallstar(kext,options,matlabd(ii),maginput(ii,:),Lm(ii),J(ii));
    end
else
    [iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);
    LstarPtr = libpointer('doublePtr',Lstar);
    maginput = maginput';
    % expand arrays
    iyear = [iyear(:)', nan(1,Nmax-ntime)];
    idoy = [idoy(:)', nan(1,Nmax-ntime)];
    maginput = [maginput, nan(25,Nmax-ntime)];
    Lm = [Lm(:)', nan(1,Nmax-ntime)];
    J = [J(:)', nan(1,Nmax-ntime)];
    
    calllib('onera_desp_lib','empiricallstar1_',ntime,kext,options,iyear,idoy,maginput,Lm,J,LstarPtr);
    % have to do this next bit because Ptr's aren't really pointers
    Lstar = get(LstarPtr,'value');
end

% the flag value is actually -1d31
Lstar(Lstar<-1e30) = nan;

% truncate arrays
Lstar = Lstar(1:ntime);

Lstar = reshape(Lstar,siz_in);

