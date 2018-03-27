function Flux = onera_desp_lib_get_ae8_ap8_flux(whichm,energy,BBo,L)
%***************************************************************************************************
% Copyright 2008, T.P. O'Brien
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
% function Flux = onera_desp_lib_get_ae8_ap8_flux(whichm,energy,BBo,L)
% returns AE/AP8 flux from NASA models
% whichm : AE8min, AE8max, AP8min, AP8max or 1, 2, 3, 4 respectively
% energy (in MeV):
% if energy is Nex1, returns differential flux at Ne energies (MeV-1 cm-2 s-1)
% if energy is Nex2, returns differential flux from energy(:,1) to energy(:,2) (MeV-1 cm-2 s-1)
% if energy(:,2) is all infinities, returns integral flux above energy(:,1) (cm-2 s-1)
% matlabd is the time in matlab date format
% L, BBo are McIlwain L and B/B0 in the appropriate magnetic field
% output: size(Flux) = length(x1) x Ne

onera_desp_lib_load;

if isnumeric(whichm)
    iwhichm = whichm;
else
    mname = lower(whichm);
    mname = mname(~ismember(mname,[' ','-']));
    switch(mname)
        case {'ae8min'}
            iwhichm = 1;
        case {'ae8max'}
            iwhichm = 2;
        case {'ap8min'}
            iwhichm = 3;
        case {'ap8max'}
            iwhichm = 4;
        otherwise
            error('Unknown model whichm="%s" in %s',whichm,mfilename);
    end
end

% whichm: long integer to select in which NASA model to fly
% 1=AE8 MIN
% 2=AE8 MAX
% 3=AP8 MIN
% 4=AP8 MAX

NE = size(energy,1); % number of energies
NEmax = 25; % maximum number of energies

ntime = length(BBo);
Nmax = onera_desp_lib_ntime_max; % maximum array size in fortran library
Flux = nan(ntime,NE);
if (ntime>Nmax) || (NE>NEmax)
    % break up the calculation into chunks the libarary can handle
    for i = 1:Nmax:ntime
        ii = i:min(i+Nmax-1,ntime);
        for ie = 1:NEmax:NE
            iie = ie:min(ie+NEmax-1,NE);
            Flux(ii,iie) = onera_desp_lib_get_ae8_ap8_flux(whichm,energy(iie,:),BBo(ii),L(ii));
        end
    end
else

    % whatf: long integer to select what flux to compute
    % 1=differential flux (MeV-1 cm-2 s-1) at energy(1)
    % 2=flux within an ernergy range (MeV-1 cm-2 s-1) - energy(1) to energy(2)
    % 3=integral flux (cm-2 s-1) at energy(1)
    %
    % energy: array(2) of double. If whatf=1 or 3 then energy(2) is not considered.
    %

    if size(energy,2)==1
        whatf = 1; % differential flux
        energy = [energy,nan(NE,1)];
    elseif (size(energy,2)==2) && any(isinf(energy(:,2)))        
        whatf = 3; % integral flux
        if ~all(isinf(energy(:,2)))
            error('%s: if any of second column of "energy" argument is infinity, all must be',mfilename);
        end
    elseif size(energy,2)==2
        whatf = 2; % wide differential flux
    else
        error('%s: "energy" argument of size %d x %d uninterpretable',mfilename,size(energy,1),size(energy,2));
    end

    nanpad = nan(Nmax-ntime,1);
    BBo = [BBo(:);nanpad];
    L = [L(:);nanpad];
    flux = nan(Nmax,NEmax);
    
    FluxPtr = libpointer('doublePtr',flux);
    calllib('onera_desp_lib','get_ae8_ap8_flux_',ntime,iwhichm,whatf,NE,energy',BBo,L,FluxPtr);
    Flux = get(FluxPtr,'value');
end
Flux = Flux(1:ntime,1:NE);
Flux(Flux<-1e30) = nan;
