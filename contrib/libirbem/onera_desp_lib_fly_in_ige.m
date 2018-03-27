function [Lower_flux,Mean_flux,Upper_flux] = onera_desp_lib_fly_in_ige(launch_year,duration,whichm,energy)
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
% [Lower_flux,Mean_flux,Upper_flux] = onera_desp_lib_fly_in_ige(launch_year,duration,whichm,energy)
% returns geosynchronous electron flux from POLE/IGE models
% launch_year: year of launch (uses integer part only)
% mission_duration: duration of mission in years (uses integer part only)
% whichm : 1 = 'POLE-V1',2 = 'POLE-V2',3 = 'IGE-2006'
% energy (in MeV):
% if energy is Nex1, returns differential flux at Ne energies (MeV-1 cm-2 s-1 sr-1)
% if energy is Nex2, returns differential flux from energy(:,1) to energy(:,2) (MeV-1 cm-2 s-1 sr-1)
% if energy(:,2) is all infinities, returns integral flux above energy(:,1) (cm-2 s-1 sr-1)
% output: all size Ne x 1
% Lower_flux: Lower flux for all energies averaged over entire mission
%   duration - This has to be considered as a lower envelop to bound
%   expected flux at GEO for any solar cycle 
% Mean_flux: Mean flux for all energies averaged over entire mission
%   duration - This spectrum is an averaged expected flux at GEO for any
%   solar cycle
% Upper_flux: Upper flux for all energies averaged over entire mission
%   duration - This has to be considered as an upper envelop for expected 
%   flux at GEO for any solar cycle

% note: this is also a wrapper for fly_in_meo_gnss for whichm>100
% in which case whichm is converted to whichm-100 and fly_in_meo_gnss1_ is
% called instead of fly_in_ige1_

onera_desp_lib_load;

if isnumeric(whichm)
    iwhichm = whichm;
else
    mname = lower(whichm); % lower case
    mname = mname(~ismember(mname,[' ','-'])); % remove spaces and dashes
    switch(mname)
        case {'polev1'}
            iwhichm = 1;
        case {'polev2'}
            iwhichm = 2;
        case {'ige2006'}
            iwhichm = 3;
        otherwise
            error('Unknown model whichm="%s" in %s',whichm,mfilename);
    end
end

iwhichm0 = iwhichm; % store this, as we're about to change it for gnss case

if iwhichm >= 100
    iwhichm = iwhichm-100;
    libfuncname = 'fly_in_meo_gnss1_';
else
    libfuncname = 'fly_in_ige1_';
end

NE = size(energy,1); % number of energies
NEmax = 50; % maximum number of energies

if NE>NEmax
    % break up the calculation into chunks the libarary can handle
    Lower_flux = nan(NE,1);
    Mean_flux = nan(NE,1);
    Upper_flux = nan(NE,1);
    for ie = 1:NEmax:NE
        iie = ie:min(ie+NEmax-1,NE);
        [lf,mf,uf] = onera_desp_lib_fly_in_ige(launch_year,duration,iwhichm0,energy(iie,:));
        Lower_flux(iie) = lf;
        Mean_flux(iie) = mf;
        Upper_flux(iie) = uf;
    end
else
    
    launch_year = floor(launch_year);
    duration = floor(duration);

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
        energy(:,2) = 0; % zeros instead of inf/nan
    elseif size(energy,2)==2
        whatf = 2; % wide differential flux
    else
        error('%s: "energy" argument of size %d x %d uninterpretable',mfilename,size(energy,1),size(energy,2));
    end
    
    energy = [energy; zeros(NEmax-NE,2)]; % zeros instead of inf/nan

    LFluxPtr = libpointer('doublePtr',zeros(NEmax,1));
    MFluxPtr = libpointer('doublePtr',zeros(NEmax,1));
    UFluxPtr = libpointer('doublePtr',zeros(NEmax,1));
    calllib('onera_desp_lib',libfuncname,launch_year,duration,iwhichm,whatf,NE,energy',LFluxPtr,MFluxPtr,UFluxPtr);
    Lower_flux = get(LFluxPtr,'value');
    Mean_flux = get(MFluxPtr,'value');
    Upper_flux = get(UFluxPtr,'value');
    Lower_flux = Lower_flux(1:NE);
    Mean_flux = Mean_flux(1:NE);
    Upper_flux = Upper_flux(1:NE);
    Lower_flux(Lower_flux<0) = nan; % flag is -1e30
    Mean_flux(Mean_flux<0) = nan; % flag is -1e30
    Upper_flux(Upper_flux<0) = nan; % flag is -1e30
end
