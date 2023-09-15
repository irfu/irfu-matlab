function [Lower_flux,Mean_flux,Upper_flux] = onera_desp_lib_fly_in_meo_gnss(launch_year,duration,whichm,energy)
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
% [Lower_flux,Mean_flux,Upper_flux] = onera_desp_lib_fly_in_meo_gnss(launch_year,duration,whichm,energy)
% This function allows one to fly any MEO GNSS type spacecraft in MEO ONERA
% models. The use of the model is limited to GPS altitude (~20000 km - 55°
% inclination).
%
% launch_year: year of launch (uses integer part only)
% mission_duration: duration of mission in years (uses integer part only)
% whichm : 1 = 'MEO-V1',2 = 'MEO-V2'
% energy (in MeV):
% if energy is Nex1, returns differential flux at Ne energies (MeV-1 cm-2 s-1 sr-1)
% if energy is Nex2, returns differential flux from energy(:,1) to energy(:,2) (MeV-1 cm-2 s-1 sr-1)
% if energy(:,2) is all infinities, returns integral flux above energy(:,1) (cm-2 s-1 sr-1)
% output: all size Ne x 1
% Lower_flux: Lower flux for all energies
% averaged over entire mission duration - This has to be considered as a
% lower envelop to bound expected flux at MEO-GNSS for any solar cycle
% Mean_flux: Mean flux for all energies
% averaged over entire mission duration - This spectrum is an averaged
% expected flux at MEO-GNSS for any solar cycle, no margins are included
% at this point.
% Upper_flux: Upper flux for all energies
% averaged over entire mission duration - This has to be considered as an
% upper envelop for expected flux at MEO-GNSS for any solar cycle, this
% spectrum can be used for any conservative approach as margins are
% included at this point. Note that the margins are energy dependent and can
% be assesed by looking at Upper_flux/Mean_flux
% Note: all flux are expressed in MeV-1 cm-2 s-1 sr-1 for differential flux
% or in MeV-1 cm-2 s-1 sr-1 for integrated flux. To derive omnidirectional
% flux at MEO-GNSS one should multiply these flux values by 4pi.

% note: this is just a wrapper for fly_in_ige1, which is called with
% 100+iwhichm to indicated meo_gnss rather than ige

onera_desp_lib_load;

if isnumeric(whichm)
    iwhichm = whichm;
else
    mname = lower(whichm); % lower case
    mname = mname(~ismember(mname,[' ','-'])); % remove spaces and dashes
    switch(mname)
        case {'v1','meov1','meo1'}
            iwhichm = 1;
        case {'v2','meov2','meo2'}
            iwhichm = 2;
        otherwise
            error('Unknown model whichm="%s" in %s',whichm,mfilename);
    end
end

[Lower_flux,Mean_flux,Upper_flux] = onera_desp_lib_fly_in_ige(launch_year,duration,100+iwhichm,energy);
