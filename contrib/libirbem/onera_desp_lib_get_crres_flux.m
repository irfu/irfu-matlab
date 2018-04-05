function Flux = onera_desp_lib_get_crres_flux(whichm,energy,BBo,L,Ap15,crres_path)
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
% function Flux = onera_desp_lib_get_crres_flux(whichm,energy,BBo,L,Ap15,crres_path)
% returns CRRESELE/CRRESPRO flux from AFRL models
% whichm :
% 1 or 'CRRESPRO QUIET'
% 2 or 'CRRESPRO ACTIVE'
% 3 or 'CRRESELE AVERAGE'
% 4 or 'CRRESELE WORST CASE'
% 5 or 'CRRESELE Ap15'
% if energy is Nex1, returns differential flux at N energies (MeV-1 cm-2 s-1)
% if energy is Nex2, returns differential flux from energy(:,1) to energy(:,2) (MeV-1 cm-2 s-1)
% if energy(:,2) is all infinities, returns integral flux above energy(:,1) (cm-2 s-1)
% L, BBo are McIlwain L and B/B0 in the appropriate magnetic field
% crres_path is path to AFRL CRRES text files
% Ap15 : array (same size as x1) that provides the 15-day average of Ap
% prior to the time of interest (only needed by the CRRESELE Ap15 model)
% if not supplied, it will be determined by searching the matlab path
% for a file named crrespro_active.txt and will supply that path
% supply [] for Ap15 if not using Ap15 model
% output: size(Flux) = length(x1) x Ne


% MODULE #        DESCRIPTION
%   0             CRRES_Ap (5.0 - 7.5)
%   1             CRRES_Ap (7.5 - 10)
%   2             CRRES_Ap (10 - 15)
%   3             CRRES_Ap (15 - 20)
%   4             CRRES_Ap (20 - 25)
%   5             CRRES_Ap (25 - 55)

if nargin < 5
    Ap15 = [];
end

if nargin < 6
    crres_path = '';
end

onera_desp_lib_load;

if isnumeric(whichm)
    iwhichm = whichm;
else
    mname = lower(whichm);
    switch(mname)
        case {'crrespro quiet'}
            iwhichm = 1;
        case {'crrespro active'}
            iwhichm = 2;
        case {'crresele average'}
            iwhichm = 3;
        case {'crresele worst case'}
            iwhichm = 4;
        case {'crresele ap15'}
            iwhichm = 5;
        otherwise
            error('Unknown model whichm="%s" in %s',whichm,mfilename);
    end
end

%
% whichm: long integer to select in which AFRL CRRES model to fly
% 1=CRRESPRO QUIET
% 2=CRRESPRO ACTIVE
% 3=CRRESELE AVERAGE
% 4=CRRESELE WORST CASE
% 5=CRRESELE Ap15
%
%

if isempty(crres_path)
    testfile = 'crrespro_quiet.txt';
    crres_path = which(testfile);
    if isempty(crres_path)
        error('Unable to locate crres files (%s) in "%s"',testfile,mfilename);
    end
    crres_path = crres_path(1:(end-length(testfile)));
end

NE = size(energy,1); % number of energies
NEmax = 25; % maximum number of energies

ntime = length(BBo);
Nmax = onera_desp_lib_ntime_max; % maximum array size in fortran library
Flux = nan(Nmax,1);
if isempty(Ap15)
    Ap15 = nan(1,Nmax);
elseif length(Ap15)==1
    Ap15 = repmat(Ap15,1,Nmax);
else
    Ap15 = [Ap15(:)', nan(1,Nmax-ntime)];
end
if (ntime>Nmax) || (NE>NEmax)
    % break up the calculation into chunks the libarary can handle
    for i = 1:Nmax:ntime
        ii = i:min(i+Nmax-1,ntime);
        for ie = 1:NEmax:NE
            iie = ie:min(ie+NEmax-1,NE);
            Flux(ii,iie) = onera_desp_lib_get_crres_flux(whichm,energy(iie,:),BBo(ii),L(ii),Ap15(ii),crres_path);
        end
    end
else
    %
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
    crres_pathlen = length(crres_path);
    calllib('onera_desp_lib','get_crres_flux_',ntime,iwhichm,whatf,NE,energy',BBo,L,Ap15,FluxPtr,crres_path,crres_pathlen);
    Flux = get(FluxPtr,'value');
end
Flux = Flux(1:ntime,1:NE);
Flux(Flux<-1e30) = nan;
