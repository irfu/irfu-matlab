function out = onera_desp_lib_msis(whichm,date,X,sysaxes,F107A,F107,Ap)
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
% out = onera_desp_lib_msis(whichm,date,X,sysaxes,F107A,F107,Ap)
% access to MSIS neutral atmosphere models
% whichm: select which MSIS to use
%   '86' = 'msis86'
%   '90' = 'msise90'
%   '00' = 'nrlmsise00'
% date: Nx1 universal time of interest (matlab date number from datenum)
% X: Nx3 matrix of positions of interes
% sysaxes, coordinate system for X, as defined in onera_desp_lib_sysaxes
% F107A: 3 month average of F10.7 flux (leading up to date)
% F107: daily F10.7 flux for previous day (i.e., before date)
% Ap: Nx1 or Nx7
%  Ap(:,1) = Daily Ap
%  Ap(:,2) = 3 hours Ap index for current time
%  Ap(:,3) = 3 hours Ap index for 3 hours before date
%  Ap(:,4) = 3 hours Ap index for 6 hours before date
%  Ap(:,5) = 3 hours Ap index for 9 hours before date
%  Ap(:,6) = Average of eight 3 hours Ap indices from 12 to 33 hours prior to date
%  Ap(:,7) = Average of eight 3 hours Ap indices from 36 to 59 hours prior to date
%
% out: structure of
%   densities in cm^-3:
%     out.He, out.O, out.N2, out.O2, out.Ar, out.H, out.N, out.AnomalousOxygen
%   out.TotalMass = mass density in g/cm^3
%   out.ExoTemp = Exosphere temperature
%   out.AltTemp = Temperature at altitude
%

out = [];

sysaxes = onera_desp_lib_sysaxes(sysaxes);
X = onera_desp_lib_coord_trans(X,[sysaxes 0],date); % get alt, lat, lon, sysaxes = 0!  GDZ

Nmax = onera_desp_lib_ntime_max;
N = length(date);
if N > Nmax
    for i = 1:Nmax:N
        ii = (i:min(i+Nmax-1,N))';
        % change the call to sysaxes below to 'gdz', since we rotated to
        % GDZ above.  
        tmpout = onera_desp_lib_msis(whichm,date(ii),X(ii,:),'gdz',F107A(ii),F107(ii),Ap(ii,:));
        if isempty(out)
            out = tmpout;
        else
            fldnames = fieldnames(tmpout);
            for ifld = 1:length(fldnames)
                var = fldnames{ifld};
                if size(tmpout.(var),1)==length(ii)
                    out.(var) = [out.(var);tmpout.(var)];
                end
            end
        end
    end
else

    switch(lower(whichm))
        case {'86','msis86'}
            libfunc_str = 'msis86_';
            DensVars = {'He','O','N2','O2','Ar','TotalMass','H','N'};
        case {'90','msise90'}
            libfunc_str = 'msise90_';
            DensVars = {'He','O','N2','O2','Ar','TotalMass','H','N'};
        case {'00','nrlmsise00'}
            libfunc_str = 'nrlmsise00_';
            % NRLMSISE-00 has one additional component, AnomalousOxygen.
            DensVars = {'He','O','N2','O2','Ar','TotalMass','H','N','AnomalousOxygen'};
        otherwise
            error('%s: Unknown "whichm" = "%s"',mfilename,whichm);
    end

    nSpecies = length(DensVars);  % generalize the number of components to accomodate NRLMSISE-00
    
    if size(Ap,2)==7
        WhichAp = 2;
    elseif size(Ap,2) ==1
        WhichAp = 1;
    else
        error('%s: Ap wrong size. Must be Nx1 or Nx7',mfilename);
    end

    nanpad = nan(Nmax-N,1);
    date = [date(:);nanpad];
    [iyear,idoy,UT] = onera_desp_lib_matlabd2yds(date);
    Ap = [Ap;repmat(nanpad,1,size(Ap,2))];
    if WhichAp==1
        Ap = [Ap,nan(Nmax,6)]; % add extra columns
    end
    % Ap in FORTRAN is 7 x Nmax
    F107A = [F107A;nanpad];
    F107 = [F107;nanpad];
    X = [X;repmat(nanpad,1,3)]; % N x 3
    DensPtr = libpointer('doublePtr',nan(Nmax,nSpecies)'); % in FORTRAN 8 x Nmax
    TempPtr = libpointer('doublePtr',nan(Nmax,2)'); % in FORTRAN 2 x Nmax

    calllib('onera_desp_lib',libfunc_str,N,WhichAp,idoy,UT,X(:,1),X(:,2),X(:,3),F107A,F107,Ap',DensPtr,TempPtr);
    Dens = get(DensPtr,'value')'; % Nmax x 8
    Temp = get(TempPtr,'value')'; % Nmax x 2
    Dens = Dens(1:N,:);
    Dens(Dens < 0) = nan; % flag is -1e30
    Temp = Temp(1:N,:);
    Temp(Temp < 0) = nan; % flag is -1e30
    for i = 1:length(DensVars)
        out.(DensVars{i}) = Dens(:,i);
    end
    out.ExoTemp = Temp(:,1);
    out.AltTemp = Temp(:,2);
    out.whichm = whichm;
end
