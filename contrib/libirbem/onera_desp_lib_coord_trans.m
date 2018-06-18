function Y = onera_desp_lib_coord_trans(X,rotation,matlabd)
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
%***************************************************************************************************
%
% Y = onera_desp_lib_coord_trans(X,rotation,matlabd)
% rotates X ([nx3]) into new coordinate system
% Y is nx3 and psi (for geo2gsm/gsm2geo) is nx1
% rotation x2y
% sph, car, geo, gsm, gse, gdz, gei (eci), sm, mag, rll, gse, hee, hae, heeq
% (car=geo, but used for generic sph2car and car2sph)
% sph/rll coordinates are ordered: r, lat, lon
% gdz coordinates are ordered: alt, lat, lon (as in sysaxes=0)
% also can specify numerical coordinate IDs like in sysaxes rotation = [num1 num2]
% matlabd (date) is optional, since many coordinate systems are time invariant

if nargin < 3
    matlabd = now;
end

matlabd = datenum(matlabd);

onera_desp_lib_load;

[r,c] = size(X);
if r==3 && c==1
    X = X';
    c = 3;
    r = 1;
end
if c ~=3
    error('Argument X must be size n x 3 in %s',mfilename);
end
ntime = r;
if length(matlabd)==1
    matlabd = repmat(matlabd,ntime,1);
end

if length(matlabd)~=r
    if r==1
        X = repmat(X,length(matlabd),1);
    else
        error('Size of X argument and size of matlabd argument are not commensurate in %s',mfilename);
    end
end

if isnumeric(rotation(1))
    istart = rotation(1);
    iend = rotation(2);
else
    rotation = lower(rotation);
    rotation = strrep(rotation,'car','geo'); % simplify
    rotation = strrep(rotation,'eci','gei'); % simplify
    rotation(rotation=='_') = '2';
    f = find(rotation=='2');
    rstart = lower(rotation(1:(f-1)));
    rend = lower(rotation((f+1):end));
    istart = onera_desp_lib_sysaxes(rstart);
    iend = onera_desp_lib_sysaxes(rend);
end

Y = coord_trans_type(istart,iend,matlabd,X);

% handle flags
Y(Y<-1e30) = nan;

function Y = coord_trans_type(istart,iend,matlabd,X)
[iyear,idoy,secs] = onera_desp_lib_matlabd2yds(matlabd);
ntime = length(iyear);
nmax = onera_desp_lib_ntime_max;
Y = nan(ntime,3);
nstart = 1;
while nstart <= ntime
    % void coord_trans_vec1_(long int *ntime, long int *sysaxesIN,long int *sysaxesOUT,
    % 		   long int *iyr,long int *idoy,double *secs,
    % 		   double *xINV,double *xOUTV);
    %
    ii = nstart:min(nstart+nmax-1,ntime);
    Ynmax = nan(3,nmax); % place holder
    YPtr = libpointer('doublePtr',Ynmax);
    calllib('onera_desp_lib','coord_trans_vec1_',length(ii),istart,iend,iyear(ii),idoy(ii),secs(ii),X(ii,:)',YPtr);
    Ynmax = get(YPtr,'value');
    Y(ii,:) = Ynmax(:,1:length(ii))';
    nstart = nstart+nmax;
end
