function out = onera_desp_lib_lstar_phi(which,options,matlabd,in)
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
% function out = onera_desp_lib_lstar_phi(which,options,matlabd,in)
% converts Lstar to Phi or vice versa
% which: 'Lstar2Phi' or 'Phi2Lstar'
% options: defined as usual for onera_desp_lib_options
% matlabd: matlab date number (N x 1)
% in: Lstar or Phi (N x 1), depending on which
% out: Lstar or Phi (N x 1), depending on which

switch(lower(which))
    case {1,'lstar2phi'}
        whichinv = 1;
    case {2,'phi2lstar'}
        whichinv = 2;
    otherwise
        if isnumeric(which)
            error('which = "%g" not supported',which);
        else
            error('which = "%s" not supported',which);
        end
end

matlabd = datenum(matlabd);

onera_desp_lib_load;

if (numel(matlabd)==1) && (numel(in)>1)
    matlabd = repmat(matlabd,size(in));
end

if (numel(in)==1) && (numel(matlabd)>1)
    in = repmat(in,size(matlabd));
end

ntime = numel(in);

options = onera_desp_lib_options(options);

Nmax = onera_desp_lib_ntime_max; % maximum array size in fortran library
if ntime > Nmax % break up into multiple calls
    out = nan(ntime,1);
    for i = 1:Nmax:ntime
        ii = i:min(i+Nmax-1,ntime);
        [out(ii)] = onera_desp_lib_get_field(which,options,matlabd(ii),in(ii));
    end
    return
end

[iyear,idoy] = onera_desp_lib_matlabd2yds(matlabd); % 3rd output, UT, not used
LstarPtr = libpointer('doublePtr',in);
PhiPtr = libpointer('doublePtr',in);

calllib('onera_desp_lib','lstar_phi1_',ntime,whichinv,options,iyear,idoy,LstarPtr,PhiPtr);
% have to do this next bit because Ptr's aren't really pointers

if (whichinv == 1)%lstar2phi
    out = get(PhiPtr,'value');
else % phi2lstar
    out = get(LstarPtr,'value');
end

% the flag value is actually -1d31
out(out<-1e30) = nan;

out = reshape(out,size(in));

