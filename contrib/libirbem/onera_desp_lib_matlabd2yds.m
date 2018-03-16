function [iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd)
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
% function [iyear,idoy,UT] = onera_desp_lib_matlabd2yds(matlabd);
% helper function to convert between matlab date number and
% date format expected by onera_desp_lib routines
matlabd = datenum(matlabd);
dvec = datevec(matlabd);
iyear = dvec(:,1);
idoy = floor(matlabd)-datenum(iyear,1,1)+1;
UT = floor(dvec(:,4)*60*60+dvec(:,5)*60+dvec(:,6));
