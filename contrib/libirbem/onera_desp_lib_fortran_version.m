function V = onera_desp_lib_fortran_version()
%***************************************************************************************************
% Copyright 2009, T.P. O'Brien
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
% V = onera_desp_lib_fortran_version()
% return version number of fortran code (repository revision number, not
% release number)

onera_desp_lib_load;

vPtr = libpointer('int32Ptr',-1);
calllib('onera_desp_lib','irbem_fortran_version1_',vPtr);
V = double(vPtr.value);

