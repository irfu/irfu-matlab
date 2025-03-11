function igrf_version = onera_desp_lib_igrf_version()
%***************************************************************************************************
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
% igrf_version = onera_desp_lib_igrf_version()
% size of ntime dimension in fortran arrays

onera_desp_lib_load;

nPtr = libpointer('int32Ptr',-1);
calllib('onera_desp_lib','get_igrf_version_',nPtr);
igrf_version = double(nPtr.value);

