function varargout = onera_desp_lib_landi2lstar(varargin)
%***************************************************************************************************
% Copyright 2006-2009, T.P. O'Brien
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
% function [Lm,Lstar,Blocal,Bmin,J,MLT] = onera_desp_lib_landi2lstar(kext,options,sysaxes,matlabd,x1,x2,x3,maginput)
% faster version of onera_desp_lib_make_lstar
% Works for Olson-Pfitzer Quiet only (ignores kext and field model options)

varargout = cell(1,nargout);
[varargout{:}] = onera_desp_lib_make_lstar_core(mfilename,varargin{:});

