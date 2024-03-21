function info = onera_desp_lib_load(libfile,headerfile)
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
%
%***************************************************************************************************
%
%function onera_desp_lib_load(libfile,@prototypefunction);
%function onera_desp_lib_load(libfile,headerfile);
% checks for the presence of the onera_desp_lib dynamic library in memory
% if not present, attempts to load it using a headerfile (libirbem.h)
% checks environment variable IRBEM_LIB_DLL to specify the file to load
%  (full file path, or a file name in the matlab search path)
% otherwise, it guesses the file name
% checks the environment variable IRBEM_THUNK_TMP_PATH to specify where to 
% create the thunkfile (64-bit systems only), otherwise let's matlab decide

if ~libisloaded('onera_desp_lib')
    if nargin < 2
        headerfile = 'libirbem.h';
    end
    
    % determine DLL extension
    if ispc
        libext = 'dll';
    elseif ismac
        libext = 'dylib';
    else
        libext = 'so';
    end

    if nargin < 1
        if ~isempty(getenv('IRBEM_LIB_DLL'))
            libfile = getenv('IRBEM_LIB_DLL');
        else
            libfile = ['libirbem_', lower(computer),'.',libext];
        end
    end
    if ~exist(libfile,'file')
        error('libfile %s not found',libfile);
    end
    %fprintf('Loading %s\n',libfile);
    if isempty(getenv('IRBEM_THUNK_TMP_PATH'))
        loadlibrary(libfile,headerfile,'alias','onera_desp_lib'); % let matlab choose the thunkfile name and location
    else % use user specifiedlocation
        old_pwd = pwd; % current wd
        cd(getenv('IRBEM_THUNK_TMP_PATH')); % cd into temp folder
        tfile = ['onera_desp_lib_thunk_',lower(computer),'.',libext];
        loadlibrary(libfile,headerfile,'alias','onera_desp_lib','thunkfilename',tfile); % load with named thunkfile
        cd(old_pwd); % cd back into original wd
    end
    
end


% command to generate proto file -- don't do this anymore as it creates problems for 64-bit machines
% loadlibrary('onera_desp_lib.dll','onera_desp_lib.h','alias','onera_desp_lib','mfilename','onera_desp_lib_proto.m'); disp('remember to move the proto file');
