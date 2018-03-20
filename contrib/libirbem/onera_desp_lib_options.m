function options  = onera_desp_lib_options(varargin)
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
% function options  = onera_desp_lib_options(options)
% converts cell array of option words into array of numbers
% leaves array of 5 numbers unchaged
% options: array(5) of long integer to set some control options on computed values
% defaults are 0
%
% options(1st element):  0 - don't compute L*;  1 - compute L*
% 'noLstar' - 0, 'doLstar' - 1, 'makePHI' - 2
% options(2nd element): 0 - initialize IGRF field once per year (year.5);
% n - n is the  frequency (in days) starting on January 1st of each year
% (i.e. if options(2nd element)=15 then IGRF will be updated on the
% following days of the year: 1, 15, 30, 45 ...)
% 'IGRFinit',n (two entries in options)
% options(3rd element): resolution to compute L* (0 to 9) where 0 is the
%   recomended value to ensure a good ratio precision/computation time
%   (i.e. an error of ~2% at L=6). The higher the value the better will
%   be the precision, the longer will be the computing time. Generally
%   there is not much improvement for values larger than 4.
% 'L*T',n (two entries in options)
% 'L*3',n (two entries in options)
% options(4th element): resolution to compute L* (0 to 9). The higher the
%   value the better will be the precision, the longer will be the computing time.
%   It is recommended to use 0 (usually sufficient) unless L* is not computed
%   on a LEO orbit. For LEO orbit higher values are recommended.
% 'L*R',n (two entries in options)
% 'L*4',n (two entries in options)
% options(5th element): allows to select an internal magnetic field model (default is set to IGRF)
% 'IGRF' 0 = IGRF
% 'ECD'  1 = Eccentric tilted dipole
% 'JC'   2 = Jensen&Cain 1960
% 'GSFC' 3 = GSFC 12/66 updated to 1970
% 'USER' 4 = user's own, compiled into library as myOwnMagField
% 'DIPOLE' 5 = centered dipole

options = zeros(5,1);
if nargin==0
    return;
end
inoptions = varargin{1};
if isnumeric(inoptions) && (length(inoptions)==5)
    options = inoptions;
    return;
end

if isnumeric(inoptions) && (nargin==5)
    for i = 1:5
        if isnumeric(varargin{i})
            options(i) = varargin{i};
        end
    end
end

if ~iscell(inoptions)
    inoptions = varargin;
end

i = 1;
while i <= length(inoptions)
    opt = inoptions{i};
    if ~isempty(opt) && ischar(opt)
        switch(upper(opt))
            case {'NOLSTAR'}, options(1) = 0;
            case {'DOLSTAR'}, options(1) = 1;
            case {'MAKEPHI'}, options(1) = 2;
            case {'IGRFINIT'}, options(2) = eval_pair(inoptions,i); i = i+1;
            case {'L*T','L*3'}, options(3) = eval_pair(inoptions,i); i = i+1;
            case {'L*R','L*4'}, options(4) = eval_pair(inoptions,i); i = i+1;
            case {'IGRF'}, options(5) = 0;
            case {'ECD'}, options(5) = 1;
            case {'JC'}, options(5) = 2;
            case {'GSFC'}, options(5) = 3;
            case {'USER'}, options(5) = 4;
            case {'DIPOLE'}, options(5) = 5;
            otherwise
                error('Unknown option word "%s" in %s',opt,mfilename);
        end
    end
    i = i+1;
end

function n = eval_pair(inoptions,i)
if i>=length(inoptions)
    error('Expected another option after "%s" in %s',inoptions{i},mfilename);
end
n = inoptions{i+1};
if ~isnumeric(n)
    n = str2num(n);
end
