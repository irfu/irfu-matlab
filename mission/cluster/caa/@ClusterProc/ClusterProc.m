function O = ClusterProc(varargin)
%ClusterProc constructor function for ClusterProc object
% O = ClusterProc(storage_path_s,ControlDB)
%
% Input:
%	storage_path_s - where to save the data, default: './'
%	ControlDB is not implemented yet
%
% ClusterProc is the main interface for producing scientific data from
% EFW LX, HX data, as well as CSDS data.
%
% ClusterProc properties: sp
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

switch nargin
case 0
% if no input arguments, create a default object
	O.sp = '.';  
	O = class(O,'ClusterProc');
case 1
% if single argument of class ClusterProc, return it
	if (isa(varargin{1},'ClusterProc'))
		O = varargin{1};
	elseif ischar(varargin{1})
		O.sp = varargin{1};
		O = class(O,'ClusterProc');
	else
		error('Wrong argument type')
	end 
otherwise
	error('Wrong number of input arguments')
end
