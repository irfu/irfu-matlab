function cdb = ClusterDB(varargin)
%ClusterDB  constructor function for ClusterDB object
%
% cdb = ClusterDB(isdat_db_s, data_path_s, storage_path_s)
% cdb = ClusterDB(storage_path_s)
% cdb = ClusterDB(isdat_db_s, data_path_s)
%
% Input:
%	isdat_db_s - ISDAT database strings separated by '|', 
%		default: from c_ctl 
%	data_path_s - path to Cluster data, default: from c_ctl 
%	storage_path_s - where to save the data, default: './'
%
% ClusterDB is the main interface for getting EFW LX, HX and auxiliary
% data from ISDAT, as well as CSDS data.
%
% ClusterDB properties: db, dp, sp
%
% $Revision$  $Date$

% Copyright 2004,2005 Yuri Khotyaintsev

switch nargin
case 0
% if no input arguments, create a default object
	cdb.db = c_ctl(0,'isdat_db');
	cdb.dp = c_ctl(0,'data_path');
	cdb.sp = '.';  
	cdb = class(cdb,'ClusterDB');
case 1
% if single argument of class ClusterDB, return it
	if (isa(varargin{1},'ClusterDB'))
		cdb = varargin{1};
	elseif ischar(varargin{1})
	  cdb.db = c_ctl(0,'isdat_db');
	  cdb.dp = c_ctl(0,'data_path');
		cdb.sp = varargin{1};
		cdb = class(cdb,'ClusterDB');
	else
		error('Wrong argument type')
	end 
case 2
	if ischar(varargin{1}) && ischar(varargin{2})
		cdb.db = varargin{1};
		cdb.dp = varargin{2};
		cdb.sp = '.';
		cdb = class(cdb,'ClusterDB');
	else
		error('Wrong argument type')
	end 
case 3
% create object using specified values
	cdb.db = varargin{1};
	cdb.dp = varargin{2};
	cdb.sp = varargin{3};
	cdb = class(cdb,'ClusterDB');
otherwise
	error('Wrong number of input arguments')
end
