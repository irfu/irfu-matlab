function cdb = ClusterDB(varargin)
%ClusterDB constructor function for ClusterDB object
% cdb = ClusterDB(isdat_db_s, data_path_s, storage_path_s)
%
% Input:
%	isdat_db_s - ISDAT database, default: 'localhost:10'
%	data_path_s - path to Cluster data, default: '/data/cluster'
%	storage_path_s - where to save the data, default: './'
%
% ClusterDB is the main interface for getting EFW LX, HX and auxiliary
% data from ISDAT, as well as CSDS data.
%
% ClusterDB properties: db, dp, sp
%
% $Revision$  $Date$

% Copyright 2004 Yuri Khotyaintsev

switch nargin
case 0
% if no input arguments, create a default object
	cdb.db = 'localhost:10';
	cdb.dp = '/data/cluster';
	cdb.sp = '.';  
	cdb = class(cdb,'ClusterDB');
case 1
% if single argument of class ClusterDB, return it
	if (isa(varargin{1},'ClusterDB'))
		cdb = varargin{1};
	elseif ischar(varargin{1})
		cdb.db = varargin{1};
		cdb.dp = '/data/cluster';
		cdb.sp = '.';
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
