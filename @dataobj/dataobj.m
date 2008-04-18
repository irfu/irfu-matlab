function dobj = dataobj(varargin)
%DATAOBJ  constructor function for dataobj object
%
% DATAOBJ properties: 
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

switch nargin
	case 0
		% if no input arguments, create a default object
		dobj.FileModDate = datestr(now);
		dobj.VariableAttributes = {};
		dobj.GlobalAttributes = {};
		dobj.Variables = {};
		dobj.vars = {};
		dobj = class(dobj,'dataobj');
	case 1
		% if single argument of class ClusterDB, return it
		if (isa(varargin{1},'dataobj'))
			dobj = varargin{1};
		elseif ischar(varargin{1})
			if ~exist(varargin{1},'file')
				error(['file ' varargin{1} ' does not exist'])
			end
			[data,info] = cdfread(varargin{1});
			dobj.FileModDate = info.FileModDate;
			dobj.VariableAttributes = info.VariableAttributes;
			dobj.GlobalAttributes = info.GlobalAttributes;
			dobj.Variables = info.Variables;
			nvars = size(info.Variables,1);
			dobj.vars = cell(nvars,2);
			if nvars>0
				dobj.vars(:,1) = info.Variables(:,1);
				dobj.vars(:,2) = info.Variables(:,1);
				for v=1:nvars
					dobj.vars{v,1}(findstr(dobj.vars{v,1},'-')) = '_'; % replace minuses with underscores
					dobj.data.(dobj.vars{v,1}).data = [data{:,v}];
					dobj.data.(dobj.vars{v,1}).dim = info.Variables{v,2};
					dobj.data.(dobj.vars{v,1}).nrec = info.Variables{v,3};
					dobj.data.(dobj.vars{v,1}).type = info.Variables{v,4};
					dobj.data.(dobj.vars{v,1}).variance = info.Variables{v,5};
					dobj.data.(dobj.vars{v,1}).sparsity = info.Variables{v,6};
					if strcmp(dobj.data.(dobj.vars{v,1}).type,'epoch')
						temp = struct(dobj.data.(dobj.vars{v,1}).data);
						temp = [temp.date]; 
						temp = temp(:);
						temp = (temp-62167219200000)/1000;
						dobj.data.(dobj.vars{v,1}).data = temp;
						clear temp
					end
				end
			else
				dobj.data = [];
			end
			dobj = class(dobj,'dataobj');
		else
			error('Wrong argument type')
		end
	otherwise
		error('Wrong number of input arguments')
end
