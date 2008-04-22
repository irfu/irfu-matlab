function dobj = dataobj(varargin)
%DATAOBJ  constructor function for DATAOBJ object
%
% DATAOBJ(FILENAME)
%    Construct dataobj form file FILENAME. FILENAME can also contain
%    wildcards ('*').
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
			
			if findstr(varargin{1},'*')
				cdf_files = dir(varargin{1});
					
				switch numel(cdf_files)
					case 0
						error('no cdf files specified')
					case 1
						cdf_file = cdf_files.name;
					otherwise
						% remove '.' and '..' from the list
						j = 1;
						while j<=numel(cdf_files)
							if strcmp(cdf_files(j).name,'.') || ...
									strcmp(cdf_files(j).name,'..')
								cdf_files(j) = [];
							else j = j + 1;
							end
						end
						for j=1:numel(cdf_files),
							disp([num2str(j) '. ' cdf_files(j).name]);
						end
						disp('Choose cdf file');
						j = irf_ask('Cdf_file? [%]>','cdf_file',1);
						cdf_file = cdf_files(j).name;
				end
				clear cdf_files
			else cdf_file = varargin{1};
			end
			
			if ~exist(cdf_file,'file')
				error(['file ' cdf_file ' does not exist'])
			end
			
			[data,info] = cdfread(cdf_file);
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
					% Replace minuses with underscores
					dobj.vars{v,1}(findstr(dobj.vars{v,1},'-')) = '_';
					% Remove training dots
					while (dobj.vars{v,1}(end) == '.')
						dobj.vars{v,1}(end) = [];
					end
					% Take care of '...'
					d3 = findstr(dobj.vars{v,1},'...');
					if d3, dobj.vars{v,1}( d3 + (1:2) ) = []; end
					% Replace dots with underscores
					dobj.vars{v,1}(findstr(dobj.vars{v,1},'.')) = '_';
					% Take care of names longer than 63 symbols (Matlab limit) 
					if length(dobj.vars{v,1})>63
						dobj.vars{v,1} = dobj.vars{v,1}(1:63);
						disp(['orig var : ' dobj.vars{v,2}])
						disp(['new var  : ' dobj.vars{v,1}])
					end
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
