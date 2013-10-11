function dobj = dataobj(varargin)
%DATAOBJ  constructor function for DATAOBJ object
%
% DATAOBJ(FILENAME)
%    Construct dataobj form file FILENAME. FILENAME can also contain
%    wildcards ('*').
%
% DATAOBJ(FILENAME,'tint',tint)
%       tint - limit dataobject to time interval (good for large files)
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
persistent usingNasaPatchCdf

if isempty(usingNasaPatchCdf), % check only once if using NASA cdf
	usingNasaPatchCdf=irf.check_if_using_nasa_cdf;
end
shouldReadAllData = true; % default read all data
noDataReturned    = 0;    % default expects data to be returned
if nargin==0, action='create_default_object'; end
if nargin==1, action='read_data_from_file'; end
if nargin==3 && ...
		ischar(varargin{2}) && strcmp(varargin{2},'tint') && ...
		isnumeric(varargin{3}) && (length(varargin{3})==2),
	tint=varargin{3};
	action='read_data_from_file';
	shouldReadAllData=false;
end
switch action
	case 'create_default_object'
		% if no input arguments, create a default object
		dobj.FileModDate = datestr(now);
		dobj.VariableAttributes = {};
		dobj.GlobalAttributes = {};
		dobj.Variables = {};
		dobj.vars = {};
		dobj = class(dobj,'dataobj');
	case 'read_data_from_file'
		% if single argument of class ClusterDB, return it
		if (isa(varargin{1},'dataobj'))
			dobj = varargin{1};			
		elseif ischar(varargin{1})
			if strfind(varargin{1},'*')
				cdf_files = dir(varargin{1});
				cdf_files([cdf_files(:).isdir])=[]; % remove directories from the list
				if strfind(varargin{1},filesep) % if there is directory in file name
					filesep_indexes=strfind(varargin{1},filesep);
					directory_name=varargin{1}(1:filesep_indexes(end));
				else
					directory_name='';
				end
				switch numel(cdf_files)
					case 0
						error('no cdf files specified')
					case 1
						cdf_file = [directory_name cdf_files.name];
					otherwise
						for j=1:numel(cdf_files),
							disp([num2str(j) '. ' cdf_files(j).name]);
						end
						j = irf_ask('Choose cdf-file? [%]>','cdf_file',1);
						cdf_file = [directory_name cdf_files(j).name];
				end
				clear cdf_files
			else cdf_file = varargin{1};
			end
			if ~exist(cdf_file,'file')
				error(['file ' cdf_file ' does not exist'])
			end
			%% read in file
			irf.log(2,['Reading: ' cdf_file]);
			% get basic info
			info   = cdfinfo(cdf_file);
			% check if cdfepoch16
			usingCdfepoch16=strcmpi('epoch16',info.Variables{1,4});
			% initialize data object
			if usingNasaPatchCdf
				[data,info] = cdfread(cdf_file,'CombineRecords',true,'KeepEpochAsIs',true);
				if usingCdfepoch16
					timeline = convert_cdfepoch16_string_to_isdat_epoch(data{1});
				else
					data{1}(data{1}==0)=NaN; % fillvalue timeline
					timeline = irf_time(data{1},'cdfepoch2epoch');
                    timeline = timeline(:); % bug fix for cdfread (time comes out as row vector)
				end
				data{1}=timeline;
				fix_order_of_array_dimensions;
				if ~shouldReadAllData
					records=(timeline > tint(1)) & (timeline < tint(2));
				end
			else
				% read in file
				if usingCdfepoch16,
					irf.log(3,['EPOCH16 time in cdf file:' cdf_file]);
					shouldReadAllData=1; % read all data
					variableNames=info.Variables(:,1);
					isCdfepoch16VariableArray=cellfun(@(x) strcmpi(x,'epoch16'), info.Variables(:,4));
					data=cell(1,size(variableNames,1));
					data(~isCdfepoch16VariableArray) = cdfread(cdf_file,'variables',variableNames(~isCdfepoch16VariableArray),'CombineRecords',true);
					iCdf16Variable = find(isCdfepoch16VariableArray);
					for i=1:length(iCdf16Variable)
						numrecs = info.Variables{iCdf16Variable(i),3};
						% get time axis
						tc=zeros(2,numrecs);
						cdfid   = cdflib.open(cdf_file);
						for jj=1:numrecs,
							tc(:,jj) = cdflib.getVarRecordData(cdfid,iCdf16Variable(i)-1,jj-1);
						end
						cdflib.close(cdfid);
						data(iCdf16Variable(i))={irf_time(tc','cdfepoch162epoch')};
					end
				else
					[data,info] = cdfread(cdf_file,'ConvertEpochToDatenum',true,'CombineRecords',true);
					data{1} = irf_time([data{:,1}],'date2epoch');
				end
				if ~shouldReadAllData,
					records=find((data{1} > tint(1)) & (data{1} < tint(2)));
				end
			end
			%% check if number of records to read is zero
			if ~shouldReadAllData && sum(records)==0,
				irf.log(2,'No data within specified time interval');
				noDataReturned=1;
			end

			%% construct data object
			dobj.FileModDate		= info.FileModDate;
			dobj.VariableAttributes = info.VariableAttributes;
			dobj.GlobalAttributes	= info.GlobalAttributes;
			dobj.Variables			= info.Variables;
            if usingCdfepoch16
                update_variable_attributes_cdfepoch16;
            else
                update_variable_attributes_cdfepoch;
            end
			% test if there are some data
			if ~(any(strcmpi(info.Variables(:,4),'epoch')==1) || ...
					any(strcmpi(info.Variables(:,4),'epoch16')==1)),
				nVariables=0; % no time variable, return nothing
				irf.log(2,'CDF FILE IS EMPTY!')
			else
				nVariables = size(info.Variables,1);
			end
			dobj.vars = cell(nVariables,2);
			if nVariables>0
				dobj.vars(:,1) = info.Variables(:,1); % new variables
				dobj.vars(:,2) = info.Variables(:,1); % original variables
				for v=1:nVariables
					dobj.vars{v,1}=variable_mat_name(dobj.vars{v,2});
					varName = dobj.vars{v,1};
					data_all_records = data{v};
					if shouldReadAllData || ...% return all data
							(usingNasaPatchCdf && strcmpi(info.Variables{v,5}(1),'F')),% fixed variable with NASA cdf patch (matlab cdfread return fixed variable as time series) 
						dobj.data.(varName).data = data_all_records;
						dobj.data.(varName).nrec = info.Variables{v,3};
					else
						nDim=numel(size(data_all_records));
						if nDim==2,
							data_records_within_interval=data_all_records(records,:);
						elseif nDim==3,
							data_records_within_interval=data_all_records(records,:,:);
						elseif nDim==4,
							data_records_within_interval=data_all_records(records,:,:,:);
						elseif nDim==5,
							data_records_within_interval=data_all_records(records,:,:,:,:);
						end
						dobj.data.(varName).data = data_records_within_interval;
						dobj.data.(varName).nrec = sum(records);
						dobj.Variables{v,3}      = sum(records);
					end
					dobj.data.(varName).dim      = info.Variables{v,2};
					dobj.data.(varName).type     = info.Variables{v,4};
					dobj.data.(varName).variance = info.Variables{v,5};
					dobj.data.(varName).sparsity = info.Variables{v,6};
				end
			end
			dobj = class(dobj,'dataobj');
		else
			error('Wrong argument type')
		end
	otherwise
		error('Wrong number of input arguments')
end
	function fix_order_of_array_dimensions
		for iDimension=3:4,
			indDatasets=find(cellfun(@(x) numel(size(x)),data(:))==iDimension); % find iDimension datasets
			for iDataset=1:numel(indDatasets)
				if iDimension==3,
					data{indDatasets(iDataset)}=permute(data{indDatasets(iDataset)},[3 1 2]);
				elseif iDimension==4,
					data{indDatasets(iDataset)}=permute(data{indDatasets(iDataset)},[4 3 1 2]);
				end
			end
		end
	end
	function update_variable_attributes_cdfepoch % nested function
		isFieldUnits        = isfield(dobj.VariableAttributes,'UNITS');
		isFieldSIConversion = isfield(dobj.VariableAttributes,'SI_CONVERSION');
		isFieldDeltaPlus    = isfield(dobj.VariableAttributes,'DELTA_PLUS');
		isFieldDeltaMinus   = isfield(dobj.VariableAttributes,'DELTA_MINUS');
		timeVariable=info.Variables{1,1};
		if isFieldUnits,
			iattr=find(strcmpi(dobj.VariableAttributes.UNITS(:,1),timeVariable));
			if iattr, dobj.VariableAttributes.UNITS(iattr,2)={'s'};end % change from ms to s UNITS of epoch if present
		end
		if isFieldSIConversion,
			iattr=find(strcmpi(dobj.VariableAttributes.SI_CONVERSION(:,1),timeVariable));
			if iattr, dobj.VariableAttributes.SI_CONVERSION(iattr,2)={'1.0>s'};end % change from ms to s SI_CONVERSION of epoch if present
		end
		if isFieldDeltaPlus,
			iattr=find(strcmpi(dobj.VariableAttributes.DELTA_PLUS(:,1),timeVariable)); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_PLUS{iattr,2}),
				dobj.VariableAttributes.DELTA_PLUS{iattr,2}=dobj.VariableAttributes.DELTA_PLUS{iattr,2}/1000;
			end
		end
		if isFieldDeltaMinus,
			iattr=find(strcmpi(dobj.VariableAttributes.DELTA_MINUS(:,1),timeVariable)); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_MINUS{iattr,2}),
				dobj.VariableAttributes.DELTA_MINUS{iattr,2}=dobj.VariableAttributes.DELTA_MINUS{iattr,2}/1000;
			end
		end
	end
	function t=convert_cdfepoch16_string_to_isdat_epoch(in) % nested function
		% the following if is because of the bug in CAA cdf files having EPOCH16
		% sometimes time variable has dimension zero and sometimes one
		% TODO: report bug to CAA team and if needed cdf team
		if numel(size(in)) == 3,
			tcdfepoch=reshape(in,size(in,1),size(in,3)); % cdfread returns (Nsamples X 1 X 2) matrix
		else 
			tcdfepoch = in; % cdfread returns (Nsamples X 2) matrix
		end
		t=irf_time(tcdfepoch,'cdfepoch162epoch');
	end
	function update_variable_attributes_cdfepoch16 % nested function
		isFieldUnits        = isfield(dobj.VariableAttributes,'UNITS');
		isFieldSIConversion = isfield(dobj.VariableAttributes,'SI_CONVERSION');
		isFieldDeltaPlus    = isfield(dobj.VariableAttributes,'DELTA_PLUS');
		isFieldDeltaMinus   = isfield(dobj.VariableAttributes,'DELTA_MINUS');
		timeVariable=info.Variables{1,1};
		if isFieldUnits,
			iattr=find(strcmpi(dobj.VariableAttributes.UNITS(:,1),timeVariable));
			if iattr, dobj.VariableAttributes.UNITS(iattr,2)={'s'};end % change from ms to s UNITS of epoch if present
		end
		if isFieldSIConversion,
			iattr=find(strcmpi(dobj.VariableAttributes.SI_CONVERSION(:,1),timeVariable));
			if iattr, dobj.VariableAttributes.SI_CONVERSION(iattr,2)={'1.0>s'};end % change from ms to s SI_CONVERSION of epoch if present
		end
		if isFieldDeltaPlus,
			iattr=find(strcmpi(dobj.VariableAttributes.DELTA_PLUS(:,1),timeVariable)); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_PLUS{iattr,2}),
				dobj.VariableAttributes.DELTA_PLUS{iattr,2}=dobj.VariableAttributes.DELTA_PLUS{iattr,2}/1e12;
			end
		end
		if isFieldDeltaMinus,
			iattr=find(strcmpi(dobj.VariableAttributes.DELTA_MINUS(:,1),timeVariable)); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_MINUS{iattr,2}),
				dobj.VariableAttributes.DELTA_MINUS{iattr,2}=dobj.VariableAttributes.DELTA_MINUS{iattr,2}/1e12;
			end
		end
	end

end % end of main functions

