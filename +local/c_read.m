function [out,dataobject]=c_read(varargin)
% LOCAL.C_READ read local cluster aux information
%	[out]=LOCAL.C_READ(variable,tint)
%		read variable for given time interval tint in matlab format (matrix)
%		tint - [tstart tend] or ISO format, see example
%	[var,dataobj]=LOCAL.C_READ(variable,tint,'caa')
%		read variable in CAA format.
%		var - variable, dataobj - dataobject
%
%	LOCAL.C_READ('list') list all datasets that are locally available and indexed
%		To see variables for any data set see LOCAL.C_CAA_META
%
% Variable can be CAA variable or shortcuts
%  'R1'  - Cluster 1 position
%  'dR1' - Cluster 1 relative position wrt center
%   'R'   - Cluster center position and position of all s/c into structure out
%
% Example variables: (question mark needs to be sustituted to Cluster number)
% sc_orbit_num__CL_SP_AUX
% sc_r_xyz_gse__CL_SP_AUX
% sc_v_xyz_gse__CL_SP_AUX
% sc_dr?_xyz_gse__CL_SP_AUX
% gse_gsm__CL_SP_AUX
% dipole_tilt__CL_SP_AUX
% sc_geom_size__CL_SP_AUX
% sc_geom_elong__CL_SP_AUX
% Invar_Lat__C1_JP_PMP
% Mag_Local_time__C1_JP_PMP
% L_value__C1_JP_PMP
%
%	Examples:
%		tint = '2005-01-01T05:00:00.000Z/2005-01-05T05:10:00.000Z';
%		   R = local.c_read('r',tint);
%		  R1 = local.c_read('R1',tint);
% DipoleTilt = local.c_read('dipole_tilt__CL_SP_AUX',tint);
%
% to update file index run LOCAL.C_UPDATE

persistent index % to make fast access read only once
persistent usingNasaPatchCdf

if isempty(usingNasaPatchCdf), % check only once if using NASA cdf
	usingNasaPatchCdf=irf.check_if_using_nasa_cdf;
end

%% Defaults
returnDataFormat = 'mat'; % default matlab format
if ispc
    caaDir='Z:\';
else
    caaDir='/data/caalocal/';
end
%% Default index is empty, read in only those indees that are used

if isempty(index), 
	index=struct('dummy',[]);
end
%% Check inputs
if nargin == 0,
	help local.c_read;
	return;
end
if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'list')
	list_indexed_datasets;
	return;
elseif nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'test')
	out = false;
	if exist(caaDir,'dir'), out = true; end
	return
elseif nargin>=2,
	varName=varargin{1};
	tint=varargin{2};
	if ischar(tint),
		tint=irf_time(tint,'iso2tint');
	end
end
if nargin ==3,
	if ischar(varargin{3}) && any(strcmpi(varargin{3},'dobj'))
		returnDataFormat = 'dobj';
	elseif ischar(varargin{3}) && any(strcmpi(varargin{3},'caa'))
		returnDataFormat = 'caa';
	elseif ischar(varargin{3}) && any(strcmpi(varargin{3},'mat'))
		returnDataFormat = 'mat';
	else
		irf_log('fcal','output data format unknown');
		out=[];
		return;
	end
end
if nargin > 3
	irf_log('fcal','max 3 arguments supported');
	return
end
%% Check if repository is there
if ~exist(caaDir,'dir')
	disp(['Local CAA data repository ' caaDir ' not available!']);
	return;
end
%% Read in data
out=[];dataobject=[]; % default return empty output
specialCaseCis=0;
switch lower(varName)
	case {'r'}
		varToRead={'sc_r_xyz_gse__CL_SP_AUX','sc_dr1_xyz_gse__CL_SP_AUX',...
			'sc_dr2_xyz_gse__CL_SP_AUX','sc_dr3_xyz_gse__CL_SP_AUX','sc_dr4_xyz_gse__CL_SP_AUX'};
		ok=readdata;
		if ok && strcmpi(returnDataFormat,'mat')
			out.R=[data{1} double(data{2})];
			c_eval('out.R?=[data{1} double(data{2}+data{2+?})];')
		end
	case {'r1','r2','r3','r4'}
		varToRead={'sc_r_xyz_gse__CL_SP_AUX',['sc_dr' varName(2) '_xyz_gse__CL_SP_AUX']};
		ok=readdata;
		if ok && strcmpi(returnDataFormat,'mat'),
			out=[data{1} double(data{2}+data{3})];
		end
	case {'dr1','dr2','dr3','dr4'}
		varToRead={['sc_dr' varName(3) '_xyz_gse__CL_SP_AUX']};
		ok=readdata;
		if ok && strcmpi(returnDataFormat,'mat'),
			out=[data{1} double(data{2})];
		end
	otherwise
		irf_log('fcal',['Reading variable (assume to exist): ' varName]);
		if strfind(varName,'CIS'),specialCaseCis=1;end
		varToRead={varName};
		ok=readdata;
		if ok && strcmpi(returnDataFormat,'mat'),
            if numel(data)==2 && numel(size(data{2}))==2,
                out=[data{1} double(data{2})];
            elseif numel(data)==1,
                out = data{1};
			else
				out=data;
            end
		  elseif ok && (strcmpi(returnDataFormat,'dobj') || strcmpi(returnDataFormat,'caa')),
			out = data;
		end
end

	function status=readdata
		status = false; % default 
		%% find index
		ii=strfind(varToRead{1},'__');
		if ii,
			dataset=varToRead{1}(ii+2:end);
			datasetIndex = strrep(dataset,'CIS-','CIS_');
			if ~isfield(index,datasetIndex) % index not yet read
				indexVarName = ['index_' datasetIndex];
				indexFile = [caaDir 'caa'];
				indexFileInfo=whos('-file',indexFile,indexVarName);
				if numel(indexFileInfo)==0, % there is no index
					irf_log('dsrc',['There is no index file:' indexVarName]);
					return;
				end
				s=load(indexFile,indexVarName);
				index.(datasetIndex)=s.(indexVarName);
			end
			index=index.(datasetIndex);
		else
			irf_log('dsrc',['Do not know how to read variable: ' varToRead{1}]);
			return
		end
		%% find files within time interval
		istart=find(index.tend>tint(1),1);
		iend=find(index.tstart<tint(2),1,'last');
		irf_log('dsrc',['Dataset: ' dataset '. Index files: ' num2str(istart) '-' num2str(iend)]);

		if isempty(istart) || isempty(iend) || istart > iend,
			return
		end
		%% read in records
		for iFile=istart:iend
			cdf_file=[caaDir index.filename(iFile,:)];
			if specialCaseCis,
				dataset=strrep(dataset,'CIS_','CIS-');
				varToRead=strrep(varToRead,'CIS_','CIS-');
			end
			switch lower(returnDataFormat)
				case 'mat'
					irf_log('dsrc',['Reading: ' cdf_file]);
					%% check if epoch16
					cdfid=cdflib.open(cdf_file);
					useCdfepoch16=strcmpi(cdflib.inquireVar(cdfid,0).datatype,'cdf_epoch16');
					if useCdfepoch16,
						irf_log('dsrc',['EPOCH16 time in cdf file:' cdf_file]);
						tmptime=readCdfepoch16(cdfid,0); % read time which has variable number 0
						timeVector=irf_time(tmptime,'cdfepoch162epoch');
						tmpdata=cell(1,numel(varToRead));
						for iVar=1:numel(varToRead),
                            tmpdata{iVar}=cdfread(cdf_file,'CombineRecords',true,...
                                'Variables',varToRead{iVar});
						end
						tmpdata = [{timeVector} tmpdata]; %#ok<AGROW>
                    else
                        % remove time variable as it is already read in
                        ii=numel(varToRead);
                        while ii,
                            if strcmp(varToRead{ii},cdflib.getVarName(cdfid,0)),
                                varToRead(ii)=[];
                            end
                            ii=ii-1;
                        end
                        % read data
                        [tmpdata,~] = cdfread(cdf_file,'ConvertEpochToDatenum',true,'CombineRecords',true,...
							'Variables', [{cdflib.getVarName(cdfid,0)},varToRead{:}]); % time and variable name
						tmpdata=fix_order_of_array_dimensions(tmpdata);
						if isnumeric(tmpdata), tmpdata={tmpdata}; end % make cell in case matrix returned
                        timeVector = irf_time(tmpdata{1},'date2epoch');
						tmpdata{1} = timeVector;
					end
					if iFile==istart, data=cell(size(tmpdata));end
					iist=1;iien=numel(timeVector);
					if iFile==istart
						iist=find(timeVector>tint(1),1);
					end
					if iFile==iend
						iien=find(timeVector<tint(2),1,'last');
					end
					%% check for NaNs
					for iVar=1:numel(varToRead),
						fillVal=value_of_variable_attribute(cdfid,varToRead{iVar},'FILLVAL');
						tmpdata{iVar+1}(tmpdata{iVar+1}==fillVal)=NaN; % +1 because first cell is time and then comes variables
					end
					%% attach to result
					for j=1:numel(data),
						nDim=numel(size(tmpdata{j}));
						if nDim==2,
							data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:));
						elseif nDim==3,
							data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:,:));
						elseif nDim==4,
							data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:,:,:));
						elseif nDim==5,
							data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:,:,:,:));
						end
					end
					cdflib.close(cdfid);
				case {'caa','dobj'}
					if iFile==istart % start of interval, initiate dataobject
                        dataobject=dataobj(cdf_file,'tint',tint);
                    elseif iFile==iend % ends of interval
						data_temp=dataobj(cdf_file,'tint',tint);
					else
						data_temp=dataobj(cdf_file);
					end
					if exist('data_temp','var') && ~isempty(data_temp)
						if isempty(dataobject)
							dataobject=data_temp;
						else
							dataobject=append(dataobject,data_temp);
						end
						clear data_temp;
					end
				otherwise
					error('unknown format');
			end
		end
		if strcmp(returnDataFormat,'caa')
			data=get(dataobject,varToRead{1}); % currently only 1 variable request implemented
		  elseif strcmp(returnDataFormat,'dobj') 
			data=dataobject; % currently only 1 variable request implemented
		end
		status = true;
	end
	function data = readCdfepoch16(cdfid,varName)
		if isnumeric(varName),
			varnum=varName;
		elseif ischar(varName)
			varnum  = cdflib.getVarNum(cdfid,varName);
		else
			error('varName should be variable number or name');
		end
		numrecs = cdflib.getVarNumRecsWritten(cdfid,varnum);
		numElements = cdflib.inquireVar(cdfid,varnum).numElements;
		data=zeros(1+numElements,numrecs);
		
		for j = 0:numrecs-1
			% This reads in the data in raw epoch 16 format
			% That implies each Epoch16 value is a 2-element double precision
			% value in MATLAB
			data(:,1+j) = cdflib.getVarRecordData(cdfid,varnum,j);
		end
		data=data';
	end
	function value=value_of_variable_attribute(cdfid,varName,attrName)
		attrnum = cdflib.getAttrNum(cdfid,attrName);
		varnum = cdflib.getVarNum(cdfid,varName);
		value = cdflib.getAttrEntry(cdfid,attrnum,varnum);
	end
	function ok=list_indexed_datasets
		ok=false;
		s=whos('-file',[caaDir 'caa']);
		sind=arrayfun(@(x) any(strfind(x.name,'index_')),s);
		for jj=1:numel(sind)
			if sind(jj)
				disp(s(jj).name(7:end));
			end
		end
		if any(sind), ok=true; end
	end
	function data=fix_order_of_array_dimensions(data)
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

end
