function [out,dataobject]=c_read(varargin)
% LOCAL.C_READ read local cluster aux information
%	[out]=LOCAL.C_READ(variableName,tint,[format])
%		read variable variableName for a given time interval tint in a specified format.
%   Default format is 'mat' but in future will be 'ts'.
%   variableName can be a cell array of names in which case cell array of
%   data is returned.
%		tint - [tstart tend] or UTC format or tint in GenericTimeArray
%
%   Possible formats:
%     'mat' - matlab array with the first column being (current default)
%     'ts'  - TSeries object (default in future)
%     'caa' - caa variable format
%    'dobj' - return dataobject instead of variable
%
%	[var,dataobj]=LOCAL.C_READ(variable,tint,'caa')
%		read variable in CAA format.
%		var - variable, dataobj - dataobject
%
%	LOCAL.C_READ('list') list all datasets that are locally available and indexed
%		To see variables for any data set see LOCAL.C_CAA_META
%
%	ok = LOCAL.C_READ('test') test if local data directory exists, output
%	is true or false.
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
% The data are loaded from local data directory. By default it is
% '/data/caalocal' but you can set it to different path by running in
% matlab:
%	> datastore('caa','localDataDirectory','/local/data/path');
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

if isempty(usingNasaPatchCdf) % check only once if using NASA cdf
  usingNasaPatchCdf=irf.check_if_using_nasa_cdf;
end

%% Defaults
returnDataFormat = 'mat'; % default matlab format
listIndexedDatasets = false;

% assign default caaDir value from datastore
caaDir = datastore('caa','localDataDirectory');
caaDirDefault = '/data/caalocal';  % value to suggest if caaDir not defined

%% Default index is empty, read in only those indees that are used
data=[];
if isempty(index)
  index=struct('dummy',[]);
end
%% Check inputs
if nargin == 0
  help local.c_read;
  return;
end
if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'list')
  listIndexedDatasets = true;
elseif nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'test')
  out = false;
  if exist(caaDir,'dir'), out = true; end
  return
elseif nargin>=2
  varName=varargin{1};
  tint=varargin{2};
  if ischar(tint)
    tint=irf_time(tint,'utc>tint');
  elseif isa(tint,'GenericTimeArray') && length(tint)==2
    tintTemp = tint.epochUnix;
    tint = tintTemp;
  end
end
if nargin ==3
  if ischar(varargin{3}) ...
      && any(strcmpi(varargin{3},{'dobj','caa','mat','ts'}))
    returnDataFormat = varargin{3};
  else
    irf.log('critical','output data format unknown');
    out=[];
    return;
  end
end
if nargin > 3
  irf.log('critical','max 3 arguments supported');
  return
end
%% Check if data directory exists
if isempty(caaDir) % not saved in datastore
  caaDir = input(['Input local caa directory [default:' ...
    caaDirDefault ']:'],'s');
  if isempty(caaDir)
    disp(['Using default data directory: ' caaDirDefault]);
    caaDir=caaDirDefault;
  end
  ok = input('Shall I save the directory location for future sessions [y/n]?','s');
  if strcmpi(ok,'y')
    datastore('caa','localDataDirectory',caaDir);
  end
end

%% Check if repository is there
if ~exist(caaDir,'dir')
  disp(['Local CAA data repository ' caaDir ' not available!']);
  disp('If you want to use other directory as default, execute:');
  disp('datastore(''caa'',''localDataDirectory'',''/your/data/directory'');');
  return;
end

%% Check if we have to only list and return
if listIndexedDatasets
  list_indexed_datasets;
  return;
end
%% Read in data
out=[];dataobject=[]; % default return empty output
specialCaseCis=0;
switch lower(varName)
  case {'r'}
    varToRead={'sc_r_xyz_gse__CL_SP_AUX','sc_dr1_xyz_gse__CL_SP_AUX',...
      'sc_dr2_xyz_gse__CL_SP_AUX','sc_dr3_xyz_gse__CL_SP_AUX','sc_dr4_xyz_gse__CL_SP_AUX'};
    ok=read_data;
    if ok && strcmpi(returnDataFormat,'mat')
      out.R=[data{1} double(data{2})];
      c_eval('out.R?=[data{1} double(data{2}+data{2+?})];')
    end
  case {'b'}
    varToRead={'B_vec_xyz_gse__C1_CP_FGM_FULL','B_vec_xyz_gse__C2_CP_FGM_FULL',...
      'B_vec_xyz_gse__C3_CP_FGM_FULL','B_vec_xyz_gse__C4_CP_FGM_FULL'};
    ok=read_data;
    if ok && strcmpi(returnDataFormat,'mat')
      c_eval('out.B?=[data{1} double(data{1+?})];')
    end
  case {'r1','r2','r3','r4'}
    varToRead={'sc_r_xyz_gse__CL_SP_AUX',['sc_dr' varName(2) '_xyz_gse__CL_SP_AUX']};
    ok=read_data;
    if ok && strcmpi(returnDataFormat,'mat')
      out=[data{1} double(data{2}+data{3})];
    end
  case {'dr1','dr2','dr3','dr4'}
    varToRead={['sc_dr' varName(3) '_xyz_gse__CL_SP_AUX']};
    ok=read_data;
    if ok && strcmpi(returnDataFormat,'mat')
      out=[data{1} double(data{2})];
    end
  otherwise
    irf.log('warning',['local.c_read() reading variable: ' varName]);
    if strfind(varName,'CIS'),specialCaseCis=1;end %#ok<STRIFCND>
    varToRead={varName};
    ok=read_data;
    if ok && iscell(data) && isscalar(data)
      out = data{1};
    elseif ok && strcmpi(returnDataFormat,'mat')
      if numel(data)==2 && numel(size(data{2}))==2
        out=[data{1} double(data{2})];
      else
        out=data;
      end
    elseif ok && (strcmpi(returnDataFormat,'dobj') || strcmpi(returnDataFormat,'caa'))
      out = data;
    end
end

%% Functions
  function status=read_data
    status = false; % default
    %% find index
    ii=strfind(varToRead{1},'__');
    if ii
      dataset=varToRead{1}(ii+2:end);
      datasetIndex = strrep(dataset,'CIS-','CIS_');
      datasetDir = [caaDir filesep datasetIndex];
      if ~isfield(index,datasetIndex) % index not yet read
        indexVarName = ['index_' datasetIndex];
        indexFileInfo=dirwhos(datasetDir,indexVarName);
        if numel(indexFileInfo)==0 % there is no index
          irf.log('critical',['There is no index file:' indexVarName]);
          irf.log('critical','Check that your localDataPath is correct, see help!');
          return;
        end
        s=dirload(datasetDir,indexVarName);
        index.(datasetIndex)=s.(indexVarName);
      end
      index=index.(datasetIndex);
      if isempty(index)
        irf.log('warning',['local.c_read: no data for dataset ' dataset]);
        return;
      end
    else
      irf.log('critical',['Do not know how to read variable: ' varToRead{1}]);
      return
    end
    %% find files within time interval
    istart=find(index.tend>tint(1),1);
    iend=find(index.tstart<tint(2),1,'last');
    irf.log('notice',['Dataset: ' dataset '. Index files: ' num2str(istart) '-' num2str(iend)]);

    if isempty(istart) || isempty(iend) || istart > iend
      return
    end
    %% read in records
    for iFile=istart:iend
      cdfFile=[caaDir filesep index.filename(iFile,:)];
      if specialCaseCis
        dataset=strrep(dataset,'CIS_','CIS-');
        varToRead=strrep(varToRead,'CIS_','CIS-');
      end
      % Get the correct CDF variables for vars>64 symbols
      for  iVar=1:numel(varToRead)
        if length(varToRead{iVar})>64
          varToRead{iVar}=[varToRead{iVar}(1:54) '...' varToRead{iVar}(end-6:end)];
        end
      end
      switch returnDataFormat
        case {'mat','ts'}
          irf.log('notice',['Reading: ' cdfFile]);
          %% check if epoch16
          cdfid=cdflib.open(cdfFile);
          useCdfepoch16=strcmpi(cdflib.inquireVar(cdfid,0).datatype,'cdf_epoch16');
          if useCdfepoch16
            irf.log('debug',['EPOCH16 time in cdf file:' cdfFile]);
            tName  = cdflib.getVarName(cdfid,0);
            tData = spdfcdfread(cdfFile,'CombineRecords',true,'KeepEpochAsIs',true,'Variables',{tName});
            if numel(size(tData)) == 3
              tcdfepoch=reshape(tData,size(tData,1),size(tData,3)); % spdfcdfread returns (Nsamples X 1 X 2) matrix
            else
              tcdfepoch = tData'; % spdfcdfread returns (2 x Nsamples) matrix
            end
            timeVector=irf_time(tcdfepoch,'cdfepoch16>epoch');
            tmpdata=cell(1,numel(varToRead));
            for iVar=1:numel(varToRead)
              tmpdata{iVar}=spdfcdfread(cdfFile,'CombineRecords',true,...
                'Variables',varToRead{iVar});
            end
            tmpdata = [{timeVector} tmpdata]; %#ok<AGROW>
          else
            % remove time variable as it is already read in
            ii=numel(varToRead);
            while ii
              if strcmp(varToRead{ii},cdflib.getVarName(cdfid,0))
                varToRead(ii)=[];
              end
              ii=ii-1;
            end
            % read data
            [tmpdata,~] = spdfcdfread(cdfFile,'ConvertEpochToDatenum',true,'CombineRecords',true,...
              'Variables', [{cdflib.getVarName(cdfid,0)},varToRead{:}]); % time and variable name TODO: get rid of conversion to datenum!
            tmpdata=fix_order_of_array_dimensions(tmpdata);
            if isnumeric(tmpdata), tmpdata={tmpdata}; end % make cell in case matrix returned
            timeVector = irf_time(tmpdata{1},'date>epoch');
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
          for iVar=1:numel(varToRead)
            fillVal=value_of_variable_attribute(cdfid,varToRead{iVar},'FILLVAL');
            tmpdata{iVar+1}(tmpdata{iVar+1}==fillVal)=NaN; % +1 because first cell is time and then comes variables
          end
          %% attach to result
          for j=1:numel(data)
            nDim=numel(size(tmpdata{j}));
            if nDim==2
              data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:));
            elseif nDim==3
              data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:,:));
            elseif nDim==4
              data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:,:,:));
            elseif nDim==5
              data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:,:,:,:));
            end
          end
          cdflib.close(cdfid);
        case {'caa','dobj'}
          if iFile==istart % start of interval, initiate dataobject
            dataobject=dataobj(cdfFile,'tint',tint);
          elseif iFile==iend % ends of interval
            data_temp=dataobj(cdfFile,'tint',tint);
          else
            data_temp=dataobj(cdfFile);
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
    if strcmp(returnDataFormat,'ts')
      Time = EpochUnix(data{1});
      dataTs = cell(numel(data)-1);
      dobj = dataobj(cdfFile);
      for iData = 1:numel(data)-1
        varTs = get_ts(dobj,varToRead{iData});
        nDataPoints = length(Time);
        varTsAll = varTs(ones(1,nDataPoints));
        varTsAll.time = Time;
        varTsAll.data = data{iData+1};
        dataTs{iData} = varTsAll;
      end
      data = dataTs;
    end
    if strcmp(returnDataFormat,'caa')
      data=get(dataobject,varToRead{1}); % currently only 1 variable request implemented
    elseif strcmp(returnDataFormat,'dobj')
      data=dataobject; % currently only 1 variable request implemented
    end
    status = true;
  end
  function value=value_of_variable_attribute(cdfid,varName,attrName)
    attrnum = cdflib.getAttrNum(cdfid,attrName);
    varnum = cdflib.getVarNum(cdfid,varName);
    value = cdflib.getAttrEntry(cdfid,attrnum,varnum);
  end
  function ok=list_indexed_datasets
    ok=false;
    tmp=dir(caaDir);
    iDir = [tmp(:).isdir]; % find directories
    dataSetArray = {tmp(iDir).name}';
    dataSetArray(ismember(dataSetArray,{'.','..'})) = []; % remove '.' and '..'
    for iDataSet=1:numel(dataSetArray)
      %% list files in data set directory
      dataSet=dataSetArray{iDataSet};
      dataSetDir = [caaDir filesep dataSet];
      if ~isempty(dirwhos(dataSetDir,['index_' dataSet]))
        disp(dataSet);
        ok = true;
      end
    end
  end
  function data=fix_order_of_array_dimensions(data)
    for iDimension=3:4
      indDatasets=find(cellfun(@(x) numel(size(x)),data(:))==iDimension); % find iDimension datasets
      for iDataset=1:numel(indDatasets)
        if iDimension==3
          data{indDatasets(iDataset)}=permute(data{indDatasets(iDataset)},[3 1 2]);
        elseif iDimension==4
          data{indDatasets(iDataset)}=permute(data{indDatasets(iDataset)},[4 3 1 2]);
        end
      end
    end
  end

end
