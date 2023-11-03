function dobj = dataobj(varargin)
%DATAOBJ  constructor function for DATAOBJ object
%
% DATAOBJ(FILENAME)
%    Construct dataobj form file FILENAME. FILENAME can also contain
%    wildcards ('*').
%
% DATAOBJ(FILENAME,'tint',tint)
%       tint = [t1, t2]        - Limit dataobj to time interval.
%                                ISDAT epoch, i.e. the number of seconds since 1-Jan-1970.
%                                See epoch2iso.

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
persistent usingNasaPatchCdf

if isempty(usingNasaPatchCdf) % check only once if using NASA cdf
  usingNasaPatchCdf=irf.check_if_using_nasa_cdf;
end
shouldReadAllData = true;  % default read all data
isDataReturned    = true;  % default expects data to be returned
if     nargin==0, action='create_default_object';
elseif nargin==1, action='read_data_from_file';
else
  action='read_data_from_file';
  if nargin==2 || nargin==4
    if isnumeric(varargin{nargin}) ||...
        ischar(varargin{nargin}) && strcmpi(varargin{nargin},'KeepTT2000')
      irf.log('warning','KeepTT2000 is deprecated (always TRUE)')
    end
  end
  if nargin >=3
    if ischar(varargin{2}) && strcmp(varargin{2},'tint') && ...
        isnumeric(varargin{3}) && (length(varargin{3})==2)
      tint=varargin{3};
      shouldReadAllData=false;
    else
      error('IRF:dataobj:dataobj:invalidInput','invalid input params');
    end
  end
end

switch action
  case 'create_default_object'
    % if no input arguments, create a default object
    dobj.FileModDate = char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"));
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
            error('No cdf files specified');
          case 1
            cdf_file = [directory_name cdf_files.name];
          otherwise
            for j=1:numel(cdf_files)
              disp([num2str(j) '. ' cdf_files(j).name]);
            end
            j = irf_ask('Choose cdf-file? [%]>','cdf_file',1);
            cdf_file = [directory_name cdf_files(j).name];
        end
        clear cdf_files
      else, cdf_file = varargin{1};
      end
      if ~exist(cdf_file,'file')
        error(['File ' cdf_file ' does not exist']);
      end
      %% read in file
      irf.log('notice',['Reading: ' cdf_file]);
      % initialize data object
      if usingNasaPatchCdf
        timeVariable = [];
        [data,info] = spdfcdfread(cdf_file,'CombineRecords',true,'KeepEpochAsIs',true);

        % Convert epoch/epoch16/tt2000 to ISDAT epoch
        isCdfEpochVariableArray=cellfun(@(x) strcmpi(x,'epoch'), info.Variables(:,4));
        if any(isCdfEpochVariableArray)
          iVar = find(isCdfEpochVariableArray);
          for i=1:length(iVar)
            if is_virtual(iVar(i))
              virtFunc = get_key('FUNCT',iVar(i));
              switch lower(virtFunc{:})
                case 'comp_themis_epoch'
                  depVarName = get_key('COMPONENT_1',iVar(i));
                  idx = get_var_idx(depVarName{:});
                  data{iVar(i)} = data{idx}(:);
                  data{iVar(i)}(data{iVar(i)}==0) = NaN; % fillvalue timeline
                  info.Variables{iVar(i),3} = length(data{iVar(i)});
                otherwise
                  errStr = sprintf('Function ''%s'' not implemented',virtFunc);
                  irf.log('error',errStr);
                  error('IRF:dataobj:dataobj:functionNotImplemented',errStr); %#ok<SPERR>
              end
            else
              data{iVar(i)}(data{iVar(i)}==0) = NaN; % fillvalue timeline
              dTmp = data{iVar(i)};
              for iCol=1:size(dTmp,2)
                dTmp(:,iCol) = irf_time(dTmp(:,iCol),'cdfepoch>epoch');
              end
              data{iVar(i)} = dTmp; clear dTmp
              % bug fix for spdfcdfread (time comes out as row vector)
              if size(data{iVar(i)},1)~=info.Variables{iVar(i),3}
                data{iVar(i)} = data{iVar(i)}';
              end
            end
            timeVariable = info.Variables{iVar(i),1};
            update_variable_attributes_time;
          end
        end
        isCdfEpoch16VariableArray=cellfun(@(x) strcmpi(x,'epoch16'), info.Variables(:,4));
        if any(isCdfEpoch16VariableArray)
          iVar = find(isCdfEpoch16VariableArray);
          for i=1:length(iVar)
            if is_virtual(iVar(i))
              virtFunc = get_key('FUNCT',iVar(i));
              switch lower(virtFunc{:})
                case 'comp_themis_epoch16'
                  depVarName = get_key('COMPONENT_1',iVar(i));
                  idx = get_var_idx(depVarName{:});
                  data{iVar(i)} = data{idx}(:);
                  data{iVar(i)}(data{iVar(i)}==0) = NaN; % fillvalue timeline
                  info.Variables{iVar(i),3} = length(data{iVar(i)});
                otherwise
                  errStr = sprintf('Function ''%s'' not implemented',virtFunc);
                  irf.log('error',errStr);
                  error('IRF:dataobj:dataobj:functionNotImplemented',errStr); %#ok<SPERR>
              end
              timeVariable = info.Variables{iVar(i),1};
              update_variable_attributes_time;
            else
              data{iVar(i)} = convert_cdfepoch16_string_to_isdat_epoch(data{iVar(i)});
            end
          end
        end
        isCdfEpochTT2000VariableArray=cellfun(@(x) strcmpi(x,'tt2000'), info.Variables(:,4));
        if (any(isCdfEpochTT2000VariableArray))
          iVar = find(isCdfEpochTT2000VariableArray);
          for i=1:length(iVar)
            if is_virtual(iVar(i))
              virtFunc = get_key('FUNCT',iVar(i));
              errStr = sprintf('Function ''%s'' not implemented',virtFunc);
              irf.log('error',errStr);
              error('IRF:dataobj:dataobj:functionNotImplemented',errStr); %#ok<SPERR>
            end
          end
        end

        fix_order_of_array_dimensions;

        if ~shouldReadAllData
          nVariables = size(info.Variables,1);
          records = cell(nVariables,1); recsTmp = {};
          for i=1:nVariables
            % Time dependent variables
            timeVarName = get_key('DEPEND_0',i);
            if isempty(timeVarName), continue, end
            timeVarName = timeVarName{:};
            if ~isempty(recsTmp)
              isRecArray=cellfun(@(x) strcmpi(x,timeVarName), recsTmp(:,1));
              idx = find(isRecArray==1);
              if ~isempty(idx), records(i) = recsTmp(idx,3); continue, end
            end
            iTimeVar = get_var_idx(timeVarName);
            timeline = data{iTimeVar};
            if strcmpi(info.Variables(iTimeVar,4),'tt2000')
              tintTmp(1) = spdfparsett2000(epoch2iso(tint(1)));
              tintTmp(2) = spdfparsett2000(epoch2iso(tint(2)));
            else, tintTmp = tint;
            end
            records{i} = (timeline >= tintTmp(1)) & (timeline < tintTmp(2));
            recsTmp = [recsTmp; {timeVarName,iTimeVar,records{i}}]; %#ok<AGROW>
            % Theat DELTA_PLUS/DELTA_MINUS for time variables
            for deltaParam={'DELTA_PLUS','DELTA_MINUS'}
              deltaVar = get_key(deltaParam{:},iTimeVar);
              if isempty(deltaVar) || ~ischar(deltaVar{:}), continue, end
              deltaVar = deltaVar{:}; iDelta = get_var_idx(deltaVar);
              recsTmp = [recsTmp; {deltaVar,iDelta,records{i}}]; %#ok<AGROW>
            end
          end
          % Time variables
          if ~isempty(recsTmp)
            for i = 1:size(recsTmp,1), records{recsTmp{i,2}} = recsTmp{i,3}; end
          end
        end
      else % i.e. NOT using NASA cdf patch
        % get basic info
        info   = cdfinfo(cdf_file);
        %check for epoch tt2000
        usingTT2000=strcmpi('tt2000',info.Variables{1,4});
        if usingTT2000
          se = 'NASA CDF patch is required to read files containing TT2000 variables';
          irf.log('critical',se);
          error('IRF:dataobj:dataobj:unsupported_data',se);
        end
        % check if cdfepoch16
        usingCdfepoch16=strcmpi('epoch16',info.Variables{1,4});
        % read in file
        if usingCdfepoch16
          irf.log('notice',['EPOCH16 time in cdf file:' cdf_file]);
          shouldReadAllData=1; % read all data
          isCdfEpoch16VariableArray=cellfun(@(x) strcmpi(x,'epoch16'), info.Variables(:,4));
          data=cell(1,size(variableNames,1));
          data(~isCdfEpoch16VariableArray) = cdfread(cdf_file,'variables',variableNames(~isCdfEpoch16VariableArray),'CombineRecords',true);
          iCdf16Variable = find(isCdfEpoch16VariableArray);
          for i=1:length(iCdf16Variable)
            numrecs = info.Variables{iCdf16Variable(i),3};
            % get time axis
            tc=zeros(2,numrecs);
            cdfid   = cdflib.open(cdf_file);
            for jj=1:numrecs
              tc(:,jj) = cdflib.getVarRecordData(cdfid,iCdf16Variable(i)-1,jj-1);
            end
            cdflib.close(cdfid);
            data(iCdf16Variable(i))={irf_time(tc','cdfepoch16>epoch')};
          end
        else
          [data,info] = cdfread(cdf_file,'ConvertEpochToDatenum',true,'CombineRecords',true);
          data{1} = irf_time([data{:,1}],'date>epoch');
        end
        if ~shouldReadAllData
          records=find((data{1} >= tint(1)) & (data{1} <= tint(2)));
        end
        timeVariable=info.Variables{1,1};
        update_variable_attributes_cdfepoch();
      end
      %% check if number of records to read is zero
      if ~shouldReadAllData
        if iscell(records)
          isDataReturned = false;
          for i=1:length(records)
            if ~isempty(records{i}) && sum(records{i}) > 0
              isDataReturned = true; break
            end
          end
        elseif sum(records)==0, isDataReturned = false;
        end
        if ~isDataReturned
          irf.log('warning','No data within specified time interval');
        end
      end
      %% construct data object
      dobj.FileModDate		= info.FileModDate;
      dobj.VariableAttributes = info.VariableAttributes;
      dobj.GlobalAttributes	= info.GlobalAttributes;
      dobj.Variables			= info.Variables;

      % check if there is time data
      if ~(any(strcmpi(info.Variables(:,4),'epoch')==1) || ...
          any(strcmpi(info.Variables(:,4),'epoch16')==1) || ...
          any(strcmpi(info.Variables(:,4),'tt2000')==1))
        irf.log('warning','No time in CDF file.');
      end
      nVariables = size(info.Variables,1);
      dobj.vars = cell(nVariables,2);
      if nVariables>0
        dobj.vars(:,1) = info.Variables(:,1); % new variables
        dobj.vars(:,2) = info.Variables(:,1); % original variables
        for v=1:nVariables
          dobj.vars{v,1}=variable_mat_name(dobj.vars{v,2});
          varName = dobj.vars{v,1};
          dataAllRecords = data{v};
          if shouldReadAllData || ...% return all data
              (usingNasaPatchCdf && strcmpi(info.Variables{v,5}(1),'F'))% fixed variable with NASA cdf patch (matlab cdfread return fixed variable as time series)
            dobj.data.(varName).data = dataAllRecords;
            dobj.data.(varName).nrec = info.Variables{v,3};
          else
            if iscell(records), recsTmp = records{v}; else, recsTmp = records; end
            nDim=numel(size(dataAllRecords));
            switch nDim
              case 2
                dataSelectedRecords=dataAllRecords(recsTmp,:);
              case 3
                dataSelectedRecords=dataAllRecords(recsTmp,:,:);
              case 4
                dataSelectedRecords=dataAllRecords(recsTmp,:,:,:);
              case 5
                dataSelectedRecords=dataAllRecords(recsTmp,:,:,:,:);
            end
            dobj.data.(varName).data = dataSelectedRecords;
            dobj.data.(varName).nrec = sum(recsTmp);
            dobj.Variables{v,3}      = sum(recsTmp);
          end
          dobj.data.(varName).dim      = info.Variables{v,2};
          dobj.data.(varName).variance = info.Variables{v,5};
          dobj.data.(varName).sparsity = info.Variables{v,6};
          typeTmp = info.Variables{v,4};
          % Update time types to catch up with conversion to unix epoch
          if strcmpi(typeTmp,'epoch16'), typeTmp = 'epoch'; end
          dobj.data.(varName).type     = typeTmp;
        end
      end
      dobj = class(dobj,'dataobj');
    else
      error('Wrong argument type');
    end
  otherwise
    error('Wrong number of input arguments');
end % Main function

%% Help functions
  function idx = get_var_idx(varName)
    isVarArray=cellfun(@(x) strcmpi(x,varName), info.Variables(:,1));
    idx = find(isVarArray==1);
    if length(idx) > 1
      strErr = sprintf('Multiple entries for variable ''%s'' in field ''Variables''',varName);
      irf.log('error',strErr);
      error('IRF:dataobj:get_var_idx:multipleEntries',strErr); %#ok<SPERR>
    end
  end

  function key = get_key(field,iVar)
    key = '';
    if ~isfield(info.VariableAttributes,field)
      strErr = sprintf('VariableAttributes does not have a field ''%s''',field);
      irf.log('error',strErr);
      error('IRF:dataobj:get_key:noField',strErr); %#ok<SPERR>
    end
    varName = info.Variables(iVar,1);
    isVarArray=cellfun(@(x) strcmpi(x,varName), info.VariableAttributes.(field)(:,1));
    idxVar = find(isVarArray==1);
    if isempty(idxVar), return, end
    if length(idxVar) > 1
      strErr = sprintf('Multiple entries for variable ''%s'' in field ''%s''',varName,field);
      irf.log('error',strErr);
      error('IRF:dataobj:get_key:multipleEntries',strErr); %#ok<SPERR>
    end
    key = info.VariableAttributes.(field)(idxVar,2);
  end

  function res = is_virtual(iVar)
    res = false;
    if ~isfield(info.VariableAttributes,'VIRTUAL'), return, end
    isVarArray=cellfun(@(x) strcmpi(x,info.Variables(iVar,1)), info.VariableAttributes.VIRTUAL(:,1));
    if any(isVarArray)
      if strcmpi(info.VariableAttributes.VIRTUAL(isVarArray==1,2),'true')
        res = true;
      end
    end
  end

  function fix_order_of_array_dimensions
    %     oldPermute=true;
    %     if oldPermute
    for iDimension=3:4   %#ok<UNRCH>
      indDatasets=find(cellfun(@(x) numel(size(x)),data(:))==iDimension); % find iDimension datasets
      for iDataset=1:numel(indDatasets)
        if iDimension==3
          data{indDatasets(iDataset)}=permute(data{indDatasets(iDataset)},[3 1 2]);
        elseif iDimension==4
          data{indDatasets(iDataset)}=permute(data{indDatasets(iDataset)},[4 3 1 2]);
        end
      end
    end
    %     else
    %       % NEW METHOD; TO BE TESTED
    %       for iDimension=3:4   %#ok<UNRCH>
    %         % Check if dimensions match expected (ignoring Depend_0) and that
    %         % it is not a single vector (which could be "Nx1").
    %         indDatasets = find(cellfun(@(x) all(x(:) > 1) && numel(x)==iDimension-1, info.Variables(:,2)));
    %         for iDataset=1:numel(indDatasets)
    %           if iDimension==3
    %             % permutate it to N x m x n, where N corrsponds to record (ie Depend_0).
    %             data{indDatasets(iDataset)} = permute(data{indDatasets(iDataset)}, [3 1 2]);
    %           elseif iDimension==4
    %             % Check if it should be simple cyclic permutation (MMS), data
    %             % was read into "recSize x nRec". Also possibly a matrix of
    %             % size "recSize" if it only contained one single record.
    %             % Or if the data should be more complexly permuted (Cluster),
    %             % ie where the recSize does not match the data.
    %             recSize = info.Variables{indDatasets(iDataset), 2}; % Expected record size, (m x n x p)
    %             nRec = info.Variables{indDatasets(iDataset), 3}; % Numer of records, (N)
    %             dataSize = size(data{indDatasets(iDataset)}); % Size of data, as it was read
    %             if isequal([recSize, nRec], dataSize) || (nRec==1 && isequal(recSize, dataSize) )
    %               % Simple cyclic permutation to get N x m x n x p
    %               data{indDatasets(iDataset)} = permute(data{indDatasets(iDataset)}, [4 1 2 3]);
    %             else
    %               % Re-order Cluster data, from (n x p x m x N) to (N x m x n x p).
    %               data{indDatasets(iDataset)} = permute(data{indDatasets(iDataset)}, [4 3 1 2]);
    %             end
    %           end
    %         end
    %       end
    %       % END OF NEW METHOD
    %     end
  end

  function update_variable_attributes_time % nested function
    isFieldUnits        = isfield(info.VariableAttributes,'UNITS');
    isFieldSIConversion = isfield(info.VariableAttributes,'SI_CONVERSION');
    isFieldDeltaPlus    = isfield(info.VariableAttributes,'DELTA_PLUS');
    isFieldDeltaMinus   = isfield(info.VariableAttributes,'DELTA_MINUS');
    if isFieldUnits
      iattr=find(strcmpi(info.VariableAttributes.UNITS(:,1),timeVariable));
      % change from ms to s UNITS of epoch if present
      if iattr, info.VariableAttributes.UNITS(iattr,2)={'s'};
      else
        info.VariableAttributes.UNITS(end+1,1)={timeVariable};
        info.VariableAttributes.UNITS(end,2)={'s'};
      end
    end
    if isFieldSIConversion
      iattr=find(strcmpi(info.VariableAttributes.SI_CONVERSION(:,1),timeVariable));
      % change from ms to s SI_CONVERSION of epoch if present
      if iattr, info.VariableAttributes.SI_CONVERSION(iattr,2)={'1.0>s'};
      else
        info.VariableAttributes.SI_CONVERSION(end+1,1)={timeVariable};
        info.VariableAttributes.SI_CONVERSION(end,2)={'1.0>s'};
      end
    end
    varType = info.Variables{cellfun(@(x) strcmpi(x,timeVariable),info.Variables(:,1)),4};
    switch varType
      case 'epoch', factor = 1e3;    % Original time in ms
      case 'epoch16', factor = 1e12; % Original time in ps
        % XXX: FIXME - implement real functionality when we get real files
      case 'tt2000', factor = 1e9;   % Original time in ns
      otherwise
    end
    % XXX: FIXME - update for VIRTUAL variables
    if isFieldDeltaPlus
      iattr=find(strcmpi(info.VariableAttributes.DELTA_PLUS(:,1),timeVariable)); % to convert DELTA_PLUS
      if any(iattr) && isnumeric(info.VariableAttributes.DELTA_PLUS{iattr,2})
        info.VariableAttributes.DELTA_PLUS{iattr,2}=info.VariableAttributes.DELTA_PLUS{iattr,2}/factor;
      end
    end
    if isFieldDeltaMinus
      iattr=find(strcmpi(info.VariableAttributes.DELTA_MINUS(:,1),timeVariable)); % to convert DELTA_PLUS
      if any(iattr) && isnumeric(info.VariableAttributes.DELTA_MINUS{iattr,2})
        info.VariableAttributes.DELTA_MINUS{iattr,2}=info.VariableAttributes.DELTA_MINUS{iattr,2}/factor;
      end
    end
  end

  function t=convert_cdfepoch16_string_to_isdat_epoch(in) % nested function
    % the following if is because of the bug in CAA cdf files having EPOCH16
    % sometimes time variable has dimension zero and sometimes one
    % TODO: report bug to CAA team and if needed cdf team
    if numel(size(in)) == 3
      tcdfepoch=reshape(in,size(in,1),size(in,3)); % cdfread returns (Nsamples X 1 X 2) matrix
    else
      tcdfepoch = in; % cdfread returns (Nsamples X 2) matrix
    end
    t=irf_time(tcdfepoch,'cdfepoch16>epoch');
  end
end % end of main functions

