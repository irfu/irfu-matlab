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
persistent flag_using_nasa_patch_cdfread

if isempty(flag_using_nasa_patch_cdfread) % check which cdfread is used, NASA may give errors
    flag_using_nasa_patch_cdfread=0; % assuming as default that matlab cdfread is used
    fid=fopen(which('cdfread'));
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if strfind(tline,'Mike Liu')
            fprintf('\n\n\n');
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            disp(' You are using NASA cdfread patch!')
            disp(' This may give errors reading in multidimensional data sets!')
            disp(' Also option ''tint'' in routine databoj is disabled.');
            disp(' We suggest you to use the MATLAB cdfread!');
            disp(' To use MATLAB cdfread please remove path to NASA cdfread patch.');
            disp(' You can execute and then continue:');
            a=which('cdfread');
            ai=strfind(a,'/');
            disp(['> rmpath ' a(1:ai(end))]);
            disp('> clear databoj');
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            fprintf('\n\n\n');
            flag_using_nasa_patch_cdfread=1;
            break;
        end
    end
    fclose(fid);
end

flag_read_all_data=1; % default read all data
if nargin==0, action='create_default_object'; end
if nargin==1, action='read_data_from_file'; end
if nargin==3 && ...
    ischar(varargin{2}) && strcmp(varargin{2},'tint') && ...
    isnumeric(varargin{3}) && (length(varargin{3})==2),
  tint=varargin{3};
  if ~flag_using_nasa_patch_cdfread,
      irf_log('fcal',['returnig time interval limited data!' irf_time(tint,'tint2iso')])
      flag_read_all_data=0;
  end
  action='read_data_from_file';
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
            cdf_file = [directory_name cdf_files(j).name];
        end
        clear cdf_files
      else cdf_file = varargin{1};
      end
      
      if ~exist(cdf_file,'file')
        error(['file ' cdf_file ' does not exist'])
      end
      
      irf_log('dsrc',['Reading: ' cdf_file]);
      %% check if epoch16 file
      
      cdfid   = cdflib.open(cdf_file);
      inq=cdflib.inquireVar(cdfid,0);
      if strcmpi(inq.datatype,'cdf_epoch16')
          flag_using_cdfepoch16=1;
      else
          flag_using_cdfepoch16=0;
      end
      % leave open cdf file
      
      %% read in file 
      if flag_using_cdfepoch16, 
          irf_log('dsrc',['EPOCH16 time in cdf file:' cdf_file]);
          flag_read_all_data=1; % read all data
          info = cdflib.inquire(cdfid);
          vars=cell(info.numVars-1,1);
          vars_i16 = ones(size(1:info.numVars)); % array indicating which of the variables are EPOCH16
          for jj=1:info.numVars
              vars{jj}=cdflib.getVarName(cdfid,jj-1);
              inq=cdflib.inquireVar(cdfid,jj-1);
              if strcmpi(inq.datatype,'cdf_epoch16')
                  vars_i16(jj)=0;
              end
          end
          data=cell(1,info.numVars);
          data(vars_i16==1) = cdfread(cdf_file,'variables',vars(vars_i16==1),'CombineRecords',true);
          info=cdfinfo(cdf_file);
          ii = find(vars_i16==0);
          numrecs = cdflib.getVarAllocRecords(cdfid,0);
          for i=1:length(ii)
              % get time axis
              tc=zeros(2,numrecs);
              for jj=1:numrecs,
                  tc(:,jj) = cdflib.getVarRecordData(cdfid,ii(i)-1,jj-1);
              end
              data(ii(i))={tc'};
          end
      else
      [data,info] = cdfread(cdf_file,...
        'ConvertEpochToDatenum',true,...
        'CombineRecords',true);
      end
      cdflib.close(cdfid);
      if flag_read_all_data==0, % check which records to return later
        info=cdfinfo(cdf_file);
        timevar=info.Variables{strcmpi(info.Variables(:,4),'epoch')==1,1};
        timeline = irf_time(cdfread(cdf_file,'Variable',{timevar},'ConvertEpochToDatenum',true,'CombineRecords',true),'date2epoch');
        records_within_interval=find((timeline > tint(1)) & (timeline < tint(2)));
      end
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
          dobj.vars{v,1}(strfind(dobj.vars{v,1},'-')) = '_';
          % Remove training dots
          while (dobj.vars{v,1}(end) == '.')
            dobj.vars{v,1}(end) = [];
          end
          % Take care of '...'
          d3 = strfind(dobj.vars{v,1},'...');
          if d3, dobj.vars{v,1}( d3 + (1:2) ) = []; end
          % Replace dots with underscores
          dobj.vars{v,1}(strfind(dobj.vars{v,1},'.')) = '_';
          % Add "x" if the varible name starts with a number
          if ~isletter(dobj.vars{v,1}(1)),
            dobj.vars{v,1}=['x' dobj.vars{v,1}];
          end
          % Take care of names longer than 63 symbols (Matlab limit)
          if length(dobj.vars{v,1})>63
            dobj.vars{v,1} = dobj.vars{v,1}(1:63);
            disp(['orig var : ' dobj.vars{v,2}])
            disp(['new var  : ' dobj.vars{v,1}])
          end
          if flag_read_all_data, % return all data
            dobj.data.(dobj.vars{v,1}).data = [data{:,v}];
            dobj.data.(dobj.vars{v,1}).nrec = info.Variables{v,3};
          else
            data_all_records=[data{:,v}];
            if numel(size(data_all_records))==2,
              data_records_within_interval=data_all_records(records_within_interval,:);
            elseif numel(size(data_all_records))==3,
              data_records_within_interval=data_all_records(records_within_interval,:,:);
            elseif numel(size(data_all_records))==4,
              data_records_within_interval=data_all_records(records_within_interval,:,:,:);
            elseif numel(size(data_all_records))==5,
              data_records_within_interval=data_all_records(records_within_interval,:,:,:,:);
            end
            dobj.data.(dobj.vars{v,1}).data = data_records_within_interval;
            dobj.data.(dobj.vars{v,1}).nrec = numel(records_within_interval);
          end
          dobj.data.(dobj.vars{v,1}).dim = info.Variables{v,2};
          dobj.data.(dobj.vars{v,1}).type = info.Variables{v,4};
          dobj.data.(dobj.vars{v,1}).variance = info.Variables{v,5};
          dobj.data.(dobj.vars{v,1}).sparsity = info.Variables{v,6};
          %Convert to isdat epoch
          if strcmp(dobj.data.(dobj.vars{v,1}).type,'epoch')
            dobj.data.(dobj.vars{v,1}).data = irf_time(dobj.data.(dobj.vars{v,1}).data,'date2epoch');
          elseif strcmp(dobj.data.(dobj.vars{v,1}).type,'epoch16')
            dobj.data.(dobj.vars{v,1}).data = irf_time(dobj.data.(dobj.vars{v,1}).data,'cdfepoch162epoch');
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
