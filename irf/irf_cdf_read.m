function [varargout]=irf_cdf_read(cdf_file,var_name,flag)
%IRF_CDF_READ   Read CDF files
%
%  [varargout]=irf_cdf_read(cdf_file,var_name,flag)
%
% Usage:
% irf_cdf_read - interactively choose file and variable
% irf_cdf_read('*.cdf','*') - allow to choose file and variables
%
% reads variable var_name from cdf_file and converts it to column vector
% with first column being time in isdat epoch
% cdf_file - file name, can include wildcards
% var_name - '*', variable name or cell structure with variable names
% flag
%   'latest' - if several cdf_files then choose latest version
%              otherwise file name is chosen interactively
%
% if no data available output is []
%
% Example:
%   irf_cdf_read('*','*.cdf')
%
% See also SPDFCDFREAD
%

flag_latest=0;
persistent old_cdfread_call;

if nargin==3
  switch flag
    case 'latest'
      flag_latest=1;
  end
end

if nargin<=1, var_name='*';end
if nargin==0
  cdf_files=dir('*.cdf');
  switch numel(cdf_files)
    case 0
      disp('no cdf files specified');return
    case 1
      cdf_file=cdf_files.name;
      disp(['Using: ' cdf_file]);
    otherwise
      D=dir('*.cdf');
      for j=1:length(D)
        disp([num2str(j) '. ' D(j).name]);
      end
      disp('Choose cdf file');
      cdf_file_number=irf_ask('Cdf_file? [%]>','cdf_file_number',1);
      cdf_file=D(cdf_file_number).name;
  end
end

if strfind(cdf_file,'*')  %#ok<STRIFCND> % use wilcard '*' expansion
  ii=strfind(cdf_file,filesep);
  if ii, cdf_directory=cdf_file(1:max(ii)); else, cdf_directory=''; end
  ff=dir(cdf_file);
  switch size(ff,1)
    case 0
      irf_log('load','No cdf files')
      if nargout>0
        varargout = cell(1,nargout);
        for i=1:nargout, varargout(i) = {[]}; end
      end
      return
    case 1
      cdf_file=[cdf_directory ff(1).name];
    otherwise
      cdf_names={};
      for j=1:size(ff,1)
        cdf_names{end+1}=ff(j).name; %#ok<AGROW>
        disp([num2str(j) '. ' ff(j).name]);
      end
      if flag_latest
        cdf_names_sort=sort(cdf_names);
        cdf_file=[cdf_directory cdf_names_sort{end}];
      else
        cdf_file_n=irf_ask('cdf_file? [%]>','cdf_file_n',1);
        cdf_file=[cdf_directory ff(cdf_file_n).name];
      end
  end
end


irf_log('load',['cdf file: ' cdf_file]);

cdf_file_info = spdfcdfinfo(cdf_file);
variable_names = cdf_file_info.Variables(:,1);
epoch_column=0; % default there is no epoch variable
for j=1:size(variable_names,1)
  if any(strfind(variable_names{j,1},'Epoch')) || ...
      any(strfind(variable_names{j,1},'time_tags__')) || ...
      any(strfind(variable_names{j,1},'_time')) || ...
      any(strfind(variable_names{j,1},'Time__'))
    epoch_variable=variable_names(j,1);
    %disp(['epoch variable: ' epoch_variable{1}]);
    epoch_column=j;
    break;
  elseif strcmp(cdf_file_info.Variables(j,4),'epoch')
    epoch_variable=variable_names(j,1);
    epoch_column=j;
    break;
  end
end

if iscell(var_name)
elseif ischar(var_name) % one specifies the name of variable
  % in case string is '*' do interactive work with cdf file including
  % reading variables
  % get variable list that have associated time
  i_time_series_variable=0; % the counter of variables that depend on time
  if epoch_column ~= 0
    for j=2:size(cdf_file_info.Variables,1)
      if cdf_file_info.Variables{j,3}==cdf_file_info.Variables{epoch_column,3}
        if isempty(strfind(cdf_file_info.Variables{j,1},'Epoch')) && ...
            isempty(strfind(variable_names{j,1},'time_tags__')) && ...
            any(strfind(variable_names{j,1},'_time')) && ...
            isempty(strfind(variable_names{j,1},'Time__'))
          i_time_series_variable=i_time_series_variable+1;
          time_series_variables{i_time_series_variable}=cdf_file_info.Variables{j,1}; %#ok<AGROW>
        end
      end
    end
  end
  if strcmp(var_name,'*')
    flag_while=1; var_name=[];
    while flag_while
      disp('======== Options ========')
      disp('q) quit, dont return anything')
      disp('v) list all variables')
      disp('vv) list all variables and their values')
      disp('r) read variable in default format')
      disp('t) read time variables in irfu format and return')
      disp('fa) list the file variable attributes')
      disp('fav) view the file variable attributes')
      disp('ga) list global attributes')
      disp('gav) view global attributes')
      var_additional_item=irf_ask('Variable? [%]>','var_additional_item','v');
      switch var_additional_item
        case 'q'
          return;
        case 'v'
          cdf_file_info.Variables
        case 'vv'
          for jj=1:size(cdf_file_info.Variables,1)
            cdf_var=cdf_file_info.Variables{jj,1};
            disp('======================================')
            disp(cdf_var);
            disp('======================================')
            dd=spdfcdfread(cdf_file,'Variables',{cdf_var});
            if size(dd,1)>10
              disp([num2str(size(dd,1)) ' samples. 1st sample below.']);
              if iscell(dd), disp(dd{1}), else, disp(dd(1,:)), end
            else
              if iscell(dd), disp(dd{1}), else, disp(dd), end
            end
          end
        case 'r'
          var_to_read=irf_ask('Variable name? [%]>','var_to_read','');
          evalin('caller',[var_to_read '= spdfcdfread(''' cdf_file ''', ''VARIABLES'', ''' var_to_read ''');']);
        case 'fa'
          cdf_file_info.VariableAttributes
        case 'fav'
          ssf=fieldnames(cdf_file_info.VariableAttributes);
          for jj=1:length(ssf)
            disp('*************************************')
            disp(ssf{jj});
            disp('*************************************')
            eval(['disp(cdf_file_info.VariableAttributes.' ssf{jj} ')'])
          end
        case 'ga'
          cdf_file_info.GlobalAttributes
        case 'gav'
          ssf=fieldnames(cdf_file_info.GlobalAttributes);
          for jj=1:length(ssf)
            disp('*************************************')
            disp(ssf{jj});
            disp('*************************************')
            eval(['disp(cdf_file_info.GlobalAttributes.' ssf{jj} ')'])
          end
        case 't'
          if i_time_series_variable == 0
            irf_log('load','there are no time variables identified')
          elseif i_time_series_variable == 1
            var_name='all';
            irf_log('load',['var: ' time_series_variables{1}])
          else
            disp('=== Choose time dependant variable ===');
            disp('0) all time time dependant variables');
            for j=1:i_time_series_variable
              disp([num2str(j) ') ' time_series_variables{j}]);
            end
            var_item=irf_ask('Variable? [%]>','var_item',0);
            if var_item==0 % read all
              var_name='all';
            elseif var_item==-2
              return;
            elseif var_item==-1
            else
              var_name={''};
              for j=1:length(var_item)
                var_name(j)=time_series_variables(var_item(j));
              end
            end
          end
          flag_while=0;
        otherwise
          eval(var_additional_item);
      end
    end
  end
  if strcmp(var_name,'all')
    var_name=time_series_variables;
  end
else
  disp('var_name should be string or cell array'); whos var_name;
  return;
end

if ischar(var_name), var_name={var_name};end
variables=[epoch_variable; var_name(:)];
% if requesting only time variable, remove the duplicate
if numel(variables)==2 && strcmp(variables{1},variables{2}), variables(2)=[];end

if isempty(old_cdfread_call)
  try
    % New, faster call (requires Matlab R2008b or higher and the CDF patch from Goddard):
    DATA = spdfcdfread(cdf_file, 'VARIABLES', variables,'CombineRecords', true,'ConvertEpochToDatenum', true);
    old_cdfread_call=0;
  catch %#ok<CTCH>
    % Old, slow call:
    DATA = cdfread(cdf_file, 'VARIABLES', variables);
    disp('**************** Using the old and slow cdfread call.')
    disp('                 Consider upgrading Matlab (R2008b) or using the Goddard CDF patch.')
    old_cdfread_call=1;
  end
else
  if(old_cdfread_call), DATA = cdfread(cdf_file, 'VARIABLES', variables);
  else, DATA = spdfcdfread(cdf_file, 'VARIABLES', variables,'CombineRecords', true,'ConvertEpochToDatenum', true); end
end

if(old_cdfread_call)
  temp=struct([DATA{:,1}]);
  t = [temp.date];
  t = t(:);
  t = (t-62167219200000)/1000; %#ok<NASGU>
else
  if iscell(DATA)
    t = DATA{1};
  else
    t=DATA;
  end
  if t<1e6 % conversion succeeded (there is bug on THEMIS files)
    t = (t-1).* (24 * 3600 * 1000);
    t = (t-62167219200000)/1000; %#ok<NASGU>
  else % THEMIS case, dont do anything, is already in epoch
  end
end

for k=2:numel(variables)
  clear var;
  if(old_cdfread_call)
    for j=1:size(DATA,1)
      temp=DATA{j,k};temp=temp(:)';var(j,:)=temp; %#ok<AGROW>
    end
    var=[DATA{:,k}]';
  else
    var=DATA{k};
  end
  
  var(var<-1e30) = NaN;
  if isfield(cdf_file_info.VariableAttributes,'FILLVAL')
    for j = 1:size(cdf_file_info.VariableAttributes.FILLVAL,1)
      if strcmp(variables{k},cdf_file_info.VariableAttributes.FILLVAL{j,1})
        %disp([variables{k} ' : FILLVAL ' ...
        %		num2str(cdf_file_info.VariableAttributes.FILLVAL{j,2})])
        if ischar(cdf_file_info.VariableAttributes.FILLVAL{j,2})
          fill_val = str2double(cdf_file_info.VariableAttributes.FILLVAL{j,2});
        else
          fill_val = cdf_file_info.VariableAttributes.FILLVAL{j,2};
        end
        var( var == fill_val ) = NaN;
      end
    end
  end
  eval([variables{k} '=[t double(var)];' ]);
end

if nargout==0
  for k=2:numel(variables)
    assignin('caller',variables{k},eval(variables{k}));
  end
else
  varargout = cell(1,nargout);
  if numel(variables)==1 % return only time variable
    varargout{1}=t;
  else
    for i=1:nargout, varargout(i) = {eval(variables{i+1})}; end
  end
end
