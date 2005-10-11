function [varargout]=irf_cdf_read(cdf_file,var_name,flag)
%IRF_CDF_READ   Read CDF files
% function [varargout]=irf_cdf_read(cdf_file,var_name,flag)
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
% $Id$

flag_latest=0;

if nargin==3,
  switch flag
  case 'latest'
    flag_latest=1;
  end
end

if nargin<=1, var_name='*';end
if nargin==0,
  cdf_files=dir('*.cdf');
  switch prod(size(cdf_files))
    case 0
      disp('no cdf files specified');return
    case 1
      cdf_file=cdf_files.name;
      disp(['Using: ' cdf_file]);
    otherwise
      D=dir('*.cdf');
      for j=1:length(D),
        disp([num2str(j) '. ' D(j).name]);
      end
      disp('Choose cdf file');
      cdf_file=irf_ask('Cdf_file? [%]>','cdf_file',1);
  end
end

if findstr(cdf_file,'*'), % use wilcard '*' expansion
  ii=findstr(cdf_file,'/');
  if ii, cdf_directory=cdf_file(1:max(ii)); else cdf_directory=''; end
  ff=dir(cdf_file);
  switch size(ff,1)
  case 0
   irf_log('load','No cdf files')
   if nargout>0,
    for i=1:nargout, varargout(i) = {[]}; end
   end
   return
  case 1
   cdf_file=[cdf_directory ff(1).name];
  otherwise
   cdf_names={};
   for j=1:size(ff,1),
    cdf_names{end+1}=ff(j).name;
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

cdf_file_info=cdfinfo(cdf_file);
variable_names=cdf_file_info.Variables(:,1);
for j=1:size(variable_names,1),
 if findstr(variable_names{j,1},'Epoch')
   epoch_variable=variable_names(j,1);
   % disp(['epoch variable: ' epoch_variable{1}]);
   epoch_column=j;
 end
end

if iscell(var_name),
elseif ischar(var_name) % one specifies the name of variable 
  % get variable list that have associated time 
  inf=cdfinfo(cdf_file); 
  i_time_series_variable=0; % the counter of variables that depend on time
  for j=2:size(inf.Variables,1)
    if inf.Variables{j,3}==inf.Variables{epoch_column,3}
      if isempty(findstr(inf.Variables{j,1},'Epoch')), % exclude epoch from variables
        i_time_series_variable=i_time_series_variable+1;
        time_series_variables{i_time_series_variable}=inf.Variables{j,1};
      end
    end
  end
  % in case string is '*' show all possibilities and allow to choose
  if strcmp(var_name,'*'),
      disp('=== Choose variable ===');
      disp('0) all variables');
      for j=1:i_time_series_variable,
        disp([num2str(j) ') ' time_series_variables{j}]);
      end
      var_item=irf_ask('Variable? [%]>','var_item',2);
      if var_item==0, % read all
        var_name='all';
      else
        var_name={''};
        for j=1:length(var_item),
          var_name(j)=time_series_variables(var_item(j));
        end
      end
  end
  if strcmp(var_name,'all'),
    var_name=time_series_variables;
  end
else
  disp('var_name should be string or cell array'); whos var_name;
  return;
end

if ischar(var_name), var_name={var_name};end

variables=[epoch_variable var_name];

[DATA, INFO] = cdfread(cdf_file, 'VARIABLES', variables);

temp=struct([DATA{:,1}]);
t=[temp.date];t=t(:);
t=(t-62167219200000)/1000;

for k=2:prod(size(variables))
  clear var;
  for j=1:size(DATA,1)
    temp=DATA{j,k};temp=temp(:)';var(j,:)=temp;
  end
  var=[DATA{:,k}]';
  i=find(var<-1e30);var(i)=NaN;
  eval([variables{k} '=[t double(var)];' ]);
end

if nargout==0,
 for k=2:prod(size(variables))
   assignin('caller',variables{k},eval(variables{k}));
 end
else
  for i=1:nargout, varargout(i) = {eval(variables{i+1})}; end
end
