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
      dir('*.cdf');
      disp('Choose cdf file');
      cdf_file=input('Cdf_file? >','s');
  end
end

if findstr(cdf_file,'*'), % use wilcard '*' expansion
  ii=findstr(cdf_file,'/');
  if ii, cdf_directory=cdf_file(1:max(ii)); else cdf_directory=''; end
  ff=dir(cdf_file);
  switch size(ff,1)
  case 0
   c_log('load','No cdf files')
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
    disp(ff(j).name);
   end
   if flag_latest
      cdf_names_sort=sort(cdf_names);
      cdf_file=[cdf_directory cdf_names_sort{end}];
   else
     cdf_file=input('cdf_file? >','s');
     cdf_file=[cdf_directory cdf_file];
   end
 end
end


c_log('load',['cdf file: ' cdf_file]);

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
elseif ischar(var_name)
  if strcmp(var_name,'*')
    disp('=== Choose variable ===');
    inf=cdfinfo(cdf_file);
    for j=2:size(inf.Variables,1)
      if inf.Variables{j,3}==inf.Variables{epoch_column,3}
        disp(inf.Variables{j,1});
      end
    end
    var_name=input('Variable? >','s');
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
