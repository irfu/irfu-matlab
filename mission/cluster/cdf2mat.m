function status=cdf2mat(cdf_wildcard,mat_file)
% cdf2mat convert cdf filel to matlab file (only vectors and scalar time series)
% function status=cdf2mat(cdf_file,mat_file);
%

if nargin<1, cdf_wildcard='*.cdf';end
if nargin<2, flag_mat_file_name=0; end

cdf_files=dir(cdf_wildcard);
switch numel(cdf_files)
  case 0
    c_log('load','no cdf files found');return
  case 1
    cdf_file{1}=cdf_files.name;
    irf_log('load',['converting cdf file: ' cdf_file{1}]);
    n_cdf_files=1;
  otherwise
    n_cdf_files=size(cdf_files,1);
    for i_file=1:n_cdf_files
      ttt=cdf_files(i_file);
      cdf_file{i_file}=ttt.name;
    end
    irf_log('load',['converting ' num2str(n_cdf_files) ' cdf files.']);
end

for i_file=1:n_cdf_files
  irf_log('load',['Loading ' num2str(i_file) '. cdf file>' cdf_file{i_file}]);
  cdf_file_info=spdfcdfinfo(cdf_file{i_file});
  variable_names=cdf_file_info.Variables(:,1);
  % keep only variables that have time axis

  for j=1:size(variable_names,1)
    if strfind(variable_names{j,1},'Epoch')
      epoch_variable=variable_names(j,1);
      % disp(['epoch variable: ' epoch_variable{1}]);
      epoch_column=j;
    end
  end
  n_time_vectors=0;
  for j=1:size(variable_names,1)
    if cdf_file_info.Variables{j,3}==cdf_file_info.Variables{epoch_column,3}
      n_time_vectors=n_time_vectors+1;
      variables_time_vectors(n_time_vectors,:)=variable_names(j,:);
      irf_log('proc',['variable: ' variable_names{j,1}]);
    end
  end

  [DATA, INFO] = spdfcdfread(cdf_file{i_file}, 'VARIABLES', variables_time_vectors);

  temp=struct([DATA{:,1}]);
  t=[temp.date];t=t(:);
  t=(t-62167219200000)/1000;

  variable_list='';
  for k=2:numel(variables_time_vectors)
    clear var;
    for j=1:size(DATA,1)
      temp=DATA{j,k};temp=temp(:)';var(j,:)=temp;
    end
    var=[DATA{:,k}]';
    var(var<-1e30)=NaN;
    eval([variables_time_vectors{k} '=[t double(var)];' ]);
    variable_list=[variable_list ' ' variables_time_vectors{k}];
  end
  irf_log('save',['Saving to matlab file: ' mat_file]);
  eval(['save maux ' variable_list]);
end
