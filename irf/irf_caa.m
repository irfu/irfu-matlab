%IRF_CAA   Work with CAA files
%
% See also SPDFCDFREAD
%

flag_show_menu=1;
flag_variable_menu=1;
flag_show_data_menu=1;
while 1
  if flag_show_menu
    disp('************************************************************************************')
    disp('********** IRF CAA *********')
    disp('m) menu')
    disp('q) quit')
    disp('a) load all caa files in cdf format')
    disp('f) load files')
    disp('l) list loaded files')
    disp('c) choose from loaded files')
  end
  flag_show_menu=0;
  var_menu=irf_ask('irf_caa>','var_menu','');
  switch var_menu
    case 'a'
      caa_load;
    case 'q'
      return;
    case 'm'
      flag_show_menu=1;
    case 'f'
      i_caa_files=[];
      dirs = dir;
      old_pwd = pwd;
      for j = 1:numel(dirs)
        if regexp(dirs(j).name,'^C[1-4]_CP_')
          i_caa_files=[i_caa_files j];
          disp([num2str(numel(i_caa_files)) '. ' dirs(j).name]);
        end
      end
      var_file_to_load=irf_ask('file to load >','var_file_to_load',[]);
      if var_file_to_load
        try
          dir_name=dirs(i_caa_files(var_file_to_load)).name;
          disp(['loading ' dir_name]);
          cd(dir_name)
          eval([dir_name '=dataobj(''*.cdf'');'])
          disp('.ok.');
        catch
          disp(['error loading ' dirs(j).name]);
        end
      end
      cd(old_pwd)
    case 'l'
      xx=whos;
      for j=1:size(xx,1)
        if strcmp(xx(j).class,'dataobj')
          if ~strcmp(xx(j).name,'current_caa_file')
            disp(xx(j).name)
          end
        end
      end
    case 'c'
      xx=whos;i_dataobj=[];name_dataobj={};
      for j=1:size(xx,1)
        if strcmp(xx(j).class,'dataobj')
          if ~strcmp(xx(j).name,'current_caa_file')
            i_dataobj=[i_dataobj j];
            disp([num2str(numel(i_dataobj)) '. ' [xx(j).name]])
          end
        end
      end
      var_caa_file_num=irf_ask('[%]>','var_caa_file_num',1);
      name_caa_file=xx(i_dataobj(var_caa_file_num)).name;
      eval(['current_caa_file=' name_caa_file ';']);
      flag_caa_file_menu=1;
      while flag_caa_file_menu
        if flag_variable_menu
          disp('************************************************************************************')
          disp(['********** IRF CAA. CAA_OBJECT: ' name_caa_file])
          disp('m) menu')
          disp('q) quit')
          disp('d) data variables')
          disp('s) support data variables')
          disp('v) list all variables')
          disp('r) read variable')
          disp('f) file variable attributes')
          disp('g) global attributes')
        end
        flag_variable_menu=0;
        var_caa_file_menu=irf_ask(['irf_caa:' name_caa_file '>'],'var_caa_file_menu','');
        switch var_caa_file_menu
          case 'm'
            flag_variable_menu=1;
          case 'q'
            flag_caa_file_menu=0;
          case 'v'
            display(current_caa_file,'full')
          case 'r'
            param=getfield(get(current_caa_file,'VariableAttributes'),'PARAMETER_TYPE');
            i_data=[];
            for j=1:size(param,1)
              yy=param(j,1);i_data=[i_data j]; data_name=yy{1};
              disp([num2str(numel(i_data)) '. ' data_name]);
            end
            var_data_num=irf_ask('Choose data>','var_data_num',1);
            yy=param(i_data(var_data_num),1);data_name=yy{1};
            disp(['Reading: ' data_name]);
            eval([data_name '= get(current_caa_file,''' data_name ''');']);
            nelem_data=eval(['numel(' data_name '.data);']);
            if nelem_data < 50
              disp('Value: ');
              eval(['xx=' data_name '.data;']);
              disp(xx);clear xx;
            else
              eval(['disp(' data_name ');']);
            end
          case 'd'
            param=getfield(get(current_caa_file,'VariableAttributes'),'PARAMETER_TYPE');
            i_data=[];
            for j=1:size(param,1)
              if strcmp(param(j,2),'Data')
                yy=param(j,1);i_data=[i_data j]; data_name=yy{1};
                disp([num2str(numel(i_data)) '. ' data_name]);
              end
            end
            var_data_num=irf_ask('Choose data>','var_data_num',1);
            yy=param(i_data(var_data_num),1);data_name=yy{1};
            flag_caa_data_menu=1;
            while flag_caa_data_menu
              if flag_show_data_menu
                disp('************************************************************************************')
                disp(['********** IRF CAA. CAA_OBJECT: ' name_caa_file '.' data_name])
                disp('m) menu')
                disp('q) quit')
                disp('p) plot data')
                disp('i) read into irfu format')
              end
              flag_show_data_menu=0;
              var_caa_data_menu=irf_ask(['irf_caa:' name_caa_file '.' data_name '>'],'var_caa_data_menu','');
              switch var_caa_data_menu
                case 'm'
                  flag_show_data_menu=1;
                case 'q'
                  flag_caa_data_menu=0;
                case 'i'
                  eval([data_name '= getmat(current_caa_file,''' data_name ''');']);
                case 'p'
                  plot(current_caa_file, data_name);
                otherwise
                  try eval(var_caa_data_menu); catch, end
              end
            end

          case 's'
            param=getfield(get(current_caa_file,'VariableAttributes'),'PARAMETER_TYPE');
            i_data=[];
            for j=1:size(param,1)
              if strcmp(param(j,2),'Support_Data')
                yy=param(j,1);i_data=[i_data j]; data_name=yy{1};
                disp([num2str(numel(i_data)) '. ' data_name]);
              end
            end
          case 'f'
            get(current_caa_file,'VariableAttributes')
            ga=get(current_caa_file,'VariableAttributes');
            ssf=fieldnames(ga);
            for jj=1:length(ssf)
              disp('*************************************')
              disp(ssf{jj});
              disp('*************************************')
              eval(['disp(ga.' ssf{jj} ')'])
            end
          case 'g'
            get(current_caa_file,'GlobalAttributes')
            ga=get(current_caa_file,'GlobalAttributes');
            ssf=fieldnames(ga);
            for jj=1:length(ssf)
              disp('*************************************')
              disp(ssf{jj});
              disp('*************************************')
              eval(['disp(ga.' ssf{jj} ')'])
            end
          otherwise
            try eval(var_caa_file_menu); catch ,end
        end
      end
    otherwise
      eval(var_menu);
  end
end

clear current_caa_file
return

irf_log('load',['cdf file: ' cdf_file]);
if 1
  if strcmp(var_name,'*')
    flag_while=1;
    while flag_while
      disp('======== Options ========')
      disp('q) quit, dont return anything')
      disp('v) list all variables')
      disp('vv) list all variables and their values')
      disp('r) read variable in default format')
      disp('t) read time variables in irfu format and return')
      var_additional_item=irf_ask('Variable? [%]>','var_additional_item','v');
      switch var_additional_item
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
          eval([var_to_read '= spdfcdfread(''' cdf_file ''', ''VARIABLES'', ''' var_to_read ''');']);
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
variables=[epoch_variable var_name];

DATA = spdfcdfread(cdf_file, 'VARIABLES', variables);

temp=struct([DATA{:,1}]);
t = [temp.date]; t = t(:);
t = (t-62167219200000)/1000; %#ok<NASGU>

for k=2:numel(variables)
  clear var;
  for j=1:size(DATA,1)
    temp=DATA{j,k};temp=temp(:)';var(j,:)=temp;
  end
  var=[DATA{:,k}]';
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
  for i=1:nargout, varargout(i) = {eval(variables{i+1})}; end
end
