function caa_load(varargin)
%CAA_LOAD Script to load data downloaded from the CAA in CDF format.
%Downloaded zip file must be unpacked, and script must be run
% from a CAA_Download_YYYYMMDD_hhmm directory  
% or directory where is subdirectories with all the CAA data object directories 
% or directory where is subdirectory CAA with all the CAA data object directories. 
%
% CAA_LOAD('string1','string2',..)
%   load only data objects that match string1 & string2 &...
% CAA_LOAD('string1','string2',...,'tint',tint)
%   load only time interval tint (good for large files, reads always from file)
% CAA_LOAD('string1','string2',...,'list')
%   list only the data objects that are on disk 
% CAA_LOAD('string1','string2',...,'nowildcard')
%   list only the data objects whos name are exactly string1, string2 etc.
%
%  Examples:
%   caa_load
%   caa_load CIS
%   caa_load C1 CIS
%   caa_load Sweep_Energy__C3_CP_PEA_PITCH_SPIN_PSD
%
%   tint=[irf_time([2007 9 2 14 25 0]) irf_time([2007 9 2 18 40 0])];
%   caa_load('C1_CP_FGM','tint',tint);
%   caa_load FGM list;
%
% $Id$

% add 'nowildcard' option to read in just specified data object

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
flag_read_all=1;        % default load everything
flag_only_list_files=0; % default do not only list files
flag_filter = 1;        % default is to filter according to names
flag_exact_match =0;    % default is to filter not according to exact match

if nargin==0, 
    flag_filter=0; % load all variables, no filtering
end

if nargin > 0, % filter which variables to load
  i=1;
  variable_filter=cell(nargin,1);
  for j=1:length(varargin),
    if ischar(varargin{j}) && ~isempty(strfind(varargin{j},'tint=')),
      flag_read_all=0;
      tint=eval(varargin{j});
      variable_filter(j:end)=[]; % remove last cells from variable file
      break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'tint') ...
        && length(varargin)>j && isnumeric(varargin{j+1}),
      flag_read_all=0;
      tint=varargin{j+1};
      variable_filter(j:end)=[]; % remove last cells from variable file
      break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'list'), % only list whats available
      flag_read_all=0;
      flag_only_list_files=1;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'nowildcard'), % load only specified names
      flag_exact_match=1;  
    elseif ischar(varargin{j})
      if strfind(varargin{j},'__') % variable name specified as input
        dd=regexp(varargin{j}, '__', 'split');
        variable_filter{i}=dd{end};
      else
        variable_filter{i}=varargin{j};
      end
      variable_filter{i}(strfind(variable_filter{i},'-'))='_'; % substitute '-' to '_'
      i=i+1;
    end
  end
  if i==1, flag_filter=0;end % no names found
  variable_filter(i:end)=[];
end

if isdir([pwd filesep 'CAA']), % check if is CAA folder, then assume data are there
  dirs = dir('CAA');
  caa_data_directory='CAA/';
else % otherwise assume one is in the CAA data folder
  dirs = dir;
  caa_data_directory='';
end

nloaded = 0;
for j = 1:numel(dirs)
  if regexp(dirs(j).name,'^C[1-4,L]_(C|P|S)(P|Q)_')
    var_name = dirs(j).name;
    var_name(strfind(var_name,'-'))='_'; % substitute '-' to '_'
    flag_load_variable=1;
    if flag_filter==1, % if there is name filtering required check if to load variable
      for jj=1:length(variable_filter),
        if flag_exact_match==1,
            if strcmpi(var_name,variable_filter{jj}),
                flag_load_variable=1;
                break; % exact match found, read the variable
            else
                flag_load_variable=0;
            end
        elseif isempty(strfind(var_name,variable_filter{jj})),
          flag_load_variable=0;
        end
      end
      if flag_load_variable && flag_only_list_files
        disp(var_name);
        flag_load_variable=0;
      end
    else % work on all variables 
        if flag_only_list_files,
            flag_load_variable=0;
            disp(var_name);
        else
            flag_load_variable=1; 
        end
    end

    if flag_load_variable,
      try
        irf_log('dsrc',['caa_load ' var_name]);
        if flag_read_all && evalin('caller',['exist(''' var_name ''',''var'')']),
          irf_log('dsrc','Variable exist in memory. NOT LOADING FROM FILE!')
        else
          if flag_read_all,
            evalin('caller',[var_name '=dataobj(''' caa_data_directory dirs(j).name filesep '*.cdf'');']);
          else
            assignin('caller','caa_load_tint_temp',tint);
            evalin('caller',[var_name '=dataobj(''' caa_data_directory dirs(j).name filesep '*.cdf'',''tint'',caa_load_tint_temp);']);
            evalin('caller','clear caa_load_tint_temp');
          end
        end
        nloaded = nloaded + 1;
      catch
        irf_log('dsrc','Did not succeed!');
        irf_log('dsrc',['error loading ' var_name]);
      end
    end
  end
end
if nloaded, 
  irf_log('dsrc',['=====> loaded ' num2str(nloaded) ' variables']);
elseif ~flag_only_list_files
  irf_log('dsrc','CAA_LOAD : nothing to load')
end