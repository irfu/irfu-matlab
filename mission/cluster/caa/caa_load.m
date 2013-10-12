function caa_load(varargin)
%CAA_LOAD Script to load data downloaded from the CAA in CDF format.
%Downloaded zip file must be unpacked, and script must be run
% from a CAA_Download_YYYYMMDD_hhmm directory  
% or directory where is subdirectories with all the CAA data object directories 
% or directory where is subdirectory CAA with all the CAA data object directories. 
%
% CAA_LOAD('string1','string2',..)
%   load from disk data objects that match string1 & string2 &...
% CAA_LOAD('string1','string2',...,'tint',tint)
%   load only time interval tint (good for large files, reads always from file)
% CAA_LOAD('string1','string2',...,'list')
%   list the data objects that are on disk 
% CAA_LOAD('string1','string2',...,'nowildcard')
%   load from disk the data objects whos name are exactly string1, string2 etc.
% CAA_LOAD('string1','string2',...,'ifnotinmemory')
%   load from disk only those data that are not in memory
% CAA_LOAD('string1','string2',...,'cfa')
%   load data primarily from CFA folder if it exists instead of CAA folder
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

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
shouldReadAllData=true;   % default load everything
listFilesOnly=false;      % default list and load files
filterNames = any(nargin);% if there is input, use it to filter names, otherwise load all
useExactNameMatch = false;% default is to filter not according to exact match
forceLoadFromFile = true; % if false, do not load object from file if it exists in memory

if nargin > 0, % filter which variables to load
  i=1;
  variable_filter=cell(nargin,1);
  for j=1:length(varargin),
    if ischar(varargin{j}) && ~isempty(strfind(varargin{j},'tint=')),
      shouldReadAllData=false;
      tint=eval(varargin{j});
      variable_filter(j:end)=[]; % remove last cells from variable file
      break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'tint') ...
        && length(varargin)>j && isnumeric(varargin{j+1}),
      shouldReadAllData=false;
      tint=varargin{j+1};
      variable_filter(j:end)=[]; % remove last cells from variable file
      break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'list'), % only list whats available
      shouldReadAllData=false;
      listFilesOnly=1;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'nowildcard'), % load only specified names
      useExactNameMatch=1;  
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'ifnotinmemory'), % load only specified names
      forceLoadFromFile=0;  
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
  if i==1, filterNames=0;end % no names found
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
  if regexp(dirs(j).name,'^C[1-4,L]_(C|J|P|S)(P|Q|T)_')
    var_name = dirs(j).name;
    var_name(strfind(var_name,'-'))='_'; % substitute '-' to '_'
    flag_load_variable=1;
    if filterNames==1, % if there is name filtering required check if to load variable
      for jj=1:length(variable_filter),
        if useExactNameMatch==1,
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
      if flag_load_variable && listFilesOnly
        disp(var_name);
        flag_load_variable=0;
      end
    else % work on all variables 
        if listFilesOnly,
            flag_load_variable=0;
            disp(var_name);
        else
            flag_load_variable=1; 
        end
    end

	if flag_load_variable,
		try
			if shouldReadAllData && ~forceLoadFromFile && evalin('caller',['exist(''' var_name ''',''var'')']),
				irf.log(2,[var_name ' exist in memory. NOT LOADING FROM FILE!'])
			else
				irf.log(2,['loading ' var_name ' from location:' caa_data_directory dirs(j).name filesep '*.cdf']);
				if shouldReadAllData,
					evalin('caller',[var_name '=dataobj(''' caa_data_directory dirs(j).name filesep '*.cdf'');']);
				else
					assignin('caller','caa_load_tint_temp',tint);
					evalin('caller',[var_name '=dataobj(''' caa_data_directory dirs(j).name filesep '*.cdf'',''tint'',caa_load_tint_temp);']);
					evalin('caller','clear caa_load_tint_temp');
				end
			end
			nloaded = nloaded + 1;
		catch
			irf.log(1,['Did not succeed! Error loading: ' var_name]);
		end
	end
  end
end
if nloaded == 0 && ~listFilesOnly
  irf.log(2,'CAA_LOAD : nothing to load')
end