function out=caa_load(varargin)
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
% ok = CAA_LOAD(..) returns true if data loaded or lister and returns false
% if nothing loaded or listed
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

%% Defaults
ok                  = false; % default routine did not succeed
%% Defaults that can be changed by input parameters
shouldOnlyListFiles = false; % default list and load files
shouldReadAllData   = true;  % default load everything
shouldLoadFromFile  = true;  % if false, do not load object from file if it exists in memory
useExactNameMatch   = false; % default is to filter not according to exact match
shouldFilterNames = any(nargin);   % if there is input, use it to filter names, otherwise load all

%% Check input parameters
if nargin > 0 % filter which datasets to load
  i=1;
  datasetFilter=cell(nargin,1);
  for j=1:length(varargin)
    if ischar(varargin{j}) && ~isempty(strfind(varargin{j},'tint='))
      shouldReadAllData=false;
      tint = eval(varargin{j});
      datasetFilter(j:end)=[]; % remove last cells from variable file
      break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'tint') ...
        && length(varargin)>j && isnumeric(varargin{j+1})
      shouldReadAllData=false;
      tint=varargin{j+1};
      datasetFilter(j:end)=[]; % remove last cells from variable file
      break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'list') % only list whats available
      shouldReadAllData   = false;
      shouldOnlyListFiles = true;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'nowildcard') % load only specified names
      useExactNameMatch = true;
    elseif ischar(varargin{j}) && strcmpi(varargin{j},'ifnotinmemory') % load only specified names
      shouldLoadFromFile = false;
    elseif ischar(varargin{j})
      if strfind(varargin{j},'__') % variable name specified as input
        dd=regexp(varargin{j}, '__', 'split');
        datasetFilter{i}=dd{end};
      else
        datasetFilter{i}=varargin{j};
      end
      datasetFilter{i}(strfind(datasetFilter{i},'-'))='_'; % substitute '-' to '_'
      i=i+1;
    end
  end
  if i==1, shouldFilterNames=0;end % no names found
  datasetFilter(i:end)=[];
end

%% Check data directory
if isdir([pwd filesep 'CAA']) % check if is CAA folder, then assume data are there
  dirs = dir('CAA');
  caaDataDirectory='CAA/';
else % otherwise assume one is in the CAA data folder
  dirs = dir;
  caaDataDirectory='';
end

%% Load the data
nloaded = 0;
for j = 1:numel(dirs)
  if regexp(dirs(j).name,'^C[1-4,L]_(C|J|P|S)(P|Q|T)_')
    datasetName = dirs(j).name;
    datasetName(strfind(datasetName,'-'))='_'; % substitute '-' to '_'
    shouldLoadVariable = true;
    if shouldFilterNames % if there is name filtering required check if to load variable
      for jj=1:length(datasetFilter)
        if useExactNameMatch
          if strcmpi(datasetName,datasetFilter{jj})
            shouldLoadVariable = true;
            break; % exact match found, read the variable
          else
            shouldLoadVariable = false;
          end
        elseif isempty(strfind(datasetName,datasetFilter{jj}))
          shouldLoadVariable = false;
        end
      end
      if shouldLoadVariable && shouldOnlyListFiles
        disp(datasetName);
        shouldLoadVariable = false;
        ok=true;
      end
    else % work on all variables
      if shouldOnlyListFiles
        shouldLoadVariable = false;
        disp(datasetName);
        if ~isempty(datasetName), ok = true; end
      else
        shouldLoadVariable = true;
      end
    end

    if shouldLoadVariable
      try
        if shouldReadAllData && ~shouldLoadFromFile && evalin('caller',['exist(''' datasetName ''',''var'')'])
          irf.log('warning',[datasetName ' exist in memory. NOT LOADING FROM FILE!'])
        else
          irf.log('warning',['loading ' datasetName ' from location:' caaDataDirectory dirs(j).name filesep '*.cdf']);
          if shouldReadAllData
            evalin('caller',[datasetName '=dataobj(''' caaDataDirectory dirs(j).name filesep '*.cdf'');']);
          else
            assignin('caller','caa_load_tint_temp',tint);
            evalin('caller',[datasetName '=dataobj(''' caaDataDirectory dirs(j).name filesep '*.cdf'',''tint'',caa_load_tint_temp);']);
            evalin('caller','clear caa_load_tint_temp');
          end
        end
        nloaded = nloaded + 1;
      catch
        irf.log('critical',['Did not succeed! Error loading: ' datasetName]);
      end
    end
  end
end
if nloaded > 0, ok = true; end
if nloaded == 0 && ~shouldOnlyListFiles
  irf.log('warning','CAA_LOAD : nothing to load')
end
if nargout > 0, out = ok; end
