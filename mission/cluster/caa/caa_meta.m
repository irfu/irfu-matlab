function out=caa_meta(varargin)
% CAA_META return meta data structure
%
%	CAA_META(string1,string2,..) returns name of datasets
%			matching string1,string2,..
%           If index file with all metadata is not accessible, it is
%           downloaded from the irfu server.
%
%	CAA_META(dataset) displays dataset informations
%	out=CAA_META(dataset) returns dataset information in structure out
%
%	CAA_META(variableName) displays CAA variableName informations
%	out=CAA_META(variableName) returns variable information in structure out
%
%	CAA_META('create') create new index file downloading all the meta data
%	from CAA. This may take time! If you want to generate new index file
%	for the irfu-matlab server, please see the comments inside the file.
%
%
%	Examples:
%		caa_meta('C1','PEA')
%		caa_meta C1 PEA
%		d=caa_meta('C4_CP_PEA_MOMENTS')
%		caa_meta('B_Vec_xyz_ISR2__C1_CP_EFW_L2_BB')
%

%% Defaults and variables
persistent s datasetNames indexFile
% If you are generating new index file for the irfu server
%   then please increase the version number
indexFileDefault = 'indexCaaMeta_v3'; % default file name, v2 added in 20130618
linkUrlFile = ['https://www.space.irfu.se/cluster/matlab/' ...
  indexFileDefault '.mat'];
%% empty arguments > show help
if nargin==0
  help caa_meta;
  return;
end

%% Create index file
if nargin>=1 && ischar(varargin{1}) && strcmp(varargin{1},'create')
  irf.log('warning','Getting all metadata from CAA, be very patient...');
  urlMetaData = 'https://csa.esac.esa.int/csa-sl-tap/data';
  csaID = get_url_identity;
  % Create temporary directory, download and unpack files
  tempDir = tempname;
  mkdir(tempDir);
  cd(tempDir);
  tempFileName  = tempname(tempDir);
  tempFileTarGz = [tempFileName '.tar.gz'];
  tempFileTar   = [tempFileName '.tar'];
  if verLessThan('matlab', '8.4')
    tempGetRequest = { 'Username', csaID.csaUser, ...
      'password', csaID.csaPwd, ...
      'RETRIEVAL_TYPE', 'HEADER', ...
      'DATASET_ID', '*'};
    [tempFileTarGz,isOk] = urlwrite(urlMetaData, tempFileTarGz, ...
      'Authentication', 'Basic', 'Get', tempGetRequest ); %#ok<URLWR> websave introduced in R2014b
  else
    tempGetRequest = {'RETRIEVAL_TYPE', 'HEADER', ...
      'DATASET_ID', '*'};
    webOpt = weboptions('Username', csaID.csaUser, ...
      'Password', csaID.csaPwd, ...
      'RequestMethod', 'get', ...
      'Timeout', Inf);
    try
      tempFileTarGz = websave(tempFileTarGz, urlMetaData, ...
        tempGetRequest{:}, webOpt);
      isOk = true;
    catch
      isOk = false;
    end
  end
  if isOk
    gunzip(tempFileTarGz);
    %untar(tempFileTar,'./');
    if( any( strcmp(version, {'8.6.0.267246 (R2015b)','9.0.0.307022 (R2016a) Prerelease'}) ) )
      if(isunix)
        % Matlab 8.6 R2015b and 9.0 R2016a Prerelease have a bug in Matlabs untar, try to use system built in tar on unix machines
        irf.log('notice','Matlab R2015b Release or R2016a Prerelease have bug in Matlab untar, trying built in system tar on Mac/Linux.');
        cmd = sprintf('tar -xf %s --directory ./', tempFileTar);
        status = system(cmd);
        if(status~=0)
          irf.log('critical','Failed untar using built in system command. Trying Matlab, but this may fail with long file names.');
          untar(tempFileTar,'./');
        end
      else
        irf.log('warning','Matlab R2015b Release or R2016a Prerelease on Windows, trying built in Matlab untar. This may fail with long file names!');
        untar(tempFileTar,'./');
      end
    else
      untar(tempFileTar,'./');
    end
  else
    irf.log('critical','Did not succeed downloading file');
    return;
  end

  % Create structures with all the metadata
  irf.log('warning','Creating structures');
  d = dir;
  ii=arrayfun(@(x) any(strfind(x.name,'CSA')),d);
  xmlDir = d(ii).name;
  d = dir([xmlDir '/*.XML']);
  isub = [d(:).isdir];
  d(isub)=[];
  nameFiles = {d.name}';
  nameDatasetList = cell(numel(nameFiles),1); % needed to construct argument for saving at the end
  datasetList = cell(numel(nameFiles),1);
  for iFile=1:numel(nameFiles)
    fileName=nameFiles{iFile};
    if fileName(1)==' ', disp(['Starts with space:' fileName]);continue;end % in case space in beginning of filename
    nameDataset=fileName(1:find(fileName=='.',1)-1);
    nameDataset(strfind(nameDataset,'-')) = '_';
    datasetList{iFile} = nameDataset;
    nameDatasetList{iFile} = [nameDataset ' ']; % need to add space at end because it is used concatenating for save
    disp(nameDataset);
    s=xml2struct([xmlDir '/' fileName]);
    if isfield(s,'DATASETS')
      s=s.DATASETS;
    end
    if isfield(s,'DATASET_METADATA')
      s=s.DATASET_METADATA;
    end
    meta.(nameDataset) = s;  %#ok<STRNU>
  end

  % Remove temporary directory, save metadata structure, display
  % information
  cd ..;
  rmdir(tempDir,'s')
  eval(['save -v7 ' indexFileDefault ' -struct meta ', horzcat(nameDatasetList{:})])
  eval(['save     ' indexFileDefault ' datasetList -append'])
  indexFileName = indexFileDefault; %#ok<NASGU> % save also indexFileName for version control
  save(indexFileDefault,'indexFileName','-append');
  disp(['Please move ''' indexFileDefault '.mat'' file from the current directory to']);
  disp('somewhere on your MATLAB path.');
  disp('If you are updating the version accessible by everybody, then')
  disp(['move file to: ' linkUrlFile]);
  return
end

%% Locating index file
if isempty(indexFile) || ~exist(indexFile,'file')
  if exist(indexFile,'file') % first usage
    indexFileName = load(indexFile,'indexFileName');
    if isempty(indexFileName) ...% indexFileName is not save (old version)
        || ~strcmpi(indexFileName,indexFileDefault) % file is not the newest version
      disp('WARNING!!! There is newer index file for caa_meta available.');
      reply = input('Shall I download? Y/N [Y]:','s');
      if isempty(reply)
        reply = 'Y';
      end
      if strcmpi(reply,'Y')
        indexFile = irf.get_file(linkUrlFile,'caa',indexFileDefault);
      else
        irf.log('warning','Using old verison of index file');
      end
    end
  else % file does not exist and has to be downloaded
    indexFile = [indexFileDefault '.mat']; % default file name
    if ~exist(indexFile,'file')
      indexFile = irf.get_file(linkUrlFile,'caa',indexFileDefault);
    end
  end
end

%% Read metadata
if 	nargin==1 && ischar(varargin{1}) && any(strfind(varargin{1},'__')) % CAA variable
  varName=strrep(varargin{1},'-','_');
  dd=regexp(varName, '__', 'split');
  %	dd{2}=strrep(dd{2},'-','_');
  if numel(varName)>2 && strcmp(varName(1:2),'x3'), varName(1)=[];end
  dd{2} = upper(dd{2});
  metaData = getfield(load(indexFile,dd{2}),dd{2});
  par=metaData.PARAMETERS.PARAMETER;
  iVar= cellfun(@(x) strcmp(strrep(x.PARAMETER_ID.Text,'-','_'),varName),par);
  parVar=par{iVar};
  display_fields(parVar);
  if nargout==1, out=parVar;end
else % dataset
  if isempty(datasetNames)
    datasetNames = getfield(load(indexFile,'datasetList'),'datasetList');
  end
  iSelected = true(numel(datasetNames),1);
  iEqual    = false(numel(datasetNames),1);
  for jInp=1:numel(varargin)
    filter=varargin{jInp};
    filter(strfind(filter,'-')) = '_';
    if ischar(filter)
      iFind=cellfun(@(x) any(strfind(lower(x),lower(filter))),datasetNames);
      iSelected=iSelected & iFind(:);
      if nargin==1 % check if there is identifcal fit
        iEqual=cellfun(@(x) strcmpi(x,filter),datasetNames);
      end
    end
  end
  if sum(iSelected)>0
    if nargout == 1
      if sum(iSelected)==1
        out = load(indexFile,datasetNames{iSelected});
      elseif any(iEqual)
        out = load(indexFile,datasetNames{iEqual});
      else
        disp('Output not assigned. There are several choices:');
        disp(vertcat(datasetNames(iSelected)));
        out=[];
      end
      return;
    end
    if sum(iSelected)==1 ||  any(iEqual)
      if any(iEqual)
        ind=find(iEqual);
      else
        ind=find(iSelected);
      end
      dataSet=getfield(load(indexFile,datasetNames{ind}),datasetNames{ind});
      try
        display_fields(dataSet);
        parameters=dataSet.PARAMETERS.PARAMETER;
        disp('PARAMETERS:');
        for j=1:numel(parameters)
          disp([num2str(j) '. ' mat_output(parameters{j}.PARAMETER_ID.Text)]);
        end
      catch
      end
      disp(' ');
    end
    if sum(iSelected) > 1
      disp('----');
      disp([num2str(sum(iSelected)) ' datasets correspond selection']);
      cellfun(@(x) fprintf('%s\n',x),vertcat(mat_output(datasetNames(iSelected),1)), 'UniformOutput',false);
      return;
    end

  end
end
end
%% Functions
function display_fields(dataSet)
disp(' ');
fn=fieldnames(dataSet);
for j=1:numel(fn)
  if isfield(dataSet.(fn{j}),'Text')
    if strfind(dataSet.(fn{j}).Text,' ') % text includes many words
      disp([fn{j} ': ' dataSet.(fn{j}).Text]);
    else
      disp([fn{j} ': ' mat_output(dataSet.(fn{j}).Text)]);
    end
  end
end
end
function outStr=mat_output(inStr,forceFlag)
% if string is caa variable output link
% if forceFlag defined and equal to one, force matlab linked output
if nargin == 1
  forceFlag = false;
end
if ~forceFlag && any(strfind(inStr,'__'))
  forceFlag = true;
end

if forceFlag
  if ischar(inStr)
    outStr=['<a href="matlab: caa_meta ' ...
      inStr '">' inStr '</a>' ];
  elseif iscell(inStr)
    outStr = cellfun(@(x) ['<a href="matlab: caa_meta ' ...
      x '">' x '</a>' ],inStr, 'UniformOutput',false);
  end
else
  outStr=inStr;
end
end

function csaID = get_url_identity()
% Users should use their own credentials as far as possible.
csaUser = datastore('csa','user');
if isempty(csaUser)
  csaUser = input('Input csa username:','s');
  if isempty(csaUser)
    disp('Please register at ESA: https://www.cosmos.esa.int/web/csa and then use your own credentials in irfu-matlab.');
  else
    datastore('csa','user',csaUser);
  end
end
csaPwd = datastore('csa','pwd');
if isempty(csaPwd)
  csaPwd = input('Input csa password:','s');
  if isempty(csaPwd) && ~isempty(csaUser)
    disp('Please register at ESA: https://www.cosmos.esa.int/web/csa and then use your own credentials in irfu-matlab.');
  else
    datastore('csa','pwd',csaPwd);
  end
end
if strcmp(csaUser, 'avaivads') && strcmp(csaPwd,'!kjUY88lm')
  % Old password used by irfu-matlab, now (2018/06/18) deprecated!
  % Every user must from now on use their own credentials with ESA.
  datastore('csa','user',[]); datastore('csa','pwd',[]);
  errStr = ['Please register at ESA: https://www.cosmos.esa.int/web/csa', ...
    ' and then use your own credentials in irfu-matlab to download data from CSA.'];
  irf.log('critical', errStr);
  error(errStr);
end
csaID = struct('csaUser', csaUser, 'csaPwd', csaPwd);
end
