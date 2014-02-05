function res = th_read_l2(varName,tint)
%TH_READ_L2  read local THEMIS L2 data
%
% res = th_read_l2(varName,tint)
%   Read variable varName from THEMIS CDF files located in local
%   repository. Examples of variable names are:
%   thX_fgs_dsl,thX_efs_dot0_dsl,thX_fgl_dsl, thX_eff_dot0_dsl
%
% Example:
%   bs = th_read_l2('the_fgs_dsl',iso2epoch('2008-04-20T22:20:00Z')+[0 3600]);

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

dataDir = '/data/themis';

thId = lower(varName(3));
switch lower(varName(5:6))
  case 'fg'
    instrumentId = 'fgm';
    levelId = 'l2';
  case 'ef'
    instrumentId = 'efi';
    levelId = 'l2';
  case 'sc'
    instrumentId = 'scm';
    levelId = 'l2';
  otherwise
    error('This instrument is not implemented')
end

DTXTRA = 600; % time to add at each side of the interval
tint = tint + DTXTRA*[-1 1];
timeVecStart = fromepoch(tint(1));
epochFileStart = toepoch([timeVecStart(1:3) 0 0 0]);
timeVecEnd = fromepoch(tint(end)+24*3600);
epochFileEnd = toepoch([timeVecEnd(1:3) 0 0 0]);
res = [];
while true
  fileToRead = find_file_to_read();
  if ~isempty(fileToRead)
    irf.log('notice',['reading ' fileToRead])
    varTmp = read_var();
    if isempty(res)
      res =  varTmp(varTmp(:,1)>tint(1) & varTmp(:,1)<tint(end),:);
    else
      % join time intervals
      res = [res; ...
        varTmp(varTmp(:,1)>res(end,1) & varTmp(:,1)<tint(end),:)]; %#ok<AGROW>
    end
  end
  epochFileStart = epochFileStart + 3600*24;
  if epochFileStart>=epochFileEnd, break, end
  timeVecStart = fromepoch(epochFileStart);
end

  function res = read_var
    tmpData = cdfread(fileToRead,'CombineRecords',true,'Variable',varName);
    depTimeVar = find_depend_time();
    if ~isempty(depTimeVar)
      tmpTime = cdfread(fileToRead,'CombineRecords',true,'Variable',depTimeVar);
      tmpData = [tmpTime double(tmpData)];
    end
    res = tmpData;
    
    function res = find_depend_time
      res = '';
      info = cdfinfo(fileToRead);
      for iVar = 1:length(info.VariableAttributes.DEPEND_TIME)
        if strcmpi(info.VariableAttributes.DEPEND_TIME{iVar,1},varName)
          res = info.VariableAttributes.DEPEND_TIME{iVar,2};
          return;
        end
      end
    end % find_depend_time()
    
  end % read_var
    
  function res = find_file_to_read
    res = '';
    fullPath = sprintf('%s%sth%s%s%s%s%s%s%02d',dataDir,filesep,...
      thId,filesep,levelId,filesep,instrumentId,filesep,timeVecStart(1));
    files = dir(sprintf('%s%sth%s_l2_%s_%d%02d%02d_v*.cdf',...
      fullPath,filesep,...
      thId,instrumentId,timeVecStart(1),timeVecStart(2),timeVecStart(3)));
    if ~isempty(files)
      maxVer = 0; fileIdx = 1;
      for iFile = length(files)
        ver = str2double(files(iFile).name(end-5:end-4));
        if ver>maxVer
          maxVer = ver;
          fileIdx = iFile;
        end
      end
      res = [fullPath,filesep,files(fileIdx).name];
    end
  end % find_file_to_read()
end