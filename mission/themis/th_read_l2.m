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
% SPDX-License-Identifier: Beerware
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
  case 'sp'
    instrumentId = 'state';
    levelId = 'l1';
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
    if ~isempty(varTmp)
      if isempty(res)
        res =  varTmp(varTmp(:,1)>tint(1) & varTmp(:,1)<tint(end),:);
      else
        % join time intervals
        res = [res; ...
          varTmp(varTmp(:,1)>res(end,1) & varTmp(:,1)<tint(end),:)]; %#ok<AGROW>
      end
    end
  end
  epochFileStart = epochFileStart + 3600*24;
  if epochFileStart>=epochFileEnd, break, end
  timeVecStart = fromepoch(epochFileStart);
end
% Remove repeating points.
% Example, THE bs&es 2007-07-07T04:11:30.000000Z -- 2007-07-07T08:39:30.000000Z
if ~isempty(res), res(diff(res(:,1))==0,:) = []; end

  function res = read_var
    res = spdfcdfread(fileToRead,'CombineRecords',true,'Variable',varName);
    if isempty(res), return, end
    depTimeVar = find_depend_time();
    if ~isempty(depTimeVar)
      tmpTime = spdfcdfread(fileToRead,'CombineRecords',true,'Variable',depTimeVar);
      res = [tmpTime double(res)];
    end
    res = res;
    
    function res = find_depend_time
      res = '';
      info = spdfcdfinfo(fileToRead);
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
    files = dir(sprintf('%s%sth%s_%s_%s_%d%02d%02d_v*.cdf',...
      fullPath,filesep,...
      thId,levelId,instrumentId,timeVecStart(1),timeVecStart(2),timeVecStart(3)));
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
