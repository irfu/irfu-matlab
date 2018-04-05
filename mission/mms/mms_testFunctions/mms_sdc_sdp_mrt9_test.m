%% Init
data_root='/data/mms/MRT9';
%data_root='/Users/yuri/Dropbox (IRFU)/Projects/MMS/DataProcessing/Data/MRT9';
%data_root='/home/thoni/MMS/MMS_cdf/MRT9';

if ~exist(data_root,'dir')
  error('DATA_ROOT(%s) does not exist',data_root)
end

if ismac
  setenv('CDF_BASE','/Applications/cdf35_0-dist/')
  outDir = '/Users/yuri/tmp_out';
else
  outDir = data_root;
end

cd(outDir)
if ~exist('log','dir'), mkdir('log'), end
if ~exist('out','dir'), mkdir('out'), end
setenv('DROPBOX_ROOT', [outDir filesep 'out'])
setenv('DATA_PATH_ROOT', [outDir filesep 'out'])
setenv('LOG_PATH_ROOT', [outDir filesep 'log'])

modes    = {'slow','brst'};
%modes   = {'slow','fast'};
versions = {'3.5.0'}; % if multiple: add extra for loop below.
procs    = {'scpot','ql','sitl','l2pre'};
dates    = {'20160101'};
defattdate = {'2015365_2016001'}; % Defatt Date corresponding to each date in "dates".
%dates   = {'20150410','20160101'};
scList   = 1;

%% Tests
for scId = scList
  for iMode=1:length(modes)
    for iDate = 1:numel(dates)
      dceFile = sprintf('%s/mms%d/mms%d_edp_%s_l1b_dce_%s*_v%s.cdf',...
        data_root, scId, scId, modes{iMode}, dates{iDate}, versions{1});
      dcvFile = sprintf('%s/mms%d/mms%d_edp_%s_l1b_dcv_%s*_v%s.cdf',...
        data_root, scId, scId, modes{iMode}, dates{iDate}, versions{1});
      hk_101File = sprintf('%s/mms%d/mms%d_fields_hk_l1b_101_%s*_v*.cdf',...
        data_root, scId, scId, dates{iDate});
      hk_10eFile = sprintf('%s/mms%d/mms%d_fields_hk_l1b_10e_%s*_v*.cdf',...
        data_root, scId, scId, dates{iDate});
      defattFile = sprintf('%s/mms%d/MMS%d_DEFATT_%s.V*', ...
        data_root, scId, scId, defattdate{iDate});
      d = dir(dceFile);
      if isempty(d), dceFile = '';
      else, dceFile = sprintf('%s/mms%d/%s', data_root, scId, d(1).name);
      end
      d = dir(dcvFile);
      if isempty(d), dcvFile = '';
      else, dcvFile = sprintf('%s/mms%d/%s', data_root, scId, d(1).name);
      end
      d = dir(hk_101File); % HK with sunpulse, use highest version found.
      if isempty(d), hk_101File='';
      else, hk_101File = sprintf('%s/mms%d/%s', data_root, scId, d(end).name);
      end
      d = dir(hk_10eFile); % HK with bias settings, use highest version found.
      if isempty(d), hk_10eFile='';
      else, hk_10eFile = sprintf('%s/mms%d/%s', data_root, scId, d(end).name);
      end
      d = dir(defattFile); % Defatt, used in L2pre, use highest version found.
      if isempty(d), defattFile = '';
      else, defattFile = sprintf('%s/mms%d/%s', data_root, scId, d(end).name);
      end
      if isempty(dcvFile) && isempty(dceFile), continue, end
      for iProc=1:length(procs)
        fprintf('TEST: %s MMS%d %s %s\n',...
          dates{iDate}, scId, modes{iMode}, procs{iProc})
        if(strcmp(procs{iProc},'l2pre'))
          % L2/L2pre use defatt
          mms_sdc_sdp_proc(procs{iProc}, dceFile, dcvFile, hk_10eFile, defattFile);
        else
          % L1B use hk_101 sunpulse
          mms_sdc_sdp_proc(procs{iProc}, dceFile, dcvFile, hk_10eFile, hk_101File);
        end
      end
    end
  end
end
