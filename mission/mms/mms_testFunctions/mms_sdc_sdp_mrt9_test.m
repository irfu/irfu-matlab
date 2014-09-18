%% Init
data_root='/data/mms/MRT9';
%data_root='/Users/yuri/Dropbox/Projects/MMS/DataProcessing/Data/MRT9';
cd(data_root)
if ~exist('log','dir'), mkdir('log'), end
if ~exist('out','dir'), mkdir('out'), end
%setenv('LOG_PATH_ROOT',[data_root filesep 'log'])
setenv('LOG_PATH_ROOT','')
setenv('DROPBOX_ROOT',[data_root filesep 'out'])
setenv('DATA_PATH_ROOT',[data_root filesep 'out'])

modes={'fast','slow','brst'};
procs={'usc','ql','sitl'};

%% Tests
for scId=2:-1:1
  hk_101File = sprintf('%s/mms%d/mms%d_fields_hk_l1b_101_20150410_v0.1.0.cdf',...
    data_root,scId,scId);
  for iMode=1:length(modes)
    dceFile = sprintf('%s/mms%d/mms%d_sdp_%s_l1b_dce_20150410_v1.0.%d.cdf',...
      data_root,scId,scId,modes{iMode},scId-1);
    dcvFile = sprintf('%s/mms%d/mms%d_sdp_%s_l1b_dcv_20150410_v1.0.%d.cdf',...
      data_root,scId,scId,modes{iMode},scId-1);
    for iProc=1:length(procs)
      fprintf('TEST: MMS%d %s %s\n',scId,modes{iMode},procs{iProc})
      mms_sdc_sdp_proc(procs{iProc},dceFile,dcvFile,hk_101File)
    end
  end
end