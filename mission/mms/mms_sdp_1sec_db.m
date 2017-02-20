function [out,Tint] = mms_sdp_1sec_db(fileName,dataPath)
% MMS_SDP_1SEC_DB  produce 1 sec data bases on one l2a edp file
%
% OUT = MMS_SDP_1SEC_DB(FILENAME,DATAPATH)
%
% Load L2a file and compute data at 1 sec resolution at fixed spin phase.
% Also load FGM, and FPI data.
%
% See also: IRF_FIXED_PHASE_EPOCH

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<2, dataPath = '/data/mms'; end

mmsIdS = fileName(1:4); mmsId = str2double(fileName(4));
l2a = dataobj([dataPath filesep fileName]);
Tint = EpochTT(l2a.data.mms1_edp_epoch_fast_l2a.data([1 end]));
% xxx - check if loaded
ScPot_sitl = irf.ts2mat(mms.get_data('V_edp_fast_sitl',Tint,mmsId));
if isempty(ScPot_sitl)
  ScPot_sitl = irf.ts2mat(mms.get_data('V_edp_fast_l2',Tint,mmsId));
end

%%
Phase = getmat(l2a,[mmsIdS '_edp_phase_fast_l2a']);
PhaseFixed = irf_fixed_phase_epoch(Phase,18,9);
tFixedPha = PhaseFixed(:,1);
dtAv = median(diff(PhaseFixed(:,1)));

dce = getmat(l2a,[mmsIdS '_edp_dce_fast_l2a']);
% ASPOC
MMS_CONST = mms_constants;
aspocBit = double(bitand(l2a.data.mms1_edp_bitmask_fast_l2a.data(:,1),...
  MMS_CONST.Bitmask.ASPOC_RUNNING));
% ScPot
[idx1,idx2] = irf_find_comm_idx(dce,ScPot_sitl);
scPot = nan(size(dce,1),2); scPot(:,1) = dce(:,1);
scPot(idx1,2) = ScPot_sitl(idx2,2);
% Resample - SLOW
dceFixedPha = irf_resamp([dce,aspocBit,scPot(:,2)],tFixedPha,...
  'window',dtAv,'median');

%% B
bfgm = irf.ts2mat(mms.get_data('B_dmpa_srvy_l2',Tint,mmsId));
bfgmFixedPha = irf_resamp(bfgm,tFixedPha,'window',dtAv,'median');

%% FPI Ions
vifpi = irf.ts2mat(mms.get_data('Vi_dbcs_fpi_fast_l2',Tint,mmsId));
nifpi = irf.ts2mat(mms.get_data('Ni_fpi_fast_l2',Tint,mmsId));
TiFPI = mms.get_data('Ti_dbcs_fpi_fast_l2',Tint,mmsId);
tifpi = [TiFPI.time.epochUnix, ...
  double((TiFPI.data(:,1,1)+TiFPI.data(:,2,2)+TiFPI.data(:,3,3)))/3];

ifpiFixedPha = irf_resamp([nifpi,tifpi(:,2),vifpi(:,2:4)],tFixedPha,...
  'window',dtAv,'linear');

%% FPI Electrons
vefpi = irf.ts2mat(mms.get_data('Vi_dbcs_fpi_fast_l2',Tint,mmsId));
nefpi = irf.ts2mat(mms.get_data('Ni_fpi_fast_l2',Tint,mmsId));
TeFPI = mms.get_data('Ti_dbcs_fpi_fast_l2',Tint,mmsId);
tefpi = [TeFPI.time.epochUnix, ...
  double((TeFPI.data(:,1,1)+TeFPI.data(:,2,2)+TeFPI.data(:,3,3)))/3];

efpiFixedPha = irf_resamp([nefpi,tefpi(:,2),vefpi(:,2:4)],tFixedPha,...
  'window',dtAv,'linear');
%%

out = struct('t',tFixedPha,'pha',PhaseFixed(:,2),...
  'e12',dceFixedPha(:,2),'e34',dceFixedPha(:,3),'e56',dceFixedPha(:,4),...
  'aspoc',logical(dceFixedPha(:,5)),'scpot',dceFixedPha(:,6),...
  'b',bfgmFixedPha(:,2:4),...
  'ni',ifpiFixedPha(:,2),'ti',ifpiFixedPha(:,3),'vi',ifpiFixedPha(:,4:6),...
  'ne',efpiFixedPha(:,2),'te',efpiFixedPha(:,3),'ve',efpiFixedPha(:,4:6));
  
end
