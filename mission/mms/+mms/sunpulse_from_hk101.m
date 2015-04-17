function sps = sunpulse_from_hk101(scId, tint)
%MMS.SUNPULSE_FROM_HK101  Get sunpulse from hk101 files
%
%  SPS = sunpulse_from_hk(scId, tint)
%
%   SPS struct :
%      .time      - Timestamp of HK packet in TT2000
%      .sunpulse  - Latest sunpulse in TT2000
%      .sunssps   - Sunpulse indicator
%      .iifsunper - CIDP sun period (in us)
%
% See also: mms_sdp_phase_2

filePref = sprintf('mms%d_fields_hk_l1b_101',scId);
varPref = sprintf('mms%d_101_',scId);
sunpulse = mms.db_get_variable(filePref,[varPref 'sunpulse'],tint);

sps.time = sunpulse.DEPEND_0.data;
sps.sunpulse = sunpulse.data;

sunssps = mms.db_get_variable(filePref,[varPref 'sunssps'],tint);
if ~(all(sps.time == sunssps.DEPEND_0.data))
  errStr = 'sunpulse and sunssps have different number of records';
  irf.log('critical', errStr), error(errStr)
end
sps.sunssps = sunssps.data;

iifsunper = mms.db_get_variable(filePref,[varPref 'iifsunper'],tint);
if ~(all(sps.time == iifsunper.DEPEND_0.data))
  errStr = 'sunpulse and iifsunper have different number of records';
  irf.log('critical', errStr), error(errStr)
end
sps.iifsunper = iifsunper.data;

