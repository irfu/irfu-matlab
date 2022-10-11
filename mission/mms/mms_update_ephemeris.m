function [R, V, RGSM, VGSM] = mms_update_ephemeris
% MMS_UPDATE_EPHEMERIS  read DEFEPH files and compute R/V in GSE & GSM
%
% [R, V, RGSM, VGSM] = mms_update_ephemeris
%
% Read DEFEPH files from the start of the mission and up to current time
% and compute spacecraft position (R) and velocity (V) in GSE and GSM
%
% NOTE 1: This function read MMS ancillary files (in GEI J2000 coordinate
% system) then transform these to GSE and GSM but does not yet handle
% transformation properly (to GSE/GSM TOD).
% Note 2: For now resample is disabled as irf_resamp uses extrapolation
% (over potentially several orbits whenever DEFEPH is missing).
%
% To save the data use:
% save /data/mms/irfu/mmsR.mat R
% save /data/mms/irfu/mmsV.mat V
% save /data/mms/irfu/mmsRGSM.mat RGSM
% save /data/mms/irfu/mmsVGSM.mat VGSM

nowUtc = char(datetime("now","Format","uuuu-MM-dd'T'HH:mm:ss'Z'","TimeZone","UTCLeapSeconds"));
tint = irf.tint(['2015-03-13T00:00:00Z/' nowUtc]);
%newTime = EpochTT((tint.start.epoch:int64(60*1e9):tint.stop.epoch)');

for mmsId = 1:4
  mmsIdS = num2str(mmsId);
  vPref = ['mms',mmsIdS,'_ancillary_defeph'];
  vars = 'rv';
  for i=1:length(vars)
    v = vars(i);
    Gei = mms.db_get_variable(vPref, v, tint);
    Gei_ts = TSeries(EpochTT(Gei.time), Gei.(v), 'vec_xyz');
    Gei_ts.coordinateSystem = 'gei';
    gse_ts = irf.geocentric_coordinate_transformation(Gei_ts,'gei>gse',false);
    %gse_ts = gse_ts.resample(newTime,'spline');
    gsm_ts = irf.geocentric_coordinate_transformation(Gei_ts,'gei>gsm',false);
    %gsm_ts = gsm_ts.resample(newTime,'spline');
    vC = upper(v);
    switch v
      case 'r'
        R.(['gse' vC mmsIdS]) = gse_ts.data;
        RGSM.(['gsm' vC mmsIdS]) = gsm_ts.data;
        if(mmsId==1)
          R.time = gse_ts.time.ttns;
          RGSM.time = gsm_ts.time.ttns;
        end
      case 'v'
        V.(['gse' vC mmsIdS]) = gse_ts.data;
        VGSM.(['gsm' vC mmsIdS]) = gsm_ts.data;
        if(mmsId==1)
          V.time = gse_ts.time.ttns;
          VGSM.time = gsm_ts.time.ttns;
        end
    end
  end
end
end