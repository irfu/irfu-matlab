function [R, V, RGSM, VGSM] = mms_update_ephemeris
% MMS_UPDATE_EPHEMERIS  read DEFEPH files and compute R/V in GSE & GSM
%
% [R, V, RGSM, VGSM] = mms_update_ephemeris
%
% Read DEFEPH files from the start of the mission and up to current time
% and compute spacecraft position (R) and velocity (V) in GSE and GSM
%
% To save the data use:
% save /data/mms/irfu/mmsR.mat R
% save /data/mms/irfu/mmsV.mat V
% save /data/mms/irfu/mmsRGSM.mat RGSM
% save /data/mms/irfu/mmsVGSM.mat VGSM

nowUtc = datestr(now,'yyyy-mm-ddTHH:MM:SSZ'); 
tint = irf.tint(['2015-03-13T00:00:00Z/' nowUtc]);

for mmsId = 1:4
  mmsIdS = num2str(mmsId);
  vPref = ['mms',mmsIdS,'_ancillary_defeph'];
  vars = 'rv';
  for i=1:length(vars)
    v = vars(i);
    Gei = mms.db_get_variable(vPref, v, tint);    
    Gei_ts = TSeries(EpochTT(Gei.time), Gei.data, 'vec_xyz');
    Gei_ts.coordinateSystem = 'gei';
    gse_ts = irf.geocentric_coordinate_transformation(Gei_ts,'gei>gse',false);
    gsm_ts = irf.geocentric_coordinate_transformation(Gei_ts,'gei>gsm',false);
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