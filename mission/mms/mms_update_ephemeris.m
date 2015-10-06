function [R, V, RGSM, VGSM] = mms_update_ephemeris
%MMS_UPDATE_EPHEMERIS  read DEFEPH files and compute R/V in GSE & GSM
%
% [R, V, RGSM, VGSM] = mms_update_ephemeris
%
% Read DEFEPH files from the start of the mission and up to current time
% and compute spacecrfat position (R) and velocity (V) in GSE and GSM
%
% To save the data use:
% save /data/mms/irfu/mmsR.mat R
% save /data/mms/irfu/mmsV.mat V
% save /data/mms/irfu/mmsRGSM.mat RGSM
% save /data/mms/irfu/mmsVGSM.mat VGSM

[~,nowUtc] = unix('TZ=UTC date "+%Y-%m-%dT%H:%M:%SZ"');
nowUtc(end)=[]; % clip the cr in the end

tint = irf.tint(['2015-03-13T00:00:00Z/' nowUtc]);
epoch0 = tint.start.epoch;
newTime = (epoch0:int64(60*1e9):tint.stop.epoch)';
R.time = newTime; V.time = newTime;
RGSM.time = newTime; VGSM.time = newTime;

for mmsId = 1:4
  mmsIdS = num2str(mmsId);
  vPref = sprintf('mms%d_ancillary_defeph',mmsId);
  
  vars = 'rv';
  for i=1:length(vars)
    v = vars(i);
    Gei=mms.db_get_variable(vPref,v,tint);
    GeiR = irf_resamp([double(Gei.time-epoch0)*1e-9 Gei.(v)],...
      double(newTime-epoch0)*1e-9);
    
    tTmp = GeiR(:,1)+EpochTT(epoch0).epochUnix;
    gse = irf.geocentric_coordinate_transformation(...
      [tTmp GeiR(:,2:4)],'gei>gse');
    
    gsm = irf.geocentric_coordinate_transformation(...
      [tTmp GeiR(:,2:4)],'gei>gsm');
    
    vC = upper(v);
   switch v
     case 'r'
       R.(['gse' vC mmsIdS]) = gse(:,2:4); 
       RGSM.(['gsm' vC mmsIdS]) = gsm(:,2:4);
     case 'v'
       V.(['gse' vC mmsIdS]) = gse(:,2:4); 
       VGSM.(['gsm' vC mmsIdS]) = gsm(:,2:4);
   end
  end
end