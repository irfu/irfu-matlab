function res = mms_saa(defatt)
%MMS_SAA  Compute Solar Aspect Angle from DEFATT
%
%  SAA = MMS_SAA(DEFATT)
%
%  returns SAA(TSeries) from DEFATT as read by mms_load_ancillary()
%
%  See also: MMS_LOAD_ANCILLARY

tt2000 = EpochTT2000(defatt.time);
[x,y,z] = sph2cart(defatt.zra*pi/180,defatt.zdec*pi/180,1);
out = irf.geocentric_coordinate_transformation(...
  [tt2000.toEpochUnix().epoch x y z],'gei>gse');

saa = atand(out(:,4)./out(:,2));
res = TSeries(tt2000,saa);