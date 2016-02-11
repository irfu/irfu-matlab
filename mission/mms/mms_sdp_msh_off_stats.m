roiList = mms_sci_roi();
clear shDslOffDb
shDslOffDb(size(roiList,1)) = struct('c1',[],'c2',[],'c3',[],'c4',[],'tint',[]);
iLast = [];
%%
if isempty(iLast), iLast = size(roiList,1); end
for i=iLast:-1:1
  iLast = i;
  Tint = irf.tint(roiList(i,:));
  disp(['Processing ' Tint.start.utc])
  shDslOffDb(i) = mms_sdp_comp_msh_dsl_off(Tint);
end