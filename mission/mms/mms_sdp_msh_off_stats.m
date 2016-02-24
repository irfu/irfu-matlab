% Script to compute MSH DSL offsets for all existing ROI 
%
% See also mms_sdp_comp_msh_dsl_off

%% Initialize
roiList = mms_sci_roi();
clear shDslOffDb
shDslOffDb(size(roiList,1)) = struct('c1',[],'c2',[],'c3',[],'c4',[],'tint',[]);
iLast = [];

%% Run the process
if isempty(iLast), iLast = size(roiList,1); end
for i=iLast:-1:1
  iLast = i;
  Tint = irf.tint(roiList(i,:));
  disp(['Processing ' Tint.start.utc])
  shDslOffDb(i) = mms_sdp_comp_msh_dsl_off(Tint);
end

%% Create TS from the result
nRec = length(shDslOffDb); clear tt
for i=1:nRec, tt(i)=shDslOffDb(i).tint.start.ttns(); end %#ok<SAGROW>

for mmsId=1:4
  mmsIdS=sprintf('c%d',mmsId);
  off = NaN(nRec,2);
  for i=1:length(shDslOffDb)
    if isempty(shDslOffDb(i).(mmsIdS)), continue, end
    off(i,:) = shDslOffDb(i).(mmsIdS);
  end
  
  OffTS = irf.ts_vec_xy(EpochTT(tt),off);
  OffTS.name = sprintf('mms%d_sh_off',mmsId); OffTS.units = 'mV/m';
  c_eval('Off?=OffTS;', mmsId)
end

%% Remove bad points
c_eval('idx = abs(Off?.data(:,2))>0.5; Off?.data(idx,:)=NaN;')

%% Plot
irf_pl_tx('Off?','x')

%% Write CAL files
for mmsId=1:4
  c_eval('OffTS=Off?;',mmsId)
  fid = fopen(['mms' num2str(mmsId) '_edp_sdp_ql_v' irf_fname(date2epoch(now),3) '00.txt'],'w');
  fprintf(fid,'%s\n','% UTC time [yyyy-mm-ddThh:mm:ss.mmmuuunnnZ]     DSL-X offset    DSL-Y offset, TAB separated');
  for i=nRec:-1:1
    if isnan(OffTS.data(i,1)), continue, end
    fprintf(fid,'%s\t%.2f\t%.2f\n',OffTS.time(i).utc(),OffTS.data(i,1),OffTS.data(i,2));
  end
  fclose(fid);
end