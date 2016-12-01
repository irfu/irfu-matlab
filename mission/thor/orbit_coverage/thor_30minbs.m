function [ts,ind] = thor_30minbs(tsQ,T)

dataQ = tsQ.data;
indCrossings = find(dataQ>0);
nCrossings = numel(indCrossings);

crossQ = tsQ(indCrossings);

t_30min = tsQ.time.start:60*T:tsQ.time.stop;
crossQ_30min = crossQ.resample(t_30min,'max');
is30minAv = find(~isnan(crossQ_30min.data));
crossQ_30minNotNaN = crossQ_30min(is30minAv);

newDataQ = zeros(tsQ.length,1);
for iCross = 1:crossQ_30minNotNaN.length 
  tintBS = crossQ_30minNotNaN(iCross).time + T*[-0.5 0.5]*60;
  indBS = tsQ.time.tlim(tintBS);
  newDataQ(indBS,1) = repmat(crossQ_30minNotNaN.data(iCross,1),numel(indBS),1);
end  
ind = find(newDataQ>0);

% Make TSeries
ts = irf.ts_scalar(tsQ.time,newDataQ);



% tsBSCpar = QBpar; allInd = 1:tsBSCpar.length;
% noCrossing = setdiff(allInd,bsCrossing);
% tsBSCpar.data(noCrossing,1) = 0;
% tsBSCper = QBperp; allInd = 1:tsBSCper.length;
% noCrossing = setdiff(allInd,bsCrossing);
% tsBSCper.data(noCrossing,1) = 0;
%
% nBS = numel(bsCrossing);
% Qpar15 = zeros(rTHOR.length,1)-1;
% Qper15 = zeros(rTHOR.length,1)-1;
% for iCross = 1:nBS  
%   Qpar15(bsCrossing(iCross)+[-15:15]) = repmat(Qpar(bsCrossing(iCross)).data,31,1);
%   Qper15(bsCrossing(iCross)+[-15:15]) = repmat(Qper(bsCrossing(iCross)).data,31,1);
% end  
% isWithin15min = find(Qpar15>0);
% iKSR_Within15min = iKSR;
% iKSR_Within15min.data(isWithin15min,1) = 6;
% 