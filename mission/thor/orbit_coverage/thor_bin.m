function [tsDist,edgesData] = thor_bin(tsData,edgesData,edgesTime)

nTints = edgesTime.length-1;
tCenter = edgesTime(1:nTints) + 0.5*(edgesTime(2:(nTints+1))-edgesTime(1:nTints));
%tmpTint.start + 0.5*[tmpTint.stop-tmpTint.start];
vecDist = zeros(nTints,numel(edgesData)-1);


for iTint = 1:nTints
  tmpTint = edgesTime(iTint+[0 1]); % time between two perigee
  
  tmpData = tsData.tlim(tmpTint);  
  
  % sample into bins
  nDist = histc(tmpData.data,edgesData);
  vecDist(iTint,:) = nDist(1:end-1);        
end

tsDist = irf.ts_scalar(tCenter,vecDist);