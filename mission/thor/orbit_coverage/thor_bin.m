function [tsDist] = thor_bin(tsData,edgesData,edgesTime)
% THOR_BIN distribution of data into time and data bins
%
% tsDist = THOR_BIN(tsData,edgesData,edgesTime)
%
% tsDist is TSeries object with time being center times of edgesTime and
% data being vector distribution counts within edgesData

nTints  = edgesTime.length-1;
tCenter = edgesTime(1:nTints) + ...
  0.5*(edgesTime(2:(nTints+1))-edgesTime(1:nTints));

vecDist = zeros(nTints,numel(edgesData)-1);

for iTint = 1:nTints
  tmpTint = edgesTime(iTint+[0 1]); % time between two perigee

  tmpData = tsData.tlim(tmpTint);

  % sample into bins
  if edgesData(end) > edgesData(1) % increasing edges
    nDist = histcounts(tmpData.data,edgesData);
    vecDist(iTint,:) = nDist(1:end);
  else
    nDist = histcounts(tmpData.data,edgesData(end:-1:1));
    vecDist(iTint,:) = nDist(end:-1:1);
  end
end

tsDist = irf.ts_scalar(tCenter,vecDist);