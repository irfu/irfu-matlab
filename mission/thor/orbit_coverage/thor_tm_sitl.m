function [tsData] = thor_tm_sitl___(data,downlink,memorySaved_,tsDownlinkDelay,tsAvTM)
% THOR_TM - calculates how much data accumulates onboard and what is downlinked THOR
%   [tsData] = thor_netdata(data,downlink,memoryLimit,downlink_delay)
%   data           - TSeries of generated data where data can be a vector
%                    of different quality data produced onboard. If data is
%                    vector the first indexes are downloaded first.
%   downlink         - how much data are downloaded per orbit
%   memorySaved      - memory limit which is saved from one time step to
%   tsData.saved      - data saved onboard after downlink
%   tsData.discarded  - data discarded	after downlink
%   tsData.downloaded - data downlinked

zeroData       = zeros(size(data.data));
savedData      = zeroData;
downloadedData = zeroData;
discardedData  = zeroData;

nDimData   = size(data.data,2);
nOrbits    = size(data.data,1);

newSavedData = zeroData(1,:);

if isnumeric(downlink)
  downlinkData=ones(nOrbits,1)*downlink;
else
  downlinkData = downlink.data;
end

for iOrbit = 1:(nOrbits-1)
  % Data that has no assigned FOM yet...
  downlinkDelay = tsDownlinkDelay(iOrbit).data;
  memorySaved = memorySaved_ - tsAvTM(iOrbit).data;

  indDelay = 0;
  if iOrbit > (downlinkDelay-1)
    indDelay = downlinkDelay-1;
  end

  % Download highest FOM data from the memorySaved
  if iOrbit > 1
    lastSavedData = savedData(iOrbit,:);
    newSavedData  = lastSavedData;
    while sum(lastSavedData) > 0
      ind=find(lastSavedData>0,1);
      if lastSavedData(ind) < downlinkData(iOrbit) - sum(downloadedData(iOrbit,:))
        downloadedData(iOrbit,ind) = downloadedData(iOrbit,ind) + lastSavedData(ind);
        newSavedData(ind) = 0;
        lastSavedData(ind) = 0;
      else
        tmpDownload = downlinkData(iOrbit) - sum(downloadedData(iOrbit,:));
        downloadedData(iOrbit,ind) = downloadedData(iOrbit,ind) + tmpDownload;
        newSavedData(ind) = newSavedData(ind) - tmpDownload;
        break
      end
    end
  end

  % Add new data to memory and check what is discarded
  discData = data.data(iOrbit-indDelay,:)+newSavedData; % default all
  newSavedData = discData*0;
  for iDim = 1:nDimData % loop over quality factors
    if sum(newSavedData) + discData(iDim)  < memorySaved
      newSavedData(iDim)  = discData(iDim);
      discData(iDim) = 0;
    else
      tmpSave = memorySaved - sum(newSavedData);
      newSavedData(iDim)  = tmpSave;
      discData(iDim) = discData(iDim) - tmpSave;
    end
  end
  discardedData(iOrbit,:) = discData;
  savedData(iOrbit+1,:) = newSavedData;
end

tsData.saved      = data; tsData.saved.data = savedData;
tsData.downloaded = data; tsData.downloaded.data = downloadedData;
tsData.discarded  = data; tsData.discarded.data = discardedData;

