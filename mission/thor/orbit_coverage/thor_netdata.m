function tsNetData = thor_netdata(data,downlink)
% THOR_NETDATA - calculates how much data accumulates onboard THOR
%   netdata = thor_netdata(data,downlink)

% there can be different subcategories, the first category is downlinked 
% first, which should for example highest quality factor

n_dim_data = size(data.data,2); 
nOrbits = size(data.data,1); 

data_stacked = irf.ts_scalar(data.time,cumsum(data.data(:,1:end),2)); % relevant for bowshock
data_stacked_T = irf.ts_scalar(data_stacked.time,cumsum(data_stacked.data,1)); % data acummululated in time

downlink_stacked = irf.ts_scalar(downlink.time,cumsum(downlink.data,2));
downlink_stacked_T = irf.ts_scalar(data.time,cumsum(data.data,1));

netData = NaN(data_stacked.length,n_dim_data);
netDataT = NaN(data_stacked.length,n_dim_data);

netData(1,:) = data_stacked.data(1,:);

for iOrbit = 2:nOrbits  
  new_data = data.data(iOrbit,:);
  if sum(new_data)>0
    1;
  end
  % add the new data
  tot_data = netData(iOrbit-1,:) + new_data;
  tot_data_stacked = cumsum(tot_data,2);
  
  % remove the downlinked data
  to_downlink = downlink_stacked.data(iOrbit,:);
  new_net_data_stacked = tot_data_stacked - to_downlink;  
  % there can be no negative data, so remove all negative values
  new_net_data_stacked(new_net_data_stacked<0) = 0;
  
  new_data = diff([0 new_net_data_stacked]);
  
  
  netData(iOrbit,:) = new_data;

%   tmpdata = data_stacked.data(iOrbit,:);
%   tmpdownlink = downlink_stacked.data(iOrbit,:);
%   tmpdiffQ = tmpdata-tmpdownlink;
%   tmpdiffQ(tmpdiffQ<0) = 0; % there's no negative data
%   netData(iOrbit,:) = tmpdiffQ;
%   
%   tmpdata = dataBSstackedT.data(iOrbit,:);
%   tmpdownlink = downlinkBSstacked.data(iOrbit,:);
%   tmpdiffQ = tmpdata-tmpdownlink;
%   tmpdiffQ(tmpdiffQ<0) = 0; % there's no negative data
%   netDataT(iOrbit,:) = tmpdiffQ;
%   if sum(tmpdata)>0
%     %iOrbit
%     %tmpdata
%    % pause;
%   end
end

tsNetData = irf.ts_scalar(data.time,netData);
%netDataBS = irf.ts_scalar(dataBSstacked.time,netData);
%netDataBSstackedT = irf.ts_scalar(netDataBS.time,cumsum(netDataBS.data,1));
%netDataBST = irf.ts_scalar(netDataBS.time,netDataT);
