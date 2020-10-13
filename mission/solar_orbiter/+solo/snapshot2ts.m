function res = snapshot2ts(d, varName)
%snapshot2ts  transform a snapshot variable to TimeSeries object
%
%  res = solo.snapshot2ts(DataObj, varName)
%
% Input: 
%   DataObj - dataobj (e.g. LFR_E)
%   varName - name of a variable
%

res = struct;
freqs = unique(d.data.SAMPLING_RATE.data);

for iFreq = 1:length(freqs)
  freq = freqs(iFreq);
  ii = d.data.SAMPLING_RATE.data == freq;
  data = d.data.(varName).data(ii,:,:);
  data(data<-1e30) = NaN;
  nData = size(data,1);
  lSnapshot = size(data,2);
  nComp = size(data,3);
  epoch = d.data.Epoch.data(ii);
  dtNs = int64((linspace(0, lSnapshot-1, lSnapshot)/double(freq))*1e9);
  tt = (repmat(epoch,1,lSnapshot)+repmat(dtNs,nData,1))';
  tt = tt(:);
  Time = EpochTT(tt);
  dataOut = zeros(Time.length,nComp);
  for i = 1:nComp
    cc = squeeze(data(:,:,i))';
    dataOut(:,i) = cc(:);
  end
  res.(sprintf('f%d',freq)) = irf.ts_scalar(Time,dataOut);
end
