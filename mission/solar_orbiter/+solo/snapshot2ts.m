function res = snapshot2ts(d, varName, varargin)
%snapshot2ts  transform a snapshot variable to TimeSeries object
%
%  res = solo.snapshot2ts(DataObj, varName, [FLAGS])
%
% Input:
%   DataObj - dataobj (e.g. LFR_E)
%   varName - name of a variable
%
% Options:
%   'NoNaN' - do not add NaN at the end of each snapshot
%   'RmMedian' - subtract median value from each of the snapshots

flagAddNan = 1; flagMedian = 0;

if nargin > 2, have_options = 1; args = varargin;
else, have_options = 0;
end
while have_options
  l = 1;
  switch lower(args{1})
    case 'nonan'
      flagAddNan = 0;
    case 'rmmedian'
      flagMedian = 1;
    otherwise
      irf.log('warning', ['Option ''' args{1} '''not recognized']);
  end
  if length(args) > l, args = args(l+1:end);
  else, break
  end
end


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
  if nComp>5 % XXX hack
    nComp = size(data,2); lSnapshot = size(data,3);
    data = permute(data,[1 3 2]);
  end
  epoch = d.data.Epoch.data(ii);
  if flagAddNan, lSnapshot = lSnapshot + 1; end
  dtNs = int64((linspace(0, lSnapshot-1, lSnapshot)/double(freq))*1e9);
  tt = (repmat(epoch,1,lSnapshot)+repmat(dtNs,nData,1))';
  tt = tt(:);
  Time = EpochTT(tt);
  dataOut = zeros(Time.length,nComp);
  for i = 1:nComp
    cc = squeeze(data(:,:,i))';
    if flagMedian
      mm = median(cc,1);
      cc = cc - repmat(mm,size(cc,1),1);
    end
    if flagAddNan,  cc(lSnapshot,:) = NaN; end
    dataOut(:,i) = cc(:);
  end
  res.(sprintf('f%d',round(freq))) = irf.ts_scalar(Time,dataOut);
end
