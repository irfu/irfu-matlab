function res = snapshot2ts(d, varName, varargin)
%snapshot2ts  transform a snapshot variable to TimeSeries object
%
%  res = jui.snapshot2ts(DataObj, varName, [FLAGS])
%
% Input:
%   DataObj - dataobj (e.g. DATA)
%   varName - name of a variable
%
% Options:
%   'NoNaN' - do not add NaN at the end of each snapshot
%   'RmMedian' - subtract median value from each of the snapshots

flagAddNan = 0; flagMedian = 0;

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

  snapNr = d.data.SNAPSHOT_NUMBER.data(ii);
  dd = diff(double(snapNr)); idxSnap = find(dd~=0);
  idx0 = find(ii>0); idx0 = idx0(1);
  idxSnap = [idx0; idxSnap+1]; nSnapshot = length(idxSnap);

  data = d.data.(varName).data(ii,:,:);
  data(data<-1e30) = NaN;
  nData = size(data,1);
  lBlock = size(data,2);
  nComp = size(data,3);
  if nComp>8 % XXX hack
    nComp = size(data,2); lBlock = size(data,3);
    data = permute(data,[1 3 2]);
  end
  epoch = d.data.Epoch.data(ii);

  lSnapshot = lBlock*nData/(nSnapshot);
  if flagAddNan, lSnapshot = lSnapshot + 1; end
  dtNs = int64((linspace(0, lSnapshot-1, lSnapshot)/double(freq))*1e9);
  tt = (repmat(epoch(idxSnap),1,lSnapshot)+repmat(dtNs,nSnapshot,1))';
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