function instfrequency = irf_instantaneousfrequency(varargin)
% IRF_INSTANTANEOUSFREQUENCY Compute the instantaneous frequency of a signal using
% change in phase from the Hilbert transformed signal
%
% Written by D. B. Graham
%
% Input:
%       inputsignal - TSeries of signal (e.g., E or B fields)
%
% Options:
%       numptsav -  Number of points averaged.
%       avmethod -  Method used in moving average: "movmean" (default) | "movmedian" | "gaussian" | "lowess" |
%                   "loess" | "rlowess" | "rloess" | "sgolay"
%       frange -    select frequency range to bandpass filter before
%                   calculating frequency
%
% Output:
%       instfrequency - Instantaneous ordinary frequency in TSeries format with the
%       same resolution as the input signal
%
% Example:
%   instfrequency = irf_instantaneousfrequency(Ehmfe,'numptsav',1000,'avmethod','movmedian','frange',[15 20]*1e3);

averageintime = false;
filterdata = false;

inputsignal = varargin{1};
inputsignal.data = double(inputsignal.data);

avmethod = 'movmean';

args=varargin(2:end);

if numel(args)>0
  haveoptions=1;
else
  haveoptions=0;
end

while haveoptions
  l = 2;
  switch(lower(args{1}))
    case 'numptsav'
      if numel(args)>1 && isnumeric(args{2})
        numptsaverage = args{2};
        averageintime = true;
      end
    case 'frange'
      if numel(args)>1 && isnumeric(args{2})
        frange = args{2};
        filterdata = true;
      end
    case 'avmethod'
      if ischar(args{2})
        avmethod = args{2};
      end
    otherwise
      irf.log('warning',['Unknown flag: ' args{1}]);
      l=1;
      break;
  end
  args = args(l+1:end);
  if isempty(args), haveoptions=0; end
end

dt = median(diff(inputsignal.time.epochUnix));

if filterdata
  inputsignal = inputsignal.filt(frange(1),frange(2),1/dt,5);
end

numfields = length(inputsignal.data(1,:));

freqdata = zeros([length(inputsignal.data(:,1))-1 numfields]);


for ii=1:numfields
  hilbertsignal = hilbert(inputsignal.data(:,ii));
  phase = unwrap(angle(hilbertsignal));
  freqdata(:,ii) = diff(phase)/dt/2/pi;
end

freqtimes = inputsignal.time(1:end-1)+dt/2;

instfrequency = irf.ts_scalar(freqtimes,freqdata);
instfrequency = instfrequency.resample(inputsignal);

if averageintime
  for ii=1:numfields
    instfrequency.data(:,ii) = smoothdata(instfrequency.data(:,ii),avmethod,numptsaverage);
  end
end

end