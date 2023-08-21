function bpfields = fftbandpass(varargin)
%
% bpfields = mms.fftbandpass(fielddata,fmin, fmax, fs)
%
% Perform simple bandpass using FFT - returns fields between with fmin < f < fmax
%
% INPUTS:
%   fielddata - Data to be bandpassed fitlered
%   fmin - min. frequency of filter, f < fmin are removed
%   fmax - max. frequency of filter, f > fmin are removed
%   fs -   sampling frequency (optional). If not submitted fs is calculated
%          from the time series.
% OUTPUT:
%   bpfields - Bandpassed filtered data in the same format as input
% CAVIATES - Can be some spurius effects near boundary. Can take longer interval
% then use .tlim to remove
% Written by D. B. Graham

if (length(varargin) < 3)
  help mms.fftbandpass;
  bpfields = NaN;
  return;
end

fielddata = varargin{1};
fmin = varargin{2};
fmax = varargin{3};

if isa(fielddata,'TSeries')
  tmpfields = fielddata.data;
  tmptime = fielddata.time;
else
  tmptime = fielddata(:,1);
  tmpfields = fielddata(:, 2:end);
end

% Make sure number of elements is an even number, if odd remove last
% element to make an even number
Nels = length(tmpfields(:,1));
if (mod(Nels,2) == 1)
  tmpfields = tmpfields(1:end-1, :);
  tmptime = tmptime(1:end-1);
  Nels = length(tmpfields(:,1));
end

numfields = length(tmpfields(1,:));

% Set NaN values to zero so FFT works
inan=isnan(tmpfields);
tmpfields(inan)=0;

%Bandpass filter field data
if (length(varargin) == 4)
  dt = 1/varargin{4};
else
  dt = tmptime(2)-tmptime(1);
end
fN = 1/(2*dt);
df = 2*fN/(Nels);
f = (-fN:df:fN-df);

% FFT and remove frequencies
for nn=1:numfields
  fieldtemp = fft(tmpfields(:,nn));
  fieldtemp = fftshift(fieldtemp);

  fieldtemp(find(abs(f) < fmin)) = 0;
  fieldtemp(find(abs(f) > fmax)) = 0;

  fieldtemp = ifftshift(fieldtemp);
  tmpfields(:,nn) = ifft(fieldtemp);
end

% Put back original NaNs
tmpfields(inan)= NaN;

% Return data in the same format as input
if isa(fielddata,'TSeries')
  bpfields = TSeries(tmptime,tmpfields,'to',fielddata.tensorOrder);
  if fielddata.tensorOrder == 1
    bpfields.coordinateSystem = fielddata.coordinateSystem;
  end
  bpfields.name = fielddata.name;
  bpfields.units = fielddata.units;
  bpfields.userData = fielddata.userData;
else
  bpfields = [tmptime, tmpfields];
end

end