function shifted_signal = dft_timeshift(sig, tau, Freq)
% shifted_signal = dft_timeshift(sig, tau, [Freq])
%
% Shifts the input signal "sig" by "tau" seconds using discrete fourier
% transform (DFT). Particularly useful when calculating the frequency-wavenumber 
% spectrum of the mms' spin-plane or axial probes. See: mms.fk_powerspectrum.
%
% Input:    sig - TSeries or array to be shifted (Note: Only tensor order 1)
%
%           tau - Applied shift in time. Positive: T->T+tau, Negative:
%                 T->T-tau.
%           Freq - Optional, Frequency of data sample.
%                  If unspecified calculated using TSeries time data.
%                  If input signal is a normal array, Freq must be
%                  specified for time shifts, otherwise tau is assumed to
%                  be a shift in index rather than time.
%
% Output:   shifted_signal - TSeries or array that has been shifted.
%
% Example 1: TSeries shifted back 4 seconds.
% t0=int64(528973090634000000);
% samples=int64(0:127);
% time=samples*10^9;
% data=sin(pi*2*(0:127)/16);
% shift = -4; %Shift the TSeries back 4 seconds.
% Input_TSeries = TSeries(EpochTT(time),data');
% Output_TSeries = mms.dft_timeshift(Input_TSeries,shift);
% h=irf_plot({Input_TSeries, Output_TSeries});
% irf_zoom(h,'x',[Output_TSeries.time(1),Input_TSeries.time(end)]);

% Input check and variable assignment:
narginchk(2,3);
if isa(sig,'TSeries')
  tseries_flag=true;
  inTime = sig.time;
  sig = sig.data;
  if nargin<3
    Freq = round(1/(median(inTime(2:end)-inTime(1:end-1))));
  end
else
  tseries_flag=false;
end

if exist('Freq','var')
  dT = 1/Freq;   % Frequency and time between samples in Hz and s.
  dS = tau/dT;   % Applied delay in samples.
else
  % If input is not a TSeries, and the sample frequency is unspecified,
  % the shift "tau" is assumed to be in samples, and not time.
  dS = tau;
end

% Performing the DFT:
U = fft(sig);
N_p = numel(sig);

if mod(N_p,2)==0, U(N_p/2+1)=0; end % Disregard Nyquist frequency for even-sized dft

f = (mod( ( (0:N_p-1) + floor(N_p/2) ), N_p) - floor(N_p/2))/N_p;
outData = ifft(U.*exp(-2i*pi*dS.*f)');

% Output:
if tseries_flag
  shifted_signal = TSeries(inTime+tau, outData);
else
  shifted_signal = outData;
end

end
