function [shifted_signal]=dft_timeshift(varargin)
% [shifted_signal]=dft_timeshift(sig, tau, Freq)
%
% Shifts the input signal "sig" by "tau" seconds using discrete fourier
% transform (DFT). Particularly useful when calculating the frequency-wavenumber 
% spectrum of the mms' spin-plane or axial probes. See: mms.fk_powerspectrum.
%
% Input:    sig - TSeries or array to be shifted (Note: Only tensor order 1)
%
%           tau - Applied shift in time. Positive: T->T+tau, Negative:
%                 T->T-tau.
%           Freq - Frequency of data sample. If unspecified calculated
%                  using TSeries time data.
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
% h=irf_plot(2,'newfigure');
% irf_plot(h(1),Input_TSeries);
% irf_plot(h(2),Output_TSeries);
% irf_zoom(h,'x',[Output_TSeries.time(1),Input_TSeries.time(end)]);

% Input check and variable assignment:
if isempty(varargin)
    help mms.dft_timeshift;
    return
end
if numel(varargin) < 2
    disp('Too few input arguments, see >help mms.dft_timeshift');
    return;
elseif numel(varargin) > 3
    disp('Too many input arguments, see >help mms.dft_timeshift');
    return;
end
if isa(varargin{1},'TSeries')
    tseries_flag=1;
else
    tseries_flag=0;
end

sig=varargin{1};
tau=varargin{2};

if numel(varargin)==3
    Freq=varargin{3};
else
    if tseries_flag
        Freq=round(1/(median(sig.time(2:end)-sig.time(1:end-1))));
    end
end

if tseries_flag
    dT=1/Freq;   % Frequency and time between samples in Hz and s.
    dS = tau/dT; % Applied delay in samples.
    a=sig.data;
else
    if numel(varargin)==3
        dT=1/Freq;
        dS=tau/dT;
        a=sig;
    else
        dS=tau; % If input is not a TSeries, and the sample frequency is unspecified,
        a=sig;  % the shift "tau" is assumed to be in samples, and not time.
    end
end

% Performing the DFT:
U=fft(a);
N_p=numel(a);
if (mod(N_p, 2) == 0)
    U(length(U)/2+1)=0; % Disregard Nyquist frequency for even-sized dft
end
f=(mod(((0:N_p-1)+floor(N_p/2)), N_p)-floor(N_p/2))/N_p;
A=ifft(U.*exp(-2i*pi*dS.*f)');

% Output:
if tseries_flag
    shifted_signal=TSeries(sig.time+tau,A);
else
    shifted_signal=A;
end

end