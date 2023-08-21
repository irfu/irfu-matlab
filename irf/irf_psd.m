function [Pxx, Pxxc, f] = irf_psd(varargin)
%IRF_PSD  Power Spectral Density estimate
%
%   Pxx = IRF_PSD(X,NFFT,Fs,WINDOW) estimates the Power Spectral Density of
%   signal vector X using Welch's averaged periodogram method.  X is
%   divided into overlapping sections, each of which is detrended, then
%   windowed by the WINDOW parameter, then zero-padded to length NFFT.
%   The magnitude squared of the length NFFT DFTs of the sections are
%   averaged to form Pxx.  Pxx is length NFFT/2+1 for NFFT even, (NFFT+1)/2
%   for NFFT odd, or NFFT if the signal X is complex.  If you specify a
%   scalar for WINDOW, a Hanning window of that length is used.  Fs is the
%   sampling frequency which doesn't affect the spectrum estimate but is
%   used for scaling of plots.
%
%   [Pxx,F] = IRF_PSD(X,NFFT,Fs,WINDOW,NOVERLAP) returns a vector of frequen-
%   cies the same size as Pxx at which the PSD is estimated, and overlaps
%   the sections of X by NOVERLAP samples.
%
%   [Pxx, Pxxc, F] = IRF_PSD(X,NFFT,Fs,WINDOW,NOVERLAP,P) where P is a scalar
%   between 0 and 1, returns the P*100% confidence interval for Pxx.
%
%   IRF_PSD(X,...,DFLAG), where DFLAG can be 'linear', 'mean' or 'none',
%   specifies a detrending mode for the prewindowed sections of X.
%   DFLAG can take the place of any parameter in the parameter list
%   (besides X) as long as it is last, e.g. PSD(X,'mean');
%
%   IRF_PSD with no output arguments plots the PSD in the current figure window,
%   with confidence intervals if you provide the P parameter.
%
%   The default values for the parameters are NFFT = 256 (or LENGTH(X),
%   whichever is smaller), NOVERLAP = 0, WINDOW = HANNING(NFFT), Fs = 2,
%   P = .95, and DFLAG = 'none'.  You can obtain a default parameter by
%   leaving it off or inserting an empty matrix [], e.g. PSD(X,[],10000).
%
%   See also CSD, COHERE, TFE, PMEM, PMTM, PMUSIC.
%   ETFE, SPA, and ARX in the Identification Toolbox.
%
%   $Id$

%   Author(s): T. Krauss, 3-26-93
%   Copyright (c) 1988-98 by The MathWorks, Inc.
%   $Revision$  $Date$

%   The units on the power spectra Pxx and Pyy are such that, using
%   Parseval's theorem:
%
%        SUM(Pxx)/LENGTH(Pxx) = SUM(X.^2)/LENGTH(X) = COV(X)
%
%   The RMS value of the signal is the square root of this.
%   If the input signal is in Volts as a function of time, then
%   the units on Pxx are Volts^2*seconds = Volt^2/Hz.
%
%   Here are the covariance, RMS, and spectral amplitude values of
%   some common functions:
%         Function   Cov=SUM(Pxx)/LENGTH(Pxx)   RMS        Pxx
%         a*sin(w*t)        a^2/2            a/sqrt(2)   a^2*LENGTH(Pxx)/4
%Normal:  a*rand(t)         a^2              a           a^2
%Uniform: a*rand(t)         a^2/12           a/sqrt(12)  a^2/12
%
%   For example, a pure sine wave with amplitude A has an RMS value
%   of A/sqrt(2), so A = SQRT(2*SUM(Pxx)/LENGTH(Pxx)).
%
%   See Page 556, A.V. Oppenheim and R.W. Schafer, Digital Signal
%   Processing, Prentice-Hall, 1975.

narginchk(1,7)
xx = varargin{1}; % xx because x will be later within the loop
if min(size(xx))==1, xx=xx(:); end % make sure x is column vector

[msg,nfft,Fs,window,noverlap,p,dflag]=psdchk(varargin(2:end),xx(:,1));
error(msg)

if Fs==2 % test sampling frequency in case Fs not given
  if xx(1,1)> 9e8 % most probably first column is isdat time
    Fs=1./(xx(2,1)-xx(1,1));
    ii_start=2;
  else
    ii_start=1;
  end
else
  % seems the sampling frequency was given
  ii_start=2;
end
for ii=ii_start:size(xx,2) % cycle over the columns of x
  x=xx(:,ii);

  % compute PSD
  window = window(:);
  n = length(x);		% Number of data points
  nwind = length(window); % length of window
  if n < nwind    % zero-pad x if it has length less than the window length
    x(nwind,1:end)=0;  n=nwind;
  end

  k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
  % (k = fix(n/nwind) for noverlap=0)

  index = 1:nwind;
  KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased
  % KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

  Spec = zeros(nfft,1);
  for i=1:k
    if strcmp(dflag,'none')
      xw = window.*(x(index));
    elseif strcmp(dflag,'linear')
      xw = window.*detrend(x(index));
    else
      xw = window.*detrend(x(index),0);
    end
    index = index + (nwind - noverlap);
    Xx = abs(fft(xw,nfft)).^2;
    Spec = Spec + Xx;
  end

  % Select first half
  if ~any(any(imag(x)~=0))   % if x is not complex
    if rem(nfft,2)    % nfft odd
      select = (1:(nfft+1)/2)';
    else
      select = (1:nfft/2+1)';
    end
    Spec = Spec(select);
  else
    select = (1:nfft)';
  end
  freq_vector = (select - 1)*Fs/nfft;

  % find confidence interval if needed
  if (nargout == 3) || ((nargout == 0) && ~isempty(p))
    if isempty(p)
      p = .95;    % default
    end
    % Confidence interval from Kay, p. 76, eqn 4.16:
    % (first column is lower edge of conf int., 2nd col is upper edge)
    confid = Spec*chi2conf(p,k)/KMU;

    if noverlap > 0
      disp('Warning: confidence intervals inaccurate for NOVERLAP > 0.')
    end
  end

  % I hate decibells, added to make correct units - AV 97.12.20
  Spec=Spec/Fs*2;
  %confid=confid/Fs*2;
  %
  % original line
  Spec = Spec*(1/KMU);   % normalize

  % set up output parameters
  if (nargout == 3)
    Pxx = Spec;
    Pxxc = confid;
    f = freq_vector;
  elseif (nargout == 2)
    Pxx = Spec;
    Pxxc = freq_vector;
  elseif (nargout == 1)
    Pxx = Spec;
  elseif (nargout == 0)
    if ~isempty(p)
      P = [Spec confid];
    else
      P = Spec;
    end
    newplot;
    %   plot(freq_vector,10*log10(abs(P))), grid on
    co=get(gca,'colororder');
    clr=co(mod(ii,size(co,1)),1:3);
    loglog(freq_vector,(abs(P)),'color',clr), grid on
    hold on
    %   xlabel('Frequency [Hz]'), ylabel('Power Spectrum Magnitude (dB)');
    xlabel('Frequency [Hz]');
    try
      q=c_desc(inputname(1));
      units=['[' q.units{1} ']^2/Hz'];
    catch
      units='(singal unit)^2/Hz';
    end
    ylabel(units);
    hold off
  end
end


function [msg,nfft,Fs,window,noverlap,p,dflag] = psdchk(P,x,y)
%PSDCHK Helper function for PSD, CSD, COHERE, and TFE.
%   [msg,nfft,Fs,window,noverlap,p,dflag]=PSDCHK(P,x,y) takes the cell
%   array P and uses each element as an input argument.  Assumes P has
%   between 0 and 7 elements which are the arguments to psd, csd, cohere
%   or tfe after the x (psd) or x and y (csd, cohere, tfe) arguments.
%   y is optional; if given, it is checked to match the size of x.
%   x must be a numeric vector.
%   Outputs:
%     msg - error message, [] if no error
%     nfft - fft length
%     Fs - sampling frequency
%     window - window vector
%     noverlap - overlap of sections, in samples
%     p - confidence interval, [] if none desired
%     dflag - detrending flag, 'linear' 'mean' or 'none'

%   Author(s): T. Krauss, 10-28-93
%   Copyright (c) 1988-98 by The MathWorks, Inc.
%       $Revision$  $Date$

msg = [];

if isempty(P)
  % psd(x)
  nfft = min(length(x),256);
  window = hanning(nfft);
  noverlap = 0;
  Fs = 2;
  p = [];
  dflag = 'none';
elseif length(P) == 1
  % psd(x,nfft)
  % psd(x,dflag)
  if isempty(P{1}),   dflag = 'none'; nfft = min(length(x),256);
  elseif ischar(P{1}), dflag = P{1};       nfft = min(length(x),256);
  else,              dflag = 'none'; nfft = P{1};
  end
  Fs = 2;
  window = hanning(nfft);
  noverlap = 0;
  p = [];
elseif length(P) == 2
  % psd(x,nfft,Fs)
  % psd(x,nfft,dflag)
  if isempty(P{1}), nfft = min(length(x),256); else, nfft=P{1};     end
  if isempty(P{2}),   dflag = 'none'; Fs = 2;
  elseif ischar(P{2}), dflag = P{2};       Fs = 2;
  else,              dflag = 'none'; Fs = P{2};
  end
  window = hanning(nfft);
  noverlap = 0;
  p = [];
elseif length(P) == 3
  % psd(x,nfft,Fs,window)
  % psd(x,nfft,Fs,dflag)
  if isempty(P{1}), nfft = min(length(x),256); else, nfft=P{1};     end
  if isempty(P{2}), Fs = 2;     else,    Fs = P{2}; end
  if ischar(P{3})
    dflag = P{3};
    window = hanning(nfft);
  else
    dflag = 'none';
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
  end
  noverlap = 0;
  p = [];
elseif length(P) == 4
  % psd(x,nfft,Fs,window,noverlap)
  % psd(x,nfft,Fs,window,dflag)
  if isempty(P{1}), nfft = min(length(x),256); else, nfft=P{1};     end
  if isempty(P{2}), Fs = 2;     else,    Fs = P{2}; end
  window = P{3};
  if length(window) == 1, window = hanning(window); end
  if isempty(window), window = hanning(nfft); end
  if ischar(P{4})
    dflag = P{4};
    noverlap = 0;
  else
    dflag = 'none';
    if isempty(P{4}), noverlap = 0; else, noverlap = P{4}; end
  end
  p = [];
elseif length(P) == 5
  % psd(x,nfft,Fs,window,noverlap,p)
  % psd(x,nfft,Fs,window,noverlap,dflag)
  if isempty(P{1}), nfft = min(length(x),256); else, nfft=P{1};     end
  if isempty(P{2}), Fs = 2;     else,    Fs = P{2}; end
  window = P{3};
  if length(window) == 1, window = hanning(window); end
  if isempty(window), window = hanning(nfft); end
  if isempty(P{4}), noverlap = 0; else, noverlap = P{4}; end
  if ischar(P{5})
    dflag = P{5};
    p = [];
  else
    dflag = 'none';
    if isempty(P{5}), p = .95;    else,    p = P{5}; end
  end
elseif length(P) == 6
  % psd(x,nfft,Fs,window,noverlap,p,dflag)
  if isempty(P{1}), nfft = min(length(x),256); else, nfft=P{1};     end
  if isempty(P{2}), Fs = 2;     else,    Fs = P{2}; end
  window = P{3};
  if length(window) == 1, window = hanning(window); end
  if isempty(window), window = hanning(nfft); end
  if isempty(P{4}), noverlap = 0; else, noverlap = P{4}; end
  if isempty(P{5}), p = .95;    else,    p = P{5}; end
  if ischar(P{6})
    dflag = P{6};
  else
    msg = 'DFLAG parameter must be a string.'; return
  end
end

% NOW do error checking
if (nfft<length(window))
  msg = 'Requires window''s length to be no greater than the FFT length.';
end
if (noverlap >= length(window))
  msg = 'Requires NOVERLAP to be strictly less than the window length.';
end
if (nfft ~= abs(round(nfft))) || (noverlap ~= abs(round(noverlap)))
  msg = 'Requires positive integer values for NFFT and NOVERLAP.';
end
if ~isempty(p)
  if (numel(p)>1) || (p(1,1)>1) || (p(1,1)<0)
    msg = 'Requires confidence parameter to be a scalar between 0 and 1.';
  end
end
if ~isnumeric(x) || length(size(x))>2
  msg = 'Requires vector (either row or column) input.';
end
if (nargin>2) && ( (min(size(y))~=1) || ~isnumeric(y) || length(size(y))>2 )
  msg = 'Requires vector (either row or column) input.';
end
if (nargin>2) && (length(x)~=length(y))
  msg = 'Requires X and Y be the same length.';
end

dflag = lower(dflag);
if strncmp(dflag,'none',1)
  dflag = 'none';
elseif strncmp(dflag,'linear',1)
  dflag = 'linear';
elseif strncmp(dflag,'mean',1)
  dflag = 'mean';
else
  msg = 'DFLAG must be ''linear'', ''mean'', or ''none''.';
end

