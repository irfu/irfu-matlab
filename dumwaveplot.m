function dumwaveplot(t,Ex)
  
%DUMWAVEPLOT "Fast" simple plot of wavelet spectrograms
%
% DUMWAVEPLOT(t,E) makes a wavelet spectrogram of the time
% series. It uses a Morlet wavelet.
% You will need to make changes in the program to suit your needs,
% but please keep this base program intact.
% 
% t = time in seconds
% E = the time series
%
%
% By Anders Tjulin (Last update 15/4-2003)
%
  
  %% Check the input

  if nargin>2
    error('Too many arguments')
  end
  if nargin<2
    error('Too few arguments')
  end
  
  %% Turn data into a column vector
  
  Ex=Ex(:);

  %% Check the sampling rate
  
  sampl=1/(t(2)-t(1));
  t0=t(1);
  
  %% Remove the last sample if the total number of samples is odd

  if length(Ex)/2 ~= floor(length(Ex)/2)
    Ex=Ex(1:end-1);
    t=t(1:end-1);
  end

  %% Find the frequencies for an FFT of all data

  nd2=length(Ex)/2; 
  nyq=1/2;
  freq=sampl*(1:nd2)/(nd2)*nyq;
  w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

  %%
  %% Set some important parameters
  %%
  
  amin=0.01; % The highest frequency to consider is 0.5*sampl/10^amin
  amax=2; % The lowest frequency to consider is 0.5*sampl/10^amax
  anumber=400; % The number of frequencies
  a=logspace(amin,amax,anumber);
%  a=logspace(0.01,2.4,100);
  w0=sampl/2; % The maximum frequency
  sigma=5.36/w0; % The width of the Morlet wavelet
  
  %% Make the FFT of all data

  Sw=fft(Ex);
%  Sw=fft(detrend(Ex));
  [aa,ww]=meshgrid(a,w); % Matrix form
  [aa,Sww]=meshgrid(a,Sw); % Matrix form

  %% Calculate the FFT of the wavelet transform
  
  Ww=sqrt(1).*Sww.*exp(-sigma*sigma*((aa.*ww-w0).^2)/2);

  %% Get the wavelet transform by IFFT of the FFT
  
  W=ifft(Ww); 

  %% Get the correct frequencies for the wavelet transform
  
  newfreq=w0./a;
  [newfreqmat,temp]=meshgrid(newfreq,w);

  %% Calculate the power spectrum
  
  power=(2*pi)*conj(W).*W./newfreqmat;
 
  %% Remove data possibly influenced by edge effects
 
  censur=floor(2*a);
  power2=power;
  for j=1:anumber;
    power2(1:censur(j),j)=NaN;power2(length(Ex)-censur(j):length(Ex),j)=NaN;
  end

 %% Plot everything

%  pcolor(t-t0,newfreq,log10(abs(power.'))) % With edge effects
  pcolor(t-t0,newfreq,log10(abs(power2.'))) % With edge effects removed
  shading flat
  ylabel('Hz')
%  irf_timeaxis(gca,t0); set(gca,'tickdir','out'); % For time in epoch
  set(gca,'yscale','log')
  set(gca,'tickdir','out')
  
%global svar
%svar=power2;
%global frekvens
%frekvens=newfreq;