function out = polarization(t,Ex,Ey,sizeFFT,avnumber)

%POLARIZATION  Calculate polarization parameters for a given equidistant signal
%
% out = POLARIZATION(t,Ex,Ey,step,av,thresh,plotflag) calculates the polarization
% parameters for a given signal. The polarisation parameters are the
% Horizontal Spectral Intensity, the Horizontal Degree of
% Polarization, the Degree of Pure Circular Polarization and the
% Horizontal Spinor Tilt Angle as defined by Carozzi et al.
% (JGR 106, pp 21395--21408, 2001).
%
% INPUT parameters are:
%  t        = time (in epoch)
%  Ex, Ey   = The x and y components of the signal
%  sizeFFT  = size of each FFT (default 256)
%  av       = number of spectra to average over (default 4)
%
% OUTPUT is structure with following fields:
%  out.tiltAngleDeg
%  out.degreeOfPolarisation
%  out.ellipticity
%  out.coherence
%  out.phaseDeg
%  out.powerSpectralDensity
%  out.frequency
%  out.time
%
% See also POLARPLOT3
%

% Based on Anders Tjulin code polarplot3
%

%% Check the input
%warning off
if nargin>7
  error('Too many arguments')
end
if nargin<5 || isempty(avnumber)
  avnumber=4;
end
if nargin<4 || isempty(sizeFFT)
  sizeFFT=512;
end
if nargin<3
  error('Too few arguments')
end

Ex=Ex(:);
Ey=Ey(:);

%% Check the sampling rate

sampl = 1/(t(2)-t(1));

%% Create the Hamming window function

j=1:sizeFFT;
w=0.54-0.46*cos(2*pi*j/(sizeFFT-1));
normalisation=w*w'/sizeFFT;

%% Calculate the spectral matrix with FFT

n=sizeFFT-1;
nd2=ceil(n/2);
for j=1:sizeFFT/2:length(Ex)-sizeFFT
  Yx=fft(detrend(Ex(j:j+sizeFFT-1)).*w.');
  Yy=fft(detrend(Ey(j:j+sizeFFT-1)).*w.');
  Yx(1)=[];
  Yy(1)=[];
  fEx(:,2*((j-1)/sizeFFT+1)-1)=Yx(1:nd2)/n;
  fEy(:,2*((j-1)/sizeFFT+1)-1)=Yy(1:nd2)/n;
  T(2*((j-1)/sizeFFT+1)-1)=t(j+sizeFFT/2);
end

%% Find the frequencies

nyq=1/2;
freq = sampl*(1:nd2)/(nd2)*nyq;

%% Calculate the polarisation parameters

I=real(fEx.*conj(fEx)+fEy.*conj(fEy))/normalisation;
Q=real(fEx.*conj(fEx)-fEy.*conj(fEy))/normalisation;
U=2*real(fEx.*conj(fEy))/normalisation;
V=2*imag(fEx.*conj(fEy))/normalisation;

%% Create mean values

Itemp=0;Qtemp=0;Utemp=0;Vtemp=0;
for j=1:avnumber
  Itemp=Itemp+[zeros(sizeFFT/2,j-1),I,zeros(sizeFFT/2,avnumber-j)];
  Qtemp=Qtemp+[zeros(sizeFFT/2,j-1),Q,zeros(sizeFFT/2,avnumber-j)];
  Utemp=Utemp+[zeros(sizeFFT/2,j-1),U,zeros(sizeFFT/2,avnumber-j)];
  Vtemp=Vtemp+[zeros(sizeFFT/2,j-1),V,zeros(sizeFFT/2,avnumber-j)];
end
Imean=Itemp(:,1:size(fEx,2))/avnumber; Imean=Imean(:,avnumber:end);
Qmean=Qtemp(:,1:size(fEx,2))/avnumber; Qmean=Qmean(:,avnumber:end);
Umean=Utemp(:,1:size(fEx,2))/avnumber; Umean=Umean(:,avnumber:end);
Vmean=Vtemp(:,1:size(fEx,2))/avnumber; Vmean=Vmean(:,avnumber:end);
T=T(avnumber:end);

root=sqrt(Qmean.*Qmean+Umean.*Umean+Vmean.*Vmean);

out.tiltAngleDeg         = atan2d(Umean,Qmean)';
out.degreeOfPolarisation = (root./Imean)';
out.ellipticity          = (Vmean./root)';
out.coherence            = ((Umean.^2+Vmean.^2)./(Imean.^2-Qmean.^2))';
out.phaseDeg             = atan2d(Vmean,Umean)';
out.powerSpectralDensity = (2*n*Imean/sampl)';
out.frequency            = freq;
out.time                 = T(:);

%suptitle(datestr(epoch2date(t0),31))
warning on
