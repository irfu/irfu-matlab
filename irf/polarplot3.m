function [h]=polarplot3(t,Ex,Ey,steplength,avnumber,threshold,plotflag)

%POLARPLOT3  Plot polarization parameters for a given signal
%
% [h]=POLARPLOT3(t,Ex,Ey,step,av,thresh,plotflag) plots the polarization
% parameters for a given signal. The polarisation parameters are the
% Horizontal Spectral Intensity, the Horizontal Degree of
% Polarization, the Degree of Pure Circular Polarization and the
% Horizontal Spinor Tilt Angle as defined by Carozzi et al.
% (JGR 106, pp 21395--21408, 2001).
%
% The input parameters are:
% t = time (in epoch)
% Ex, Ey = The x and y components of the signal
% step = size of each FFT
% av = number of spectra to average over
% thresh = degree of polarization above which we consider the
%          signal to be polarized.
% plotflag = set equal to 1 if we want the phase difference
%            and the coherence, instead of the usual output.
%
% See also POLARPLOT2
%

% By Anders Tjulin (Last update 23/10-2002)
%

%% Check the input
%warning off
if nargin>7
  error('Too many arguments')
end
if nargin<7
  plotflag=0;
end
if nargin<6
  threshold=0.7;
end
if nargin<5
  avnumber=4;
end
if nargin<4
  steplength=512;
end
if nargin<3
  error('Too few arguments')
end

Ex=Ex(:);
Ey=Ey(:);

%% Check the sampling rate

sampl=1/(t(2)-t(1));
t0=floor(t(1));

%% Create the Hamming window function

j=1:steplength;
w=0.54-0.46*cos(2*pi*j/(steplength-1));
normalisation=w*w'/steplength;

%% Calculate the spectral matrix with FFT

n=steplength-1;
nd2=ceil(n/2);
for j=1:steplength/2:length(Ex)-steplength
  Yx=fft(detrend(Ex(j:j+steplength-1)).*w.');
  Yy=fft(detrend(Ey(j:j+steplength-1)).*w.');
  %    Yx=fft(Ex(j:j+steplength-1).*w.');
  %    Yy=fft(Ey(j:j+steplength-1).*w.');
  Yx(1)=[];
  Yy(1)=[];
  fEx(:,2*((j-1)/steplength+1)-1)=Yx(1:nd2)/n;
  fEy(:,2*((j-1)/steplength+1)-1)=Yy(1:nd2)/n;
  T(2*((j-1)/steplength+1)-1)=t(j+steplength/2);
end

%% Find the frequencies

nyq=1/2;
freq=sampl*(1:nd2)/(nd2)*nyq;

%% Calculate the polarisation parameters

I=real(fEx.*conj(fEx)+fEy.*conj(fEy))/normalisation;
Q=real(fEx.*conj(fEx)-fEy.*conj(fEy))/normalisation;
U=2*real(fEx.*conj(fEy))/normalisation;
V=2*imag(fEx.*conj(fEy))/normalisation;

%% Create mean values

Itemp=0;Qtemp=0;Utemp=0;Vtemp=0;
for j=1:avnumber
  Itemp=Itemp+[zeros(steplength/2,j-1),I,zeros(steplength/2,avnumber-j)];
  Qtemp=Qtemp+[zeros(steplength/2,j-1),Q,zeros(steplength/2,avnumber-j)];
  Utemp=Utemp+[zeros(steplength/2,j-1),U,zeros(steplength/2,avnumber-j)];
  Vtemp=Vtemp+[zeros(steplength/2,j-1),V,zeros(steplength/2,avnumber-j)];
end
Imean=Itemp(:,1:size(fEx,2))/avnumber;Imean=Imean(:,avnumber:end);
Qmean=Qtemp(:,1:size(fEx,2))/avnumber;Qmean=Qmean(:,avnumber:end);
Umean=Utemp(:,1:size(fEx,2))/avnumber;Umean=Umean(:,avnumber:end);
Vmean=Vtemp(:,1:size(fEx,2))/avnumber;Vmean=Vmean(:,avnumber:end);
T=T(avnumber:end);

root=sqrt(Qmean.*Qmean+Umean.*Umean+Vmean.*Vmean);
p=root./Imean;
vp=Vmean./root;
phi=atan2(Umean,Qmean)*180/pi;
coh=(Umean.^2+Vmean.^2)./(Imean.^2-Qmean.^2);
phase=atan2(Vmean,Umean)*180/pi;

%% Plot everything

if plotflag == 1 %%% The case with phase difference and coherence
  phase(coh<threshold)=NaN;

  h(1)=subplot(4,1,1);
  pcolor(T-t0,freq,log10(2*n*Imean/sampl))
  shading flat
  colorbar('vert');
  ylabel('Hz')
  title('Power spectral density')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out','XTickLabel','');
  xlabel('')
  %set(gca,'yscale','log')

  h(2)=subplot(4,1,2);
  pcolor(T-t0,freq,coh)
  shading flat
  caxis([0,1])
  colorbar('vert')
  ylabel('Hz')
  title('Coherence')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out','XTickLabel','');
  xlabel('')
  %set(gca,'yscale','log')

  h(3)=subplot(4,1,3);
  pcolor(T-t0,freq,phase)
  shading flat
  caxis([-180,180])
  colorbar('vert')
  ylabel('Hz')
  title('Phase difference')
  temp1=axis;
  irf_timeaxis(gca,t0); set(gca,'tickdir','out','XTickLabel','');
  xlabel('')
  %set(gca,'yscale','log')

  h(4)=subplot(4,1,4);
  plot(t-t0,Ex,t-t0,Ey)
  xlabel('UT')
  ylabel('mV/m')
  title('Time series')
  set(gca,'position',[0.13 0.11 0.7023 0.1642])
  grid
  temp2=axis;
  axis([temp1(1:2),temp2(3:4)])
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')

else %%% The standard case
  temp=find(p<threshold);
  vp(temp)=NaN;
  phi(temp)=NaN;

  h(1)=subplot(4,1,1);
  pcolor(T-t0,freq,log10(Imean))
  shading flat
  colorbar('vert');
  ylabel('Hz')
  title('Spectral intensity')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out','XTickLabel','');
  xlabel('')
  %set(gca,'yscale','log')

  h(2)=subplot(4,1,2);
  pcolor(T-t0,freq,p)
  shading flat
  caxis([0,1])
  colorbar('vert')
  ylabel('Hz')
  title('Degree of polarisation')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out','XTickLabel','');
  xlabel('')
  %set(gca,'yscale','log')

  h(3)=subplot(4,1,3);
  pcolor(T-t0,freq,vp)
  shading flat
  caxis([-1,1])
  colorbar('vert')
  ylabel('Hz')
  title('Ellipticity')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out','XTickLabel','');
  xlabel('')
  %set(gca,'yscale','log')

  h(4)=subplot(4,1,4);
  pcolor(T-t0,freq,phi)
  shading flat
  caxis([-180,180])
  colorbar('vert');
  xlabel('UT')
  ylabel('Hz')
  title('Tilt angle')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
end

%suptitle(datestr(epoch2date(t0),31))
warning on
