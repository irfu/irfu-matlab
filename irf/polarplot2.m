function h=polarplot2(t,Ex1,Ex2,Ey1,Ey2,steplength,avnumber,threshold,plotflag,distance)

clear polarplot2_manager
%POLARPLOT2 Plot polarization parameters for a given signals
%
% POLARPLOT2(t,Ex1,Ex2,Ey1,Ey2,steplength,avnumber,threshold,plotflag,distance) plots the polarization
% parameters for a given signals. The polarisation parameters are the
% Horizontal Spectral Intensity, the Horizontal Degree of
% Polarization, the Degree of Pure Circular Polarization and the
% Horizontal Spinor Tilt Angle as defined by Carozzi et al.
% (JGR 106, pp 21395--21408, 2001).
%
% CLICK ON THE PHASE SPECTRA TO GET LINE-PLOTS
%
% The input parameters are:
% t          = time (in epoch)
% Ex1, Ex2   = Two signals separated by distance along X axis
% Ey1, Ey2   = Two signals separated by distance along Y axis
% steplength = size of each FFT
% avnumber   = number of spectra to average over
% threshold  = degree of polarization above which we consider the
%              signal to be polarized (values between 0 and 1).
% plotflag   = set equal to 1 if we want the phase difference
%              and the coherence, instead of the usual output.
% distance   = distance between signals
%
% By Anders Tjulin (Last update 17/10-2002)
%

%% Check the input
warning off
if nargin>10
  error('Too many arguments')
end
if nargin<10,  distance=44;                end
if nargin<9,   plotflag=1;                 end
if nargin<8,   threshold=0.7;              end
if nargin<7,   avnumber=3;                 end
if nargin<6,   steplength=512;             end
if nargin<5,   error('Too few arguments'); end

Ex1=Ex1(:);
Ex2=Ex2(:);
Ey1=Ey1(:);
Ey2=Ey2(:);

%% Check the sampling rate

sampl=1/(t(2)-t(1));
%  t0=floor(t(1));
t0=t(1);

%% Create the Hamming window function

j=1:steplength;
w=0.54-0.46*cos(2*pi*j/(steplength-1));
normalisation=w*w'/steplength;

%% Calculate the spectral matrix with FFT

n=steplength-1;
nd2=ceil(n/2);
for j=1:steplength/2:length(Ex1)-steplength
  Yx1=fft(detrend(Ex1(j:j+steplength-1)).*w.');
  Yx2=fft(detrend(Ex2(j:j+steplength-1)).*w.');
  Yy1=fft(detrend(Ey1(j:j+steplength-1)).*w.');
  Yy2=fft(detrend(Ey2(j:j+steplength-1)).*w.');
  %    Yx=fft(Ex(j:j+steplength-1).*w.');
  %    Yy=fft(Ey(j:j+steplength-1).*w.');
  Yx1(1)=[];Yx2(1)=[];
  Yy1(1)=[];Yy2(1)=[];
  fEx1(:,2*((j-1)/steplength+1)-1)=Yx1(1:nd2)/n;
  fEx2(:,2*((j-1)/steplength+1)-1)=Yx2(1:nd2)/n;
  fEy1(:,2*((j-1)/steplength+1)-1)=Yy1(1:nd2)/n;
  fEy2(:,2*((j-1)/steplength+1)-1)=Yy2(1:nd2)/n;
  T(2*((j-1)/steplength+1)-1)=t(j+steplength/2);
end

%% Find the frequencies

nyq=1/2;
freq=sampl*(1:nd2)/(nd2)*nyq;

%% Calculate the polarisation parameters

I1=real(fEx1.*conj(fEx1)+fEx2.*conj(fEx2))/normalisation;
Q1=real(fEx1.*conj(fEx1)-fEx2.*conj(fEx2))/normalisation;
U1=2*real(fEx1.*conj(fEx2))/normalisation;
V1=2*imag(fEx1.*conj(fEx2))/normalisation;

I2=real(fEy1.*conj(fEy1)+fEy2.*conj(fEy2))/normalisation;
Q2=real(fEy1.*conj(fEy1)-fEy2.*conj(fEy2))/normalisation;
U2=2*real(fEy1.*conj(fEy2))/normalisation;
V2=2*imag(fEy1.*conj(fEy2))/normalisation;

%% Create mean values

Itemp1=0;Qtemp1=0;Utemp1=0;Vtemp1=0;
Itemp2=0;Qtemp2=0;Utemp2=0;Vtemp2=0;
for j=1:avnumber
  Itemp1=Itemp1+[zeros(steplength/2,j),I1,zeros(steplength/2,avnumber-j)];
  Qtemp1=Qtemp1+[zeros(steplength/2,j),Q1,zeros(steplength/2,avnumber-j)];
  Utemp1=Utemp1+[zeros(steplength/2,j),U1,zeros(steplength/2,avnumber-j)];
  Vtemp1=Vtemp1+[zeros(steplength/2,j),V1,zeros(steplength/2,avnumber-j)];
  Itemp2=Itemp2+[zeros(steplength/2,j),I2,zeros(steplength/2,avnumber-j)];
  Qtemp2=Qtemp2+[zeros(steplength/2,j),Q2,zeros(steplength/2,avnumber-j)];
  Utemp2=Utemp2+[zeros(steplength/2,j),U2,zeros(steplength/2,avnumber-j)];
  Vtemp2=Vtemp2+[zeros(steplength/2,j),V2,zeros(steplength/2,avnumber-j)];
end
Imean1=Itemp1(:,1:size(fEx1,2))/avnumber;Imean1=Imean1(:,avnumber:end);
Qmean1=Qtemp1(:,1:size(fEx1,2))/avnumber;Qmean1=Qmean1(:,avnumber:end);
Umean1=Utemp1(:,1:size(fEx1,2))/avnumber;Umean1=Umean1(:,avnumber:end);
Vmean1=Vtemp1(:,1:size(fEx1,2))/avnumber;Vmean1=Vmean1(:,avnumber:end);
Imean2=Itemp2(:,1:size(fEy1,2))/avnumber;Imean2=Imean2(:,avnumber:end);
Qmean2=Qtemp2(:,1:size(fEy1,2))/avnumber;Qmean2=Qmean2(:,avnumber:end);
Umean2=Utemp2(:,1:size(fEy1,2))/avnumber;Umean2=Umean2(:,avnumber:end);
Vmean2=Vtemp2(:,1:size(fEy1,2))/avnumber;Vmean2=Vmean2(:,avnumber:end);
T=T(avnumber:end);

root1=sqrt(Qmean1.*Qmean1+Umean1.*Umean1+Vmean1.*Vmean1);
p1=root1./Imean1;
vp1=Vmean1./root1;
phi1=atan2(Umean1,Qmean1)*180/pi;
coh1=(Umean1.^2+Vmean1.^2)./(Imean1.^2-Qmean1.^2);
phase1=atan2(Vmean1,Umean1)*180/pi;

root2=sqrt(Qmean2.*Qmean2+Umean2.*Umean2+Vmean2.*Vmean2);
p2=root2./Imean2;
vp2=Vmean2./root2;
phi2=atan2(Umean2,Qmean2)*180/pi;
coh2=(Umean2.^2+Vmean2.^2)./(Imean2.^2-Qmean2.^2);
phase2=atan2(Vmean2,Umean2)*180/pi;

%% Plot everything

if plotflag == 1 %%% The case with phase difference and coherence
  % 1. Time series
  % 2. Power spectra
  % 3. Coherence
  % 4. Phase Ex1,Ex2
  % 5. Phase Ey1,Ey2
  %
  temp1=find(coh1<threshold);phase_all1=phase1;phase1(temp1)=NaN;
  temp2=find(coh2<threshold);phase_all2=phase2;phase2(temp2)=NaN;
  
  clf;
  h(1)=irf_subplot(5,1,-1);
  plot(t-t0,Ex1,t-t0,Ex2);colorbar;axpos=get(gca,'position');
  plot(t-t0,Ex1,t-t0,Ex2);set(gca,'position',axpos);
  grid on;hold on;
  xlabel('UT')
  ylabel('mV/m')
  title([' nfft=' num2str(steplength) ', av=' num2str(avnumber)] )
  ht=irf_pl_info(['polarplot2() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],gca,[0,1 ]); set(ht,'interpreter','none');
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  
  h(2)=irf_subplot(5,1,-2);
  pcolor(T-t0,freq,log10(2*n*Imean1/sampl))
  shading flat
  colorbar('vert');
  ht=text(0,0,'Power density');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  ylabel('Hz')
  %    title('Power spectral density')
  %    irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  
  h(3)=irf_subplot(5,1,-3);
  pcolor(T-t0,freq,coh1)
  shading flat
  caxis([0,1])
  colorbar('vert')
  ylabel('Hz')
  %    title('Coherence')
  ht=text(0,0,'Coherence');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  
  h(4)=irf_subplot(5,1,-4);
  hp=pcolor(T-t0,freq,phase1);
  shading flat
  caxis([-180,180])
  colorbar('vert')
  ylabel('Hz')
  ht=text(0,0,'\Delta \phi_1');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  %    title('Phase difference')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  set(hp,    'buttondownfcn', 'polarplot2_manager(''ja'')');
  
  h(5)=irf_subplot(5,1,-5);
  hp=pcolor(T-t0,freq,phase2);
  shading flat
  caxis([-180,180])
  colorbar('vert')
  ylabel('Hz')
  ht=text(0,0,'\Delta \phi_2');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  %    title('Phase difference')
  temp1=axis;
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  set(hp,    'buttondownfcn', 'polarplot2_manager(''ja'')');
  
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
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  
  h(2)=subplot(4,1,2);
  pcolor(T-t0,freq,p)
  shading flat
  caxis([0,1])
  colorbar('vert')
  ylabel('Hz')
  title('Degree of polarisation')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
  %set(gca,'yscale','log')
  
  h(3)=subplot(4,1,3);
  pcolor(T-t0,freq,vp)
  shading flat
  caxis([-1,1])
  colorbar('vert')
  ylabel('Hz')
  title('Ellipticity')
  irf_timeaxis(gca,t0); set(gca,'tickdir','out');
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

%  suptitle([datestring(fromepoch(t0)) ', nfft=' num2str(steplength) ', av=' num2str(avnumber)] )
warning on

irf_zoom([0 t(end)-t(1)],'x',h);irf_timeaxis(h,t0);
global POLARPLOT_RESULTS
POLARPLOT_RESULTS={t0, T-t0,freq,2*n*Imean1/sampl,coh1,phase_all1,2*n*Imean2/sampl,coh2,phase_all2,threshold};


