function h=av_ebBplot(e,b,B,parameters)

%AV_ebBplot(e,b,B,parameters) plot of your choice
% E, B, E/B, Poynting flux, E/B  wavelet spectra
% developed from Anders Tjulins DUMWAVEPLOT
% assumes equidistant time spacing
%
% It uses a Morlet wavelet.
% e = wave electric field, columns (t ex ey ez)
% b = wave magnetic field, columns (t bx by bz)
% B = background magnetic field, columns (t Bx By Bz)
% parameters=[spec_width freq_min freq_max freq_n Morlet_width detrend colorbar plot_type panel]
%  spec -  Maximum pixel width of spectrogram (0 for all points)
%  freq_min, freq_max  -     Frequency interval to analyze
%  Morlet_width -  The width of the Morlet wavelet (original 5.36)
%  detrend  - 1) do linear detrend, 0) do nothing
%  colorbar -  1) green to red with white in the middle  2) default 'jet'
%  plot_type - 0) one panel with ...
%              1) spectra e/b/Spar/e2b
%              2) time series e/b/B and spectra e/b/Spar/e2b
%              3) spectra Ex/Ey/Ez/Etot
%  panel -  1) E spectra                (necessary only if plot_type==0)
%           2) B spectra
%           3) Poyn flux
%           4) E/B spectra
%
% h=AV_WAVEPLOT_ebB(e,b,B,parameters)
% h=AV_WAVEPLOT_ebB(e,b,B)   % specify manually parameters
% Returns in h handles to subplots
%
% Example h=av_waveplot_ebB(e,b,B,[1000 2 180 50 5.36 1 1 2 0]);
% Modified by Andris Vaivads 15 April 2003
% By Anders Tjulin (Last update 15/4-2003)
%

% parameter definitions
global AV_DEBUG; if isempty(AV_DEBUG), debug=0;else, debug=AV_DEBUG;end
persistent freq_int freq_number Morlet_width q_colormap plot_type colorbar_scale panel q_detrend

if nargin == 4,
  q_spectra_width=parameters(1);
  freq_int=[parameters(2) parameters(3)];
  freq_number=parameters(4);
  Morlet_width=parameters(5);
  q_detrend=parameters(6);if q_detrend==1, q_detrend='y';end
  colorbar_scale=parameters(7);
  plot_type=parameters(8);
  panel=parameters(9);
  if plot_type == 0,
    if panel == 1, plot_param='e'; end
    if panel == 2, plot_param='b'; end
    if panel == 3, plot_param='s'; end
    if panel == 4, plot_param='eb'; end
  end
elseif nargin == 3,
  q_spectra_width=av_q('Maximum pixel width of spectrogram (0 for all points) [%]>','q_spectra_width',1000);
  freq_int=av_q('Frequency interval to analyze [fmin fmax]? [%]>','freq_int',[2 10]);
  freq_number=av_q('Number of frequency steps? [%]>','freq_number',100);
  Morlet_width=av_q('The width of the Morlet wavelet (original 5.36)? [%]>','Morlet_width',5.36);
  q_detrend=av_q('Shall I detrend the data y/n? [%]>','q_detrend','y');
  colorbar_scale=av_q('Colorbar 1) green to red with the white in the middle 2) default ''jet''. [%]','colorbar_scale',2);
  plot_type=av_q('Plot: \n 0) one panel with ... \n 1) spectra e/b/Spar/e2b \n 2) time series e/b/B and spectra e/b/Spar/e2b \n 3) spectra Ex/Ey/Ez/Etot \n [%]>','plot_type',1);
  if plot_type == 0,
    plot_param=av_q(' s) Poyn flux \n e) E spectra \n b) B spectra \n eb) E/B spectra \n [%]>','plot_param','s');
  end
else
  help av_waveplot_ebB;return;
end


  sampl_e=1/(e(2,1)-e(1,1));
  sampl_b=1/(b(2,1)-b(1,1));
  if     sampl_b > 1.5*sampl_e, e=av_interp(e,b); sampl=sampl_b; disp('av_ebBplot: interpolating e to b');
  elseif sampl_e > 1.5*sampl_b, b=av_interp(b,e); sampl=sampl_e; disp('av_ebBplot: interpolating b to e');
  else   sampl=2*sampl_e; t=max(e(1,1),b(1,1)):1/sampl:min(e(end,1),b(end,1)); t=t'; e=av_interp(e,t); b=av_interp(b,t); disp('av_ebBplot: interpolating b and e to 2x e sampling');
  end
  %% Check the sampling rate
  disp(['Fs=' num2str(sampl) 'Fs_e=' num2str(sampl_e) 'Fs_b=' num2str(sampl_b)]);
  t0=e(1,1);

  %% Remove the last sample if the total number of samples is odd

  if size(e,1)/2 ~= floor(size(e,1)/2)
    e=e(1:end-1,:);
    b=b(1:end-1,:);
  end

%% the direction of background magnetic field
bn=av_norm(av_interp(B,e));
t=e(:,1);

  %% Find the frequencies for an FFT of all data

  nd2=size(e,1)/2;
  nyq=1/2;
  freq=sampl*(1:nd2)/(nd2)*nyq;
  w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

  %%
  %% Set some important parameters
  %%

  amin=log10(0.5*sampl/freq_int(2));amax=log10(0.5*sampl/freq_int(1));anumber=freq_number;
%  amin=0.01; % The highest frequency to consider is 0.5*sampl/10^amin
%  amax=2; % The lowest frequency to consider is 0.5*sampl/10^amax
%  anumber=400; % The number of frequencies
  a=logspace(amin,amax,anumber);
%  a=logspace(0.01,2.4,100);
w0=sampl/2; % The maximum frequency
%  sigma=5.36/w0; % The width of the Morlet wavelet
sigma=Morlet_width/w0; % The width of the Morlet wavelet
if strcmp(q_detrend,'y'),
  e(:,2:4)=detrend(e(:,2:4));
  b(:,2:4)=detrend(b(:,2:4));
end

disp('av_ebBplot ... calculate e and b wavelet transform ....');
%% Make the FFT of all data
Swex=fft(e(:,2));Swey=fft(e(:,3));Swez=fft(e(:,4));
Swbx=fft(b(:,2));Swby=fft(b(:,3));Swbz=fft(b(:,4));

%% Get the correct frequencies for the wavelet transform
  newfreq=w0./a; 
  newfreqmat=w0/a(ind_a);

  %% Loop through all frequencies 

for ind_a=1:length(a),
  if debug, disp([num2str(ind_a) '. frequency, ' num2str(newfreq(ind_a)) ' Hz.']);end
  Wwex=sqrt(1).*Swex.*exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  Wwey=sqrt(1).*Swey.*exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  Wwez=sqrt(1).*Swez.*exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  Wwbx=sqrt(1).*Swbx.*exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  Wwby=sqrt(1).*Swby.*exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  Wwbz=sqrt(1).*Swbz.*exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  %% Get the wavelet transform by IFFT of the FFT
  %W=ifft(Ww);
  Wex=ifft(Wwex);Wey=ifft(Wwey);Wez=ifft(Wwez);
  Wbx=ifft(Wwbx);Wby=ifft(Wwby);Wbz=ifft(Wwbz);
  
  %  clear Wwex Wwey Wwez Wwbx Wwby Wwbz;
  % disp('av_waveplot_ebB ... calculate all spectra ....');
  
  % Poynting flux calculations, assume E and b units mV/m and nT, get  S in uW/m^2
  coef_poynt=10/4/pi/2;
  Sx= real(Wey.*conj(Wbz)-Wez.*conj(Wby));
  Sy=-real(Wez.*conj(Wbx)-Wex.*conj(Wbz));
  Sz= real(Wex.*conj(Wby)-Wey.*conj(Wbx));
  
  Spar=Sx.*bn(:,2)+Sy.*bn(:,3)+Sz.*bn(:,4);
  
  %% Calculate the power spectrum
  
  %  power=(2*pi)*conj(W).*W./newfreqmat;
  powerEx=2*pi*(Wex.*conj(Wex))./newfreqmat;
  powerEy=2*pi*(Wey.*conj(Wey))./newfreqmat;
  powerEz=2*pi*(Wez.*conj(Wez))./newfreqmat;
  powerE=powerEx+powerEy+powerEz;
  
  powerBx=2*pi*(Wbx.*conj(Wbx))./newfreqmat;
  powerBy=2*pi*(Wby.*conj(Wby))./newfreqmat;
  powerBz=2*pi*(Wbz.*conj(Wbz))./newfreqmat;
  powerB=powerBx+powerBy+powerBz;
  
  %% Remove data possibly influenced by edge effects
  censur=floor(2*a);
  power2E=powerE;
  power2B=powerB;
  %  for j=1:anumber;
  censur_indexes=[1:min(censur(ind_a),size(e,1)) max(1,size(e,1)-censur(ind_a)):size(e,1)];
  power2E(censur_indexes)=NaN;
  power2B(censur_indexes)=NaN;
  Spar(censur_indexes)=NaN;
  
  if q_spectra_width,
    if q_spectra_width<length(e(:,1)),
      nav=ceil(length(e(:,1))/q_spectra_width);
      t=e(1:nav:end,1);
      for jj=1:length(t)-1
        power2E_plot(jj,ind_a)=sum(power2E((jj-1)*nav+1:jj*nav))/nav;
        power2B_plot(jj,ind_a)=sum(power2B((jj-1)*nav+1:jj*nav))/nav;
        Spar_plot(jj,ind_a)=sum(Spar((jj-1)*nav+1:jj*nav))/nav;
      end
      last_point=length(t);n_points=length(power2E);
      last_interval=(last_point-1)*nav+1:max(jj*nav,n_points);
      n_in_last_interval=length(last_interval);
      power2E_plot(last_point,ind_a)=sum(power2E(last_interval))/n_in_last_interval;
      power2B_plot(last_point,ind_a)=sum(power2B(last_interval))/n_in_last_interval;
      Spar_plot(last_point,ind_a)   =sum(Spar(last_interval))/n_in_last_interval;
    end
  else,
    power2E_plot(:,ind_a)=power2E;
    power2B_plot(:,ind_a)=power2B;
    Spar_plot(:,ind_a)=Spar;
  end
  EtoB_plot=sqrt(power2E_plot./power2B_plot);
  %% Plot everything
end
clear h;
switch plot_type
case 0,   npl=NaN;
case 1,   npl=4;clf;
case 2,   npl=7;clf;
case 3,   npl=4;clf;
end
ipl=1;

if colorbar_scale==1,
      it=0:.02:1;it=it'; xcm=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
else
      colormap('default');xcm=colormap;
end

if plot_type == 2,
%%%%%% E time series %%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(e);ylabel('E_{wave} [mV/m]');
hh=get(gca,'position');set(gca,'position',[hh(1) hh(2) hh(3)*.8550 hh(4)]); % to get the same size as colorbar plots
%  colorbar;
%%%%%% E time series %%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(b);ylabel('B_{wave} [nT]');
hh=get(gca,'position');set(gca,'position',[hh(1) hh(2) hh(3)*.8550 hh(4)]); % to get the same size as colorbar plots
%  colorbar;
%%%%%% E time series %%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(B);ylabel('B [nT]');
hh=get(gca,'position');set(gca,'position',[hh(1) hh(2) hh(3)*.8550 hh(4)]); % to get the same size as colorbar plots
%  colorbar;
  add_timeaxis(h(1:3),'nolabels');
  add_timeaxis(gca,'nodate');
end

if plot_type == 1 | plot_type == 2 | plot_type == 0,
%%%%%%%%% E spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 | (plot_type == 0 & strcmp(plot_param,'e')),
%    pcolor(t-t0,newfreq,log10(abs(power2E.'))) % With edge effects removed
    pcolor(t,newfreq,log10(abs(power2E_plot.'))) % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    ht=text(0,0,'E [(mV/m)^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log');set(gca,'tickdir','out');
    cmean=mean(mean(log10(abs(powerE))));
    caxis(floor(cmean)+[-3.5 3.5]);
%      caxis([-5 2]);
    colormap(xcm);colorbar
  end
%%%%%%%%% B spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 | (plot_type == 0 & strcmp(plot_param,'b')),
    pcolor(t,newfreq,log10(abs(power2B_plot.'))) % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    ht=text(0,0,'B [nT^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log');set(gca,'tickdir','out');
      cmean=mean(mean(log10(abs(powerB))));
      caxis(floor(cmean)+[-3.5 3.5]);
%      caxis([-8 -1]);
    colormap(xcm);colorbar;
%    add_timeaxis(h,t0); % For time in epoch
  end
%%%%%%%%% S spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 | (plot_type == 0 & strcmp(plot_param,'s')),
    pcolor(t,newfreq,(sign(Spar_plot).*sqrt(abs(Spar_plot))).') % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    ht=text(0,0,'S_{II} [\mu W/m^2Hz]^{1/2}');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log');set(gca,'tickdir','out');

    cc = [-max(max(sqrt(abs(Spar)))) max(max(sqrt(abs(Spar))))];
  	caxis(cc);
  	colormap(xcm);colorbar;
    if plot_type ~= 0, axes(h(ipl-2));colorbar; axes(h(ipl-3));colorbar; end
  end
%%%%%%%%% E/B spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 | (plot_type == 0 & strcmp(plot_param,'eb')),
    pcolor(t,newfreq,log10(abs(EtoB_plot.'))) % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    ht=text(0,0,'log10(E/B) [(1000 km/s)]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log');set(gca,'tickdir','out');
    caxis([-1 4]);
    colormap(xcm);colorbar
  end
elseif plot_type == 3,
%%%%%%%%% Ex spectra %%%%%%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t0,newfreq,log10(abs(powerEx.'))) % Without edge effects removed
  shading flat
  ylabel('f [Hz]')
  ht=text(0,0,'Ex [(mV/m)^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  set(gca,'yscale','log');set(gca,'tickdir','out');
  caxis([-5 2]);colorbar
%%%%%%%%% Ey spectra %%%%%%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t0,newfreq,log10(abs(powerEy.'))) % Without edge effects removed
  shading flat
  ylabel('f [Hz]')
  ht=text(0,0,'Ey [(mV/m)^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  set(gca,'yscale','log');set(gca,'tickdir','out');
  caxis([-5 2]);colorbar
%%%%%%%%% Ez spectra %%%%%%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t0,newfreq,log10(abs(powerEz.'))) % Without edge effects removed
  shading flat
  ylabel('f [Hz]')
  ht=text(0,0,'Ez [(mV/m)^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  set(gca,'yscale','log');set(gca,'tickdir','out');
  caxis([-5 2]);colorbar
%%%%%%%%% E spectra %%%%%%%%%%%%
  h(ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t0,newfreq,log10(abs(powerE.'))) % Without edge effects removed
  shading flat
  ylabel('f [Hz]')
  ht=text(0,0,'E [(mV/m)^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
  set(gca,'yscale','log');set(gca,'tickdir','out');
  caxis([-5 2]);colorbar

  add_timeaxis(h,t0);
end

if plot_type ~=0,
  axes(h(1));
  title(['Width Morlet wavelet = ' num2str(Morlet_width)]);
  ht=av_pl_info([mfilename '  ' datestr(now)]); set(ht,'interpreter','none'); % add information to the plot
  av_zoom([min(t) max(t)],'x',h);
  add_timeaxis(h)
end

