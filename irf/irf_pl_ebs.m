function h=irf_pl_ebs(e,b,B,parameters)
%IRF_PL_EBS   Plot E, B, E/B, Poynting flux, E/B  wavelet spectra
%
% h=irf_pl_ebs(e,b,B,parameters)
% Plot E, B, E/B, Poynting flux, E/B  wavelet spectra
% developed from Anders Tjulins DUMWAVEPLOT
% assumes equidistant time spacing
%
% It uses a Morlet wavelet.
% e = wave electric field, columns (t ex ey ez)
% b = wave magnetic field, columns (t bx by bz)
% B = background magnetic field, columns (t Bx By Bz)
% parameters=[spec_width...
%    freq_min freq_max freq_n Morlet_width detrend colorbar plot_type panel]
%  spec -  Maximum pixel width of spectrogram (0 for all points) [CURRENTLY IGNORED]
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
% h=IRF_PL_EBS(e,b,B,parameters)
% h=IRF_PL_EBS(e,b,B)   % specify manually parameters
% Returns in h handles to subplots
%
% Example:
%    h=IRF_PL_EBS(e,b,B,[1000 2 180 50 5.36 1 1 2 0]);
%
% $Id$

% Modified by Andris Vaivads 15 April 2003
% By Anders Tjulin (Last update 15/4-2003)
%

% parameter definitions
persistent freq_int freq_number Morlet_width q_colormap plot_type colorbar_scale panel q_detrend q_spectra_width

if nargin == 4
  q_spectra_width=parameters(1);
  freq_int=[parameters(2) parameters(3)];
  freq_number=parameters(4);
  Morlet_width=parameters(5);
  q_detrend=parameters(6);if q_detrend==1, q_detrend='y';end
  colorbar_scale=parameters(7);
  plot_type=parameters(8);
  panel=parameters(9);
  if plot_type == 0
    if panel == 1, plot_param='e'; end
    if panel == 2, plot_param='b'; end
    if panel == 3, plot_param='s'; end
    if panel == 4, plot_param='eb'; end
  end
elseif nargin == 3
  q_spectra_width=irf_ask('Maximum pixel width of spectrogram (0 for all points) [%]>','q_spectra_width',1000);
  freq_int=irf_ask('Frequency interval to analyze [fmin fmax]? [%]>','freq_int',[2 10]);
  freq_number=irf_ask('Number of frequency steps? [%]>','freq_number',100);
  Morlet_width=irf_ask('The width of the Morlet wavelet (original 5.36)? [%]>','Morlet_width',5.36);
  q_detrend=irf_ask('Shall I detrend the data y/n? [%]>','q_detrend','y');
  colorbar_scale=irf_ask('Colorbar 1) green to red with the white in the middle 2) default ''jet''. [%]','colorbar_scale',2);
  plot_type=irf_ask('Plot: \n 0) one panel with ... \n 1) spectra e/b/Spar/e2b \n 2) time series e/b/B and spectra e/b/Spar/e2b \n 3) spectra Ex/Ey/Ez/Etot \n [%]>','plot_type',1);
  if plot_type == 0
    plot_param=irf_ask(' s) Poyn flux \n e) E spectra \n b) B spectra \n eb) E/B spectra \n [%]>','plot_param','s');
  end
else
  help irf_pl_ebs;return;
end


sampl_e=1/(e(2,1)-e(1,1));
sampl_b=1/(b(2,1)-b(1,1));
if     sampl_b > 1.5*sampl_e, e=irf_resamp(e,b); sampl=sampl_b; disp('irf_pl_ebs: interpolating e to b');
elseif sampl_e > 1.5*sampl_b, b=irf_resamp(b,e); sampl=sampl_e; disp('irf_pl_ebs: interpolating b to e');
elseif sampl_e == sampl_b && size(e)==size(b),   sampl=sampl_e;
else, sampl=2*sampl_e;
  t=max(e(1,1),b(1,1)):1/sampl:min(e(end,1),b(end,1)); t=t';
  e=irf_resamp(e,t); b=irf_resamp(b,t);
  irf_log('proc','interpolating b and e to 2x e sampling');
end
%% Check the sampling rate
disp(['Fs=' num2str(sampl) 'Fs_e=' num2str(sampl_e) 'Fs_b=' num2str(sampl_b)]);

%% Remove the last sample if the total number of samples is odd

if size(e,1)/2 ~= floor(size(e,1)/2)
  e=e(1:end-1,:);
  b=b(1:end-1,:);
end

% set to zero NaNs
ind_nan_e=isnan(e); e(ind_nan_e)=0;
ind_nan_b=isnan(b); b(ind_nan_b)=0;
ind_nan_B=isnan(B); B(ind_nan_B)=0;


%% the direction of background magnetic field
bn=irf_norm(irf_resamp(B,e));
bn(:,1) = [];
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
if strcmp(q_detrend,'y')
  e(:,2:4)=detrend(e(:,2:4));
  b(:,2:4)=detrend(b(:,2:4));
end

disp('irf_pl_ebs ... calculate e and b wavelet transform ....');
%% Make the FFT of all data
Swe=fft(e(:,2:4),[],1);
Swb=fft(b(:,2:4),[],1);

%% Get the correct frequencies for the wavelet transform
newfreq=w0./a;

%% Loop through all frequencies
ndata = size(e,1); nfreq = length(a);
powerEx_plot = zeros(ndata,nfreq);
powerEy_plot = zeros(ndata,nfreq);
powerEz_plot = zeros(ndata,nfreq);
power2E_plot = zeros(ndata,nfreq);
powerBx_plot = zeros(ndata,nfreq);
powerBy_plot = zeros(ndata,nfreq);
powerBz_plot = zeros(ndata,nfreq);
power2B_plot = zeros(ndata,nfreq);
Spar_plot = zeros(ndata,nfreq);
parfor ind_a=1:length(a)
  % if debug, disp([num2str(ind_a) '. frequency, ' num2str(newfreq(ind_a)) ' Hz.']);end
  mWexp = exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  mWexp = repmat(mWexp,1,3);
  Wwe = sqrt(1).*Swe.*mWexp;
  Wwb = sqrt(1).*Swb.*mWexp;
  
  %% Get the wavelet transform by IFFT of the FFT
  We = ifft(Wwe,[],1);
  Wb = ifft(Wwb,[],1);
  
  %% Calculate the power spectrum
  newfreqmat=w0/a(ind_a);
  
  %  power=(2*pi)*conj(W).*W./newfreqmat;
  powerE = 2*pi*(We.*conj(We))./newfreqmat;
  powerE(:,4) = sum(powerE,2);
  powerB = 2*pi*(Wb.*conj(Wb))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
  
  %% Poynting flux calculations, assume E and b units mV/m and nT, get  S in uW/m^2
  coef_poynt=10/4/pi*(1/4)*(4*pi); % 4pi from wavelets, see A. Tjulins power estimates a few lines above
  S = zeros(ndata,3);
  Wex=We(:,1);Wey=We(:,2);Wez=We(:,3);
  Wbx=Wb(:,1);Wby=Wb(:,2);Wbz=Wb(:,3);
  S(:,1)= coef_poynt*real(Wey.*conj(Wbz)+conj(Wey).*Wbz-Wez.*conj(Wby)-conj(Wez).*Wby)./newfreqmat;
  S(:,2)= coef_poynt*real(Wez.*conj(Wbx)+conj(Wez).*Wbx-Wex.*conj(Wbz)-conj(Wex).*Wbz)./newfreqmat;
  S(:,3)= coef_poynt*real(Wex.*conj(Wby)+conj(Wex).*Wby-Wey.*conj(Wbx)-conj(Wey).*Wbx)./newfreqmat;
  
  %For some reason this code works 20% slower than the above one
  %conjWe = conj(We); conjWb = conj(Wb);
  %S = coef_poynt*real( ...
  %    + We(:,[2 3 1]).*conjWb(:,[3 1 2]) + conjWe(:,[2 3 1]).*Wb(:,[3 1 2])...
  %    - We(:,[3 1 2]).*conjWb(:,[2 3 1]) - conjWe(:,[3 1 2]).*Wb(:,[2 3 1])...
  %    )./newfreqmat;
  Spar=sum(S.*bn,2);
  
  %% Remove data possibly influenced by edge effects
  censur=floor(2*a);
  censur_indexes=[1:min(censur(ind_a),size(e,1)) max(1,size(e,1)-censur(ind_a)):size(e,1)];
  powerE(censur_indexes,:) = NaN;
  powerB(censur_indexes,:) = NaN;
  Spar(censur_indexes) = NaN;
  
  powerEx_plot(:,ind_a) = powerE(:,1);
  powerEy_plot(:,ind_a) = powerE(:,2);
  powerEz_plot(:,ind_a) = powerE(:,3);
  power2E_plot(:,ind_a) = powerE(:,4);
  powerBx_plot(:,ind_a) = powerB(:,1);
  powerBy_plot(:,ind_a) = powerB(:,2);
  powerBz_plot(:,ind_a) = powerB(:,3);
  power2B_plot(:,ind_a) = powerB(:,4);
  Spar_plot(:,ind_a) = Spar;
end
idx_nan_e = sum(ind_nan_e,2)>0;
powerEx_plot(idx_nan_e,:) = NaN;
powerEy_plot(idx_nan_e,:) = NaN;
powerEz_plot(idx_nan_e,:) = NaN;
power2E_plot(idx_nan_e,:) = NaN;
Spar_plot(idx_nan_e,:) = NaN;
EtoB_plot=sqrt(power2E_plot./power2B_plot);

% if q_spectra_width,
%     error('q_spectra_width not impelemted')
%     if q_spectra_width<length(e(:,1)),
%         nav=ceil(length(e(:,1))/q_spectra_width);
%         t=e(1:nav:end,1);
%         for jj=1:length(t)-1
%             power2E_plot(jj,ind_a)=sum(power2E((jj-1)*nav+1:jj*nav))/nav;
%             power2B_plot(jj,ind_a)=sum(power2B((jj-1)*nav+1:jj*nav))/nav;
%             Spar_plot(jj,ind_a)=sum(Spar((jj-1)*nav+1:jj*nav))/nav;
%         end
%         last_point=length(t);n_points=length(power2E);
%         last_interval=(last_point-1)*nav+1:max(jj*nav,n_points);
%         n_in_last_interval=length(last_interval);
%         power2E_plot(last_point,ind_a)=sum(power2E(last_interval))/n_in_last_interval;
%         power2B_plot(last_point,ind_a)=sum(power2B(last_interval))/n_in_last_interval;
%         Spar_plot(last_point,ind_a)   =sum(Spar(last_interval))/n_in_last_interval;
%     end
% end

%% Plot everything
clear h;
switch plot_type
  case 0,   npl=NaN;
  case 1,   npl=4;clf;
  case 2,   npl=7;clf;
  case 3,   npl=4;clf;
end
ipl=1;

if colorbar_scale==1
  it=0:.02:1;it=it'; xcm=[ [0*it flipud(it) it];[it it it*0+1];[it*0+1 flipud(it) flipud(it)]; [flipud(it) 0*it 0*it]]; clear it;
else
  colormap('default');xcm=colormap;
end

if plot_type == 2
  %%%%%% E time series %%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(e);ylabel('E_{wave} [mV/m]');
  set(gca,'tickdir','out');colorbar;hh=get(gca,'position');colorbar off;set(gca,'position',hh); set(gca,'tickdir','in'); % to get the same size as colorbar plots
  %%%%%% B time series %%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(b);ylabel('B_{wave} [nT]');
  set(gca,'tickdir','out');colorbar;hh=get(gca,'position');colorbar off;set(gca,'position',hh); set(gca,'tickdir','in'); % to get the same size as colorbar plots
  %%%%%% Bo time series %%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(B);ylabel('B [nT]');
  set(gca,'tickdir','out');colorbar;hh=get(gca,'position');colorbar off;set(gca,'position',hh); set(gca,'tickdir','in'); % to get the same size as colorbar plots
  
  irf_timeaxis(h(1:3),'nolabels');
  irf_timeaxis(gca,'nodate');
end

if plot_type == 1 || plot_type == 2 || plot_type == 0
  t_start_epoch=get_t_start_epoch(t(1,1));
  %%%%%%%%% E spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 || (plot_type == 0 && strcmp(plot_param,'e'))
    %    pcolor(t-t0,newfreq,log10(abs(power2E.'))) % With edge effects removed
    pcolor(t-t_start_epoch,newfreq,log10(abs(power2E_plot.'))) % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    set(gca,'yscale','log','tickdir','out');
    cmean=nanmean(nanmean(log10(abs(power2E_plot))));
    caxis(floor(cmean)+[-3.5 3.5]);
    colormap(xcm);
    hca = colorbar;
    ylabel(hca,'E [(mV/m)^2/Hz]');
  end
  %%%%%%%%% B spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 || (plot_type == 0 && strcmp(plot_param,'b'))
    pcolor(t-t_start_epoch,newfreq,log10(abs(power2B_plot.'))) % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    %ht=text(0,0,'B [nT^2/Hz]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log','tickdir','out');
    cmean=nanmean(nanmean(log10(abs(power2B_plot))));
    caxis(floor(cmean)+[-3.5 3.5]);
    colormap(xcm);
    hca = colorbar;
    ylabel(hca,'B [nT^2/Hz]');
    %    irf_timeaxis(h,t0); % For time in epoch
  end
  %%%%%%%%% S spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1; end
  if plot_type ~= 0 || (plot_type == 0 && strcmp(plot_param,'s'))
    pcolor(t-t_start_epoch,newfreq,(sign(Spar_plot).*sqrt(abs(Spar_plot))).') % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    %ht=text(0,0,'S_{II} [\mu W/m^2Hz]^{1/2}');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log');set(gca,'tickdir','out');
    
    cc = [-max(max(sqrt(abs(Spar_plot)))) max(max(sqrt(abs(Spar_plot))))];
    caxis(cc);
    colormap(xcm);
    hca = colorbar;
    ylabel(hca,'S_{II} [\mu W/m^2Hz]^{1/2}');
  end
  %%%%%%%%% E/B spectra %%%%%%%%%%%%
  if plot_type ~= 0, h(ipl)=irf_subplot(npl,1,-ipl); end
  if plot_type ~= 0 || (plot_type == 0 && strcmp(plot_param,'eb'))
    pcolor(t-t_start_epoch,newfreq,log10(abs(EtoB_plot.'))) % With edge effects removed
    shading flat
    ylabel('f [Hz]')
    %ht=text(0,0,'log10(E/B) [(1000 km/s)]');set(ht,'units','normalized','position',[1 0.5],'rotation',90,'verticalalignment','top','horizontalalignment','center')
    set(gca,'yscale','log');set(gca,'tickdir','out');
    caxis([-1 4]);
    colormap(xcm);
    hca = colorbar;
    ylabel(hca,'log10(E/B) [(1e3 km/s)]');
  end
elseif plot_type == 3
  t_start_epoch=get_t_start_epoch(t(1,1));
  %%%%%%%%% Ex spectra %%%%%%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t_start_epoch,newfreq,log10(abs(powerEx_plot.')))
  shading flat
  ylabel('f [Hz]')
  set(gca,'yscale','log','tickdir','out');
  cmean=nanmean(nanmean(log10(abs(powerEx_plot))));
  caxis(floor(cmean)+[-3.5 3.5]);
  hca = colorbar;
  ylabel(hca,'Ex [(mV/m)^2/Hz]');
  %%%%%%%%% Ey spectra %%%%%%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t_start_epoch,newfreq,log10(abs(powerEy_plot.')))
  shading flat
  ylabel('f [Hz]')
  set(gca,'yscale','log','tickdir','out');
  cmean=nanmean(nanmean(log10(abs(powerEy_plot))));
  caxis(floor(cmean)+[-3.5 3.5]);
  hca = colorbar;
  ylabel(hca,'Ey [(mV/m)^2/Hz]');
  %%%%%%%%% Ez spectra %%%%%%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t_start_epoch,newfreq,log10(abs(powerEz_plot.')))
  shading flat
  ylabel('f [Hz]')
  set(gca,'yscale','log','tickdir','out');
  cmean=nanmean(nanmean(log10(abs(powerEz_plot))));
  caxis(floor(cmean)+[-3.5 3.5]);
  hca = colorbar;
  ylabel(hca,'Ez [(mV/m)^2/Hz]');
  %%%%%%%%% E spectra %%%%%%%%%%%%
  h(ipl)=irf_subplot(npl,1,-ipl);
  pcolor(t-t_start_epoch,newfreq,log10(abs(power2E_plot.')))
  shading flat
  ylabel('f [Hz]')
  set(gca,'yscale','log','tickdir','out');
  cmean=nanmean(nanmean(log10(abs(power2E_plot))));
  caxis(floor(cmean)+[-3.5 3.5]);
  hca = colorbar;
  ylabel(hca,'E [(mV/m)^2/Hz]');
  
  irf_timeaxis(h);
end

if plot_type ~=0
  axes(h(1));
  title(['Width Morlet wavelet = ' num2str(Morlet_width)]);
  %ht=irf_pl_info([mfilename '  ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))]); set(ht,'interpreter','none'); % add information to the plot
  irf_zoom(h,'x',[min(t) max(t)]);
  irf_timeaxis(h)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_start_epoch=get_t_start_epoch(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives back the value of t_start_epoch of the figure
% if not  set, sets t_start_epoch of the figure
ud=get(gcf,'userdata');
ii = find(~isnan(t));
if ii
  valid_time_stamp=t(ii(1));
else
  valid_time_stamp=[];
end

if isfield(ud,'t_start_epoch')
  t_start_epoch=ud.t_start_epoch;
elseif valid_time_stamp
  if valid_time_stamp > 1e8 % set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
    t_start_epoch=valid_time_stamp;
    ud.t_start_epoch=t_start_epoch;
    set(gcf,'userdata',ud);
    irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
  else
    t_start_epoch=0;
  end
else
  t_start_epoch=0;
end
end

function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   Revision: 1.1.8.1   Date: 2010/03/16 00:15:50

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
  % Count up non-NaNs.
  n = sum(~nans);
  n(n==0) = NaN; % prevent divideByZero warnings
  % Sum up non-NaNs, and divide by the number of non-NaNs.
  m = sum(x) ./ n;
else
  % Count up non-NaNs.
  n = sum(~nans,dim);
  n(n==0) = NaN; % prevent divideByZero warnings
  % Sum up non-NaNs, and divide by the number of non-NaNs.
  m = sum(x,dim) ./ n;
end
end

