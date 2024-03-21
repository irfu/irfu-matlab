function h=irf_wave_detection_algorithm(tint,scId,varargin)
%IRF_WAVE_DETECTION_ALGORITHM  Find banded EMIC waves
%
% [t,newfreq,powerCrossCov_SM_plot,hCyclFreq,heCyclFreq,oCyclFreq,...
%     power_median_removed,waveFrequencies] = ...
%     IRF_WAVE_DETECTION_ALGORITHM(tint,varargin)
%
% IRF_WAVE_DETECTION_ALGORITHM duplicates the method of Bortnik et al. 2007
% to detect EMIC waves
%
% 	IRF_WAVE_DETECTION_ALGORITHM(b) Produces a plot showing magnetic field spectrogram in the upper panel
%    and the spectrogram with the background median removed in the lower
%    panel. Overplotted are points picking out wave events, and also lines
%    for the cyclotron frequency (H, solid) (He, dashed) (O dash-dot).
%   A file of the time table is exported.
%
%   The conditions for finding events are not always perfect.
%   It is important to confirm any event by eye!
%
%   Example:
%    irf_wave_detection_algorithm(tint,cl_id,'freq',[.02 5]);
%      OR (as used within c_ulf_process.m) give ebsp structure and
%      background B field
%    irf_wave_detection_algorithm2(ebsp,cl_id,bf);
%

save_plot=0;
Units = irf_units;
if isnumeric(scId), scId_s = sprintf('C%d',scId);
else, scId_s = ['TH' upper(scId)];
end
if isstruct(tint)
  ebsp=tint;
  powerCrossCov_SM_plot = ebsp.bb_xxyyzzss(:,:,4);
  t = ebsp.t;
  tint = [t(1) t(end)];
  deltaT=30;
  sampl=1/30;
  outTime = (t(1):deltaT:t(end))' + deltaT/2; outTime(end) = [];
  powerCrossCov_SM_plot = AverageData(powerCrossCov_SM_plot,t,outTime,deltaT);
  t = outTime;
  newfreq = ebsp.f;
  nfreq = length(newfreq);
  ndata = length(t);
  Btot = ebsp.B0;
  save_plot = 1;
else

  %% Check input
  [~,args,~] = axescheck(varargin{:});
  b=local.c_read(['B_vec_xyz_gse__' scId_s '_CP_FGM_5VPS'],tint);

  %% get background magnetic field
  bf=irf_filt(b,1/600,0,[],5);
  b0=b;
  b0(:,2:4)=b(:,2:4)-bf(:,2:4);
  B=b0;

  sampl=10;
  t=b(1,1):1/sampl:b(end,1); t=t';
  B=irf_resamp(B,t);
  b=irf_resamp(b,t); disp('resampling to 10 Hz');
  b(:,2:4)=b(:,2:4)-B(:,2:4);

  %% Remove the last sample if the total number of samples is odd

  if size(b,1)/2 ~= floor(size(b,1)/2)
    b=b(1:end-1,:);
    B=B(1:end-1,:);
    t=t(1:end-1,:);
  end


  % set to zero NaNs
  ind_nan_b=isnan(b); b(ind_nan_b)=0;
  ind_nan_B=isnan(B); B(ind_nan_B)=0;

  Btot = B(:,1:2);
  Btot(:,2)=sqrt(B(:,2).*B(:,2)+B(:,3).*B(:,3)+B(:,4).*B(:,4));

  %% Find the frequencies for an FFT of all data

  nd2=size(b,1)/2;
  nyq=1/2;
  freq=sampl*(1:nd2)/(nd2)*nyq;
  w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

  %% Set some important parameters
  freq_int=[.02 5];
  freq_number=21;
  Morlet_width=5.36;

  if isnumeric(args{end})
    freq_int=args{end};
  elseif ismember({'freq'},args)
    disp('frequency interval values missing. using default')
  end


  amin=log10(0.5*sampl/freq_int(2));amax=log10(0.5*sampl/freq_int(1));anumber=freq_number;
  a=logspace(amin,amax,anumber);
  w0=sampl/2; % The maximum frequency
  sigma=Morlet_width/w0; % The width of the Morlet wavelet

  Swb=fft(b(:,2:4),[],1);


  %% Get the correct frequencies for the wavelet transform
  newfreq=w0./a;
  ndata = size(b,1); nfreq = length(a);
  powerCrossCov_SM_plot = zeros(ndata,nfreq);

  parfor ind_a=1:length(a)
    mWexp = exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
    mWexp = repmat(mWexp,1,3);
    Wwb = sqrt(1).*Swb.*mWexp;

    %% Get the wavelet transform by IFFT of the FFT
    Wb = ifft(Wwb,[],1);

    %% Calculate the power spectrum
    newfreqmat=w0/a(ind_a);

    %% spectral matrix
    spectralMatrix = zeros(3,3,ndata);

    spectralMatrix(1,1,:) = 2*pi*(Wb(:,1).*conj(Wb(:,1)))./newfreqmat;
    spectralMatrix(1,2,:) = 2*pi*(Wb(:,1).*conj(Wb(:,2)))./newfreqmat;
    spectralMatrix(1,3,:) = 2*pi*(Wb(:,1).*conj(Wb(:,3)))./newfreqmat;
    spectralMatrix(2,1,:) = 2*pi*(Wb(:,2).*conj(Wb(:,1)))./newfreqmat;
    spectralMatrix(2,2,:) = 2*pi*(Wb(:,2).*conj(Wb(:,2)))./newfreqmat;
    spectralMatrix(2,3,:) = 2*pi*(Wb(:,2).*conj(Wb(:,3)))./newfreqmat;
    spectralMatrix(3,1,:) = 2*pi*(Wb(:,3).*conj(Wb(:,1)))./newfreqmat;
    spectralMatrix(3,2,:) = 2*pi*(Wb(:,3).*conj(Wb(:,2)))./newfreqmat;
    spectralMatrix(3,3,:) = 2*pi*(Wb(:,3).*conj(Wb(:,3)))./newfreqmat;

    SMpermute = permute(spectralMatrix,[3,1,2]);


    %% power
    powerCrossCov_SM_plot(:,ind_a) = SMpermute(:,3,3)+SMpermute(:,1,1)+SMpermute(:,2,2);
  end


  %% average down the signal


  %% smooth data without decreasing sampling frequency

  powerCrossCov_avg = zeros(ndata,nfreq);
  smoothwidth = 200;
  halfw = 100;
  for j=1:nfreq
    sumPower = nansum(powerCrossCov_SM_plot(1:smoothwidth,j));
    for i=1:ndata-smoothwidth
      powerCrossCov_avg(i+halfw-1,j) = sumPower;
      sumPower=sumPower-powerCrossCov_SM_plot(i,j);
      sumPower=sumPower+powerCrossCov_SM_plot(i+smoothwidth,j);
    end
  end

  powerCrossCov_SM_plot = powerCrossCov_avg./smoothwidth;

  %% put NaNs back in place

  idx_nan_b = sum(ind_nan_b,2)>0;
  powerCrossCov_SM_plot(idx_nan_b,:) = NaN;

  %% Remove data possibly influenced by edge effects
  for ind_a=1:length(a)
    censur=floor(2*a);
    censur_indexes=[1:min(censur(ind_a),size(b,1)) max(1,size(b,1)-censur(ind_a)):size(b,1)];
    SMpermute(censur_indexes,:,:) = NaN;
    powerCrossCov_SM_plot(censur_indexes,:,:) = NaN;
  end

  %% down-sample

  sampl1 = 1/30; %2 samples per minute
  t1 = b(1,1):1/sampl1:b(end,1); t1=t1';
  powerCrossCov_SM_plot = interp1(t,powerCrossCov_SM_plot,t1,'linear','extrap');
  sampl = sampl1;
  t = t1;
  ndata = size(powerCrossCov_SM_plot,1);

end

%% remove background median over a window
toa = .5; %time (in hours) of average for background median

toa = toa*3600;
toaind = toa*sampl;
halft = floor(toa/2);
halftind = floor(toaind/2);
totalt = t(end)-t(1);
totaltind = ceil(totalt*sampl);

if totalt <= halft
  medianPower = nanmedian(powerCrossCov_SM_plot);
  medianPower = repmat(medianPower,ndata,1);
elseif totalt <= toa && totalt > halft
  edgetind = totaltind-halftind;
  medianPower = zeros(size(powerCrossCov_SM_plot));
  medianPower(edgetind+1:ndata-edgetind,:) = repmat(nanmedian(powerCrossCov_SM_plot),ndata-2*edgetind,1);
  for i = 1:edgetind
    medianPower(i,:) = nanmedian(powerCrossCov_SM_plot(1:i+halftind,:));
  end
  for i = ndata-edgetind+1:ndata
    medianPower(i,:) = nanmedian(powerCrossCov_SM_plot(i-halftind:end,:));
  end
else
  medianPower = zeros(size(powerCrossCov_SM_plot));
  for i = 1:halftind
    medianPower(i,:) = nanmedian(powerCrossCov_SM_plot(1:i+halftind,:));
  end
  for i = halftind+1:ndata-halftind
    medianPower(i,:) = nanmedian(powerCrossCov_SM_plot(i-halftind:i+halftind,:));
  end
  for i = ndata-halftind+1:ndata
    medianPower(i,:) = nanmedian(powerCrossCov_SM_plot(i-halftind:end,:));
  end

end


%power_median_removed = log10(abs(powerCrossCov_SM_plot)) - log10(abs(medianPower));
power_median_removed = powerCrossCov_SM_plot - medianPower;

ind_nan_pmr=isnan(power_median_removed); power_median_removed(ind_nan_pmr)=0;


%% set cutoff frequencies and put a limit on the width of the emic event
magB = irf_resamp(Btot,t);
hCyclFreq  = magB(:,2).*1e-9 .* Units.e/(2*pi*Units.mp);
heCyclFreq = magB(:,2).*1e-9 .* Units.e/(2*pi*4*Units.mp);
oCyclFreq  = magB(:,2).*1e-9 .* Units.e/(2*pi*16*Units.mp);

%% wave detection
P=[0 0 0 0 0]; %array for peaks [time, freq, peakValue, lowerBound, upperBound]
peak = 1;
for i=1:ndata
  Y = power_median_removed(i,:);
  dY = zeros(size(Y));
  for j=1:nfreq-1
    dY(j) = (Y(j+1)-Y(j))/(newfreq(j+1)-newfreq(j));
  end
  for j=2:nfreq
    if sign(dY(j)) > sign(dY(j-1))
      peakValue = Y(j);
      peakFreq = newfreq(j);
      lowerBound = peakFreq;
      upperBound = peakFreq;
      k=0;
      eventThresh = (powerCrossCov_SM_plot(i,:)./medianPower(i,:)) > 4;
      while (j-k > 1) && eventThresh(j-k) == 1
        upperBound = newfreq((j-k));
        k=k+1;
      end
      k=0;
      while eventThresh(j+k) == 1 && j+k < nfreq
        lowerBound = newfreq((j+k));
        k=k+1;
        if j+k == nfreq
          lowerBound = upperBound;
        end
      end
      if j+k+1 < numel(eventThresh) && eventThresh(j+k+1) == 1
        lowerBound = upperBound;
      end
      if lowerBound == newfreq(end-1) && eventThresh(end) == 1
        lowerBound = upperBound;
      end

      if eventThresh(j) == 1 && upperBound < 4*lowerBound && ...
          upperBound > 1.25*lowerBound && peakFreq < hCyclFreq(i) ...
          && peakFreq ~= newfreq(end) && peakFreq > oCyclFreq(i)/4
        P(peak,:) = [t(i) peakFreq peakValue lowerBound upperBound];
        peak = peak+1;
      end
    end
  end

end

power_median_removed(ind_nan_pmr)=NaN;


%% wave association
%less than 10 minutes between consecutive points,
%wave activity detected 15% of the time (at least 18 points per hour)
%(or 20% of the time/ 24 points per hour),
%at least 5 points total

nPeaks = size(P,1);
waveEvent = [0 0]; %start time, end time
nevents=1;
counter=1;
try
  TTemic=irf.TimeTable([scId_s '_MAARBLE_PC12_wave_events']);
  createTTemic = 0;
catch
  createTTemic = 1;
end
if createTTemic
  TTemic = irf.TimeTable;
  TTemic.Header={[scId_s ' EMIC events for Maarble']};
end
while counter < nPeaks
  startTime = P(counter,1);
  lowFreq = P(counter,4);
  highFreq = P(counter,5);
  endTime = P(counter,1);
  npeaks=1;
  while counter < nPeaks && P(counter+1,1)-P(counter,1) < 10*60
    endTime = P(counter+1,1);
    if P(counter+1,4) < lowFreq
      lowFreq = P(counter+1,4);
    end
    if P(counter+1,5) > highFreq
      highFreq = P(counter+1,5);
    end
    counter=counter+1;
    npeaks=npeaks+1;
  end
  if (npeaks)/(endTime-startTime) > 24/3600 && npeaks > 4
    waveEvent(nevents,:) = [startTime endTime];
    timeInt = [startTime-300 endTime+300]; %include 5 minutes on either side for the time table
    TTemic=add(TTemic,timeInt,{num2str(lowFreq),num2str(highFreq)});
    nevents = nevents+1;
  else
    counter=counter+1;
  end

end

if waveEvent(1) > 0
  ascii(TTemic)
  TTemic=unique(TTemic);
  export_ascii(TTemic,[scId_s '_MAARBLE_PC12_wave_events'])


  %% plot

  figure(318), clf
  clear h;
  npl=2;clf;
  ipl=1;
  cmapSpace = irf_colormap('space');

  t_start_epoch=get_t_start_epoch(t(1,1));


  %%%%% B spectra from spectral matrix

  h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  pcolor(t-t_start_epoch,newfreq,log10(abs(powerCrossCov_SM_plot.')))
  shading flat
  ylabel('f [Hz]')
  set(gca,'yscale','log','tickdir','out');
  axis(axis)
  hold on
  plot(t-t_start_epoch,hCyclFreq,'color','k');
  plot(t-t_start_epoch,heCyclFreq,'--','color','k');
  plot(t-t_start_epoch,oCyclFreq,'-.','color','k');
  hold off
  caxis([-6.5 2.5]);
  colormap(cmapSpace);
  if isa(h(1),'handle'), hca2 = colorbar(h(1)); % HG2
  else, hca2 = colorbar('peer',h(1));
  end
  ylabel(hca2,{'log(B)'; '[nT^2/Hz]'});

  %%%%% B spectra with median removed

  h(ipl)=irf_subplot(npl,1,-ipl);
  pcolor(t-t_start_epoch,newfreq,log10(abs(power_median_removed.')))
  shading flat
  ylabel('f [Hz]')
  set(gca,'yscale','log','tickdir','out');
  axis(axis)
  hold on
  plot(P(:,1)-t_start_epoch,P(:,2),'d','markerfacecolor','k','markeredgecolor','k','markersize',6)
  plot(P(:,1)-t_start_epoch,P(:,4),'o','markerfacecolor','k','markeredgecolor','k','markersize',4)
  plot(P(:,1)-t_start_epoch,P(:,5),'o','markerfacecolor','k','markeredgecolor','k','markersize',4)

  plot(t-t_start_epoch,hCyclFreq,'color','k');
  plot(t-t_start_epoch,heCyclFreq,'--','color','k');
  plot(t-t_start_epoch,oCyclFreq,'-.','color','k');

  hold off
  caxis([-3.5 2.5]);
  colormap(cmapSpace);
  if isa(h(1),'handle'), hca2 = colorbar(h(1)); % HG2
  else, hca2 = colorbar('peer',h(1));
  end
  ylabel(hca2,{'log(B)'; '[nT^2/Hz]'});


  %% Add figure menu
  irf_figmenu;


  axes(h(1));
  if exist('timeInt','var')
    irf_zoom(h,'x',[timeInt(1)-600 timeInt(2)+600]);
  else
    irf_zoom(h,'x',[min(t) max(t)]);
  end
  xlimlast=get(h(end),'xlim');
  start_time = irf_time(xlimlast(1) + t_start_epoch,'vector');
  time_label = datestr( datenum(start_time),1 );

  irf_legend(h(1),[scId_s '           ' time_label],[0 1.05],'fontsize',10,'color','cluster');
  irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);

  if save_plot
    print('-dpng',[scId_s '_MAARBLE_PC12_wave_detection_' irf_fname(tint,5)])
  end
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


function y = nanmedian(x)
%NANMEDIAN NaN protected median value.
%   NANMEDIAN(X) returns the median treating NaNs as missing values.
%   For vectors, NANMEDIAN(X) is the median value of the non-NaN
%   elements in X.  For matrices, NANMEDIAN(X) is a row vector
%   containing the median value of each column, ignoring NaNs.
%
%   See also NANMEAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2013/02/07 09:13:00 $

[m,n] = size(x);
x = sort(x); % NaNs are forced to the bottom of each column

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));
if min(size(x))==1
  n = length(x)-sum(nans);
  if n == 0
    y = NaN;
  else
    if rem(n,2)     % n is odd
      y = x((n+1)/2);
    else            % n is even
      y = (x(n/2) + x(n/2+1))/2;
    end
  end
else
  n = size(x,1)-sum(nans);
  y = zeros(size(n));

  % Odd columns
  odd = find(rem(n,2)==1 & n>0);
  idx =(n(odd)+1)/2 + (odd-1)*m;
  y(odd) = x(idx);

  % Even columns
  even = find(rem(n,2)==0 & n>0);
  idx1 = n(even)/2 + (even-1)*m;
  idx2 = n(even)/2+1 + (even-1)*m;
  y(even) = (x(idx1)+x(idx2))/2;

  % All NaN columns
  i = find(n==0);
  y(i) = i + nan;
end
end

function y = nansum(x)
%NANSUM Sum ignoring NaNs.
%   NANSUM(X) returns the sum treating NaNs as missing values.
%   For vectors, NANSUM(X) is the sum of the non-NaN elements in
%   X. For matrices, NANSUM(X) is a row vector containing the sum
%   of the non-NaN elements in each column of X.
%
%    See also NANMEDIAN, NANSTD, NANMIN, NANMAX, NANMEAN.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 2.10 $  $Date: 2002/01/17 21:31:14 $

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

% Protect against an entire column of NaNs
y = sum(x);
i = find(all(nans));
y(i) = i + NaN;
end


function out = AverageData(data,x,y,avWindow,flagSerial)
% average data with time x to time y using window
dtx = median(diff(x));
if nargin<4, avWindow = median(diff(y)); end
if nargin<5, flagSerial = 0; end
dt2 = avWindow/2;
ndataOut = length(y);

% Pad data with NaNs from each side
nPointToAdd = ceil(dt2/dtx);
padNan = zeros(nPointToAdd,size(data,2))*NaN;
data = [padNan; data; padNan];
padTime = dtx*(1:nPointToAdd);
x = [x(1)-fliplr(padTime)'; x; x(end)+padTime'];

out = zeros(ndataOut,size(data,2));
if flagSerial % Serial execution
  for i=1:length(y)
    out(i,:) = FastNanMean(data,x>=y(i)-dt2 & x<y(i)+dt2);
  end
else % Parallel execution
  parfor i=1:length(y)
    out(i,:) = FastNanMean(data,x>=y(i)-dt2 & x<y(i)+dt2);
  end
end
end
function m = FastNanMean(x,idx)
% Faster version of nanmean()
xx = x(idx,:);
% Find NaNs and set them to zero
nans = isnan(xx); xx(nans) = 0;
% Count up non-NaNs.
n = sum(~nans,1);
n(n==0) = NaN; % prevent divideByZero warnings
% Sum up non-NaNs, and divide by the number of non-NaNs.
m = sum(xx,1) ./ n;
m(n<size(xx,1)*0.75) = NaN; % minDataFrac = .075
end
