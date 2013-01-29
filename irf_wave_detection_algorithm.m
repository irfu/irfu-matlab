function h=irf_wave_detection_algorithm(b,varargin)

% To duplicate the method of Bortnik et al. 2007

  %% Check input 
  [ax,args,nargs] = axescheck(varargin{:});

  %% get background magnetic field
  sampl_low=5;
  t_low=b(1,1):1/sampl_low:b(end,1); t_low=t_low';
  b_low=irf_resamp(b,t_low); %low sample frequency to avoid filter issues
  bf=irf_filt(b_low,1/600,0,[],5);
  b0=b_low;
  b0(:,2:4)=b_low(:,2:4)-bf(:,2:4);
  %b=b_low;
  %t=t_low;
  B=b0;
  sampl=sampl_low;
% %% resample to 25 Hz
  sampl_b=1/(b(2,1)-b(1,1));
  sampl=10;
  t=b(1,1):1/sampl:b(end,1); t=t'; 
  B=irf_resamp(B,t);
  b=irf_resamp(b,t); disp('resampling to 10 Hz');
  b(:,2:4)=b(:,2:4)-B(:,2:4);


  %% Remove the last sample if the total number of samples is odd

  if size(b,1)/2 ~= floor(size(b,1)/2)
    b=b(1:end-1,:);
    B=B(1:end-1,:);
  end
  
    % set to zero NaNs
  ind_nan_b=isnan(b); b(ind_nan_b)=0;
  ind_nan_B=isnan(B); B(ind_nan_B)=0;


  %% Find the frequencies for an FFT of all data

  nd2=size(b,1)/2;
  nyq=1/2;
  freq=sampl*(1:nd2)/(nd2)*nyq;
  w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

  %% Set some important parameters
  freq_int=[.01 5];
  freq_number=45;
  Morlet_width=5.36;

  if isnumeric(args{end})
      freq_int=args{end};
  elseif ismember({'freq'},args)
      disp('frequency interval values missing. using default')
  end
  
% Btot=irf_resamp(B(:,1),t);
% B2=irf_resamp(B,t);
% Btot(:,2)=sqrt(B2(:,2).*B2(:,2)+B2(:,3).*B2(:,3)+B2(:,4).*B2(:,4));
% %topCutoffFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi;
% %middleCutoffFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./4;
% %bottomCutoffFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./16;
% hCyclFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi;
% heCyclFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./4;
% oCyclFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./16;
% freq_int=[min(oCyclFreq)-.1*min(oCyclFreq) max(heCyclFreq)+.1*max(heCyclFreq)];
% freq_int=[.02 max(heCyclFreq)+.1*max(heCyclFreq)];


  amin=log10(0.5*sampl/freq_int(2));amax=log10(0.5*sampl/freq_int(1));anumber=freq_number;
  a=logspace(amin,amax,anumber);
  w0=sampl/2; % The maximum frequency
  sigma=Morlet_width/w0; % The width of the Morlet wavelet

  Swb=fft(b(:,2:4),[],1);
  
%   %frequency bin half width
%   freqBins=w0./logspace(amin,amax,2*freq_number-1);
%   fbhw=zeros(size(freq_number));
%   j=2;
%   for i=2:2:2*freq_number-3,
%       fbhw(j)=(abs(freqBins(i)-freqBins(i+1))+abs(freqBins(i+1)-freqBins(i+2)))/2;
%       j=j+1;
%   end
%   fbhw(1)=freqBins(1)-freqBins(3)-fbhw(2);
%   fbhw(freq_number)=freqBins(end-2)-freqBins(end)-fbhw(freq_number-1);
%   format shortE
%   display(fbhw')
%   
%% Get the correct frequencies for the wavelet transform
newfreq=w0./a;
%display(newfreq')
ndata = size(b,1); nfreq = length(a);
powerBx_SM_plot = zeros(ndata,nfreq);
powerBy_SM_plot = zeros(ndata,nfreq);
powerBz_SM_plot = zeros(ndata,nfreq);
powerCrossCov_SM_plot = zeros(ndata,nfreq);

parfor ind_a=1:length(a),
  mWexp = exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  mWexp = repmat(mWexp,1,3);
  Wwb = sqrt(1).*Swb.*mWexp;
  
  %% Get the wavelet transform by IFFT of the FFT
  Wb = ifft(Wwb,[],1);
  
  %% Calculate the power spectrum
  newfreqmat=w0/a(ind_a);
  powerB = 2*pi*(Wb.*conj(Wb))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
    
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

%   %% average spectral matrix over one wave period
% 
%   sampl_av = fix(sampl/a(ind_a));
%   if sampl_av/2 == floor(sampl_av/2)
%     sampl_av=sampl_av+1;
%   end
% 
%   for i = sampl_av:ndata-sampl_av,
%       if sampl_av < 2
%           SMpermute(i,:,:)=SMpermute(i,:,:);
%       elseif sampl_av > 2 && sampl_av < 4
%           SMpermute(i,:,:) = (SMpermute(i-1,:,:)+SMpermute(i,:,:)+...
%               SMpermute(i+1,:,:))./sampl_av;
%       else 
%           SMpermute(i,:,:) = sum(SMpermute(i-((sampl_av-1)/2):i+...
%               ((sampl_av-1)/2),:,:),1)./sampl_av;
%       end
%   end       

  %% Remove data possibly influenced by edge effects
  censur=floor(2*a);
  censur_indexes=[1:min(censur(ind_a),size(b,1)) max(1,size(b,1)-censur(ind_a)):size(b,1)];
  SMpermute(censur_indexes,:,:) = NaN;
    
  %% power
  powerBx_SM_plot(:,ind_a) = SMpermute(:,1,1);
  powerBy_SM_plot(:,ind_a) = SMpermute(:,2,2);
  powerBz_SM_plot(:,ind_a) = SMpermute(:,3,3)+SMpermute(:,1,1)+SMpermute(:,2,2);
  powerCrossCov_SM_plot(:,ind_a) = SMpermute(:,3,3)+SMpermute(:,1,1)+SMpermute(:,2,2);
  %powerCrossCov_SM_plot(:,ind_a) = real(SMpermute(:,1,2))+real(SMpermute(:,1,3))+real(SMpermute(:,2,3));
  %powerCrossCov_SM_plot(:,ind_a) = (SMpermute(:,1,2))+(SMpermute(:,1,3))+(SMpermute(:,2,3));
end

idx_nan_b = sum(ind_nan_b,2)>0;
powerCrossCov_SM_plot(idx_nan_b,:) = NaN;


%% average down the signal
ndata2=floor(ndata./300);
powerCrossCov_avg = zeros(ndata2,nfreq);
ind_nan_b_avg = zeros(size(powerCrossCov_avg,1));
for i=1:ndata2,
  powerCrossCov_avg(i,:) = mean(powerCrossCov_SM_plot(i.*300-299:i.*300,:),1);
end
%powerCrossCov_avg = powerCrossCov_SM_plot(1:ndata2,:);
%powerCrossCov_avg = powerCrossCov_SM_plot(30:30:end,:);
%powerCrossCov_SM_plot = flipdim(powerCrossCov_avg,2);
powerCrossCov_SM_plot = powerCrossCov_avg;
%powerCrossCov_SM_plot(40,:)
ndata=ndata2;
sampl=10/300;
t=b(1,1):1/sampl:b(end,1); t=t';

%% remove background median over a window
t_start_epoch=get_t_start_epoch(t(1,1));
toa = 4; %time of average for background median

if t(end)-t(1) <= toa*3600,
    medianPower = nanmedian(powerCrossCov_SM_plot);
    medianPower = repmat(medianPower,ndata,1);
elseif t(end)-t(1) <= (toa+.5*toa)*3600 & t(end)-t(1) > toa*3600,
    medianPower = zeros(size(powerCrossCov_SM_plot));
    medianPower(1:toa*3600*sampl,:) = repmat(nanmedian(...
        powerCrossCov_SM_plot(1:toa*3600*sampl,:)),toa*3600*sampl,1);
    medianPower(end-toa*3600*sampl:end,:) = repmat(nanmedian(...
        powerCrossCov_SM_plot(end-toa*3600*sampl:end,:)),...
        max(size(medianPower,1)-toa*3600*sampl-1:size(medianPower,1))...
        -min(size(medianPower,1)-toa*3600*sampl-1:size(medianPower,1)),1);
else
    medianPower = zeros(size(powerCrossCov_SM_plot));
    medianPower(1:.5*toa*3600*sampl,:) = repmat(nanmedian(...
        powerCrossCov_SM_plot(1:toa*3600*sampl,:)),.5*toa*3600*sampl,1);
%     medianPower(2*3600*sampl:end-2*3600*sampl,:) = repmat(nanmedian(...
%         powerCrossCov_SM_plot(2*3600*sampl:end-2*3600*sampl,:)),...
%         2*3600*sampl:size(powerCrossCov_SM_plot,1)-2*3600*sampl,1);
    medianPower(end-.5*toa*3600*sampl+1:end,:) = repmat(nanmedian(...
        powerCrossCov_SM_plot(end-toa*3600*sampl:end,:)),...
        max(size(powerCrossCov_SM_plot,1)-.5*toa*3600*sampl:size(powerCrossCov_SM_plot,1))...
        -min(size(powerCrossCov_SM_plot,1)-.5*toa*3600*sampl:size(powerCrossCov_SM_plot,1)),1);
    for i = .5*toa*3600*sampl+1:size(powerCrossCov_SM_plot,1)-.5*toa*3600*sampl,
        medianPower(i,:) = nanmedian(powerCrossCov_SM_plot(i-.5*toa*3600*sampl+1:i+.5*toa*3600*sampl,:));
    end
        
end
power_median_removed = log10(abs(powerCrossCov_SM_plot)) - log10(abs(medianPower));

%idx_nan_b = sum(ind_nan_b,2)>0;
%powerCrossCov_SM_plot(idx_nan_b,:) = NaN;
%power_median_removed(idx_nan_b,:) = NaN;
powerBx_SM_plot(idx_nan_b,:) = NaN;
powerBz_SM_plot(idx_nan_b,:) = NaN;
peakInd = power_median_removed < ones(size(power_median_removed));
%power_median_removed(peakInd) = NaN;

% waveFrequencies=zeros(ndata,6);
% for i=1:ndata,
%     peakInd = power_median_removed(i,:)>1;
%     peakFreq = newfreq(peakInd);
%     spectralPeak = max(power_median_removed(i,peakInd));
%     spectralPeakInd = power_median_removed(i,peakInd) == max(power_median_removed(i,peakInd));
%     if numel(peakFreq) > 0,
%       bottomFreq = peakFreq(1);
%       maxPowerFreq = peakFreq(spectralPeakInd);
%       topFreq = peakFreq(end);
%       doublePeak = power_median_removed(i,min(find(peakInd)):max(find(peakInd)));
%       if any(doublePeak < 1),
%           topFreq2 = peakFreq(end);
%           lowPoints = find(doublePeak < 1);
%           topFreq = newfreq(lowPoints(1));
%           bottomFreq2 = newfreq(lowPoints(end));
%           maxPowerFreq2 = NaN;
%       else
%           bottomFreq2 = NaN;
%           maxPowerFreq2 = NaN;
%           topFreq2 = NaN;
%       end    
%     else
%       bottomFreq = NaN;
%       maxPowerFreq = NaN;
%       topFreq = NaN;
%     end
%     waveFrequencies(i,1) = bottomFreq;
%     waveFrequencies(i,2) = maxPowerFreq;
%     waveFrequencies(i,3) = topFreq;
%     waveFrequencies(i,4) = bottomFreq2;
%     waveFrequencies(i,5) = maxPowerFreq2;
%     waveFrequencies(i,6) = topFreq2;
% end

waveFrequencies=zeros(ndata,6);
threshold=1;
for i=1:ndata,
    peakInd = power_median_removed(i,:)>threshold;
    findIndices = find(peakInd);
    peakFreq = newfreq(peakInd);
    spectralPeak = max(power_median_removed(i,peakInd));
    spectralPeakInd = power_median_removed(i,peakInd) == max(power_median_removed(i,peakInd));
    if numel(peakFreq) > 0,
      topFreq = peakFreq(1);
      maxPowerFreq = peakFreq(spectralPeakInd);
      bottomFreq = peakFreq(end);
      if any(power_median_removed(i,findIndices(1):findIndices(end)) < threshold),
          topFreq2 = topFreq;
          bottomFreq2 = newfreq(min(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))-1+findIndices(1));
          topFreq = newfreq(max(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))+findIndices(1));
%           if maxPowerFreq > bottomFreq2,
%               spectralPeak2 = max(power_median_removed(i,findIndices(1):min(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))-1+findIndices(1)));
%               %%%spectralPeak2 = max(power_median_removed(i,bottomFreq:topFreq));
%               maxPowerFreq2 = maxPowerFreq;
%               spectralPeakInd2 = power_median_removed(i,findIndices(1):min(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))-1+findIndices(1)) == spectralPeak2;
%               %%%spectralPeakInd2 = power_median_removed(i,bottomFreq:topFreq) == spectralPeak2;
%               maxPowerFreq = peakFreq(spectralPeakInd2);
%               %maxPowerFreq2 = peakFreq(spectralPeakInd2);
%           else
%               spectralPeak3 = max(power_median_removed(i,max(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))+findIndices(1)):findIndices(end));
%               %spectralPeakInd3 = find(max(power_median_removed(i,max(find(power_median_removed(i,findIndices(1):findIndices(end)) < 1))+findIndices(1)):findIndices(end)));
%               spectralPeakInd3 = (power_median_removed(i,max(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))+findIndices(1)):findIndices(end)) == spectralPeak3;
%               %maxPowerFreq2 = peakFreq(spectralPeakInd3);
%               %maxPowerFreq2 = maxPowerFreq;
%               maxPowerFreq2 = newfreq(spectralPeakInd3);
%               %maxPowerFreq = newfreq(spectralPeakInd3);
%               %spectralPeak3 = power_median_removed(i,findIndices(1):min(find(power_median_removed(i,findIndices(1):findIndices(end)) < threshold))+findIndices(1));
%               %spectralPeakInd3 = max(spectralPeak3) == spectralPeak3;
%               %maxPowerFreq2 = newfreq(spectralPeakInd3);
%               
%           end
          maxPowerFreq2 = NaN;
      else
          bottomFreq2 = NaN;
          maxPowerFreq2 = NaN;
          topFreq2 = NaN;
      end    
         
    else
      bottomFreq = NaN;
      maxPowerFreq = NaN;
      topFreq = NaN;
      bottomFreq2 = NaN;
      maxPowerFreq2 = NaN;
      topFreq2 = NaN;
    end
    waveFrequencies(i,1) = bottomFreq;
    waveFrequencies(i,2) = maxPowerFreq;
    waveFrequencies(i,3) = topFreq;
    waveFrequencies(i,4) = bottomFreq2;
    waveFrequencies(i,5) = maxPowerFreq2;
    waveFrequencies(i,6) = topFreq2;
end

%% set cutoff frequencies and put a limit on the width of the emic event
Btot=irf_resamp(B(:,1),t);
B2=irf_resamp(B,t);
Btot(:,2)=sqrt(B2(:,2).*B2(:,2)+B2(:,3).*B2(:,3)+B2(:,4).*B2(:,4));
%topCutoffFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi;
%middleCutoffFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./4;
%bottomCutoffFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./16;
hCyclFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi;
heCyclFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./4;
oCyclFreq = Btot(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./16;
bottomCutoffFreq = heCyclFreq;   
topCutoffFreq = hCyclFreq;  %or should I use 10?
deltaPeak = .03;
 validInd = waveFrequencies(:,3) > topCutoffFreq;
 validInd2 = waveFrequencies(:,1) < bottomCutoffFreq;
%validInd = waveFrequencies(:,3) > newfreq(2);
%validInd2 = waveFrequencies(:,1) < newfreq(end-1);
 validInd3 = abs(waveFrequencies(:,1) - waveFrequencies(:,3)) < deltaPeak;
% waveFrequencies(validInd3,1:3) = NaN;
waveFrequencies(validInd2,1:3) = NaN;
waveFrequencies(validInd,1:3) = NaN;
 validInd4 = waveFrequencies(:,6) > topCutoffFreq;
 validInd5 = waveFrequencies(:,4) < bottomCutoffFreq;
%validInd4 = waveFrequencies(:,6) > newfreq(2);
%validInd5 = waveFrequencies(:,4) < newfreq(end-1);
 validInd6 = abs(waveFrequencies(:,4) - waveFrequencies(:,6)) < deltaPeak;
waveFrequencies(validInd4,4:6) = NaN;
waveFrequencies(validInd5,4:6) = NaN;
% waveFrequencies(validInd6,4:6) = NaN;
 validInd7 = waveFrequencies(:,3)./waveFrequencies(:,1) > 4;
 waveFrequencies(validInd7,1:3) = NaN;
 validInd8 = waveFrequencies(:,6)./waveFrequencies(:,4) > 4;
 waveFrequencies(validInd8,4:6) = NaN;

nanmean(abs(waveFrequencies(:,1) - waveFrequencies(:,3)))
nanmean(waveFrequencies(:,2))

% %discard impulsive bursts
% for i=4000:ndata-4000,
%     if waveFrequencies(i,1) ~= NaN,
%         samplesInWindow = length(find(waveFrequencies(i-3999:i+3999,1) ~= NaN));
%         if samplesInWindow < 4000,
%             waveFrequencies(i-3999:i+3999,:) = NaN;
%         end
%     end
% end

%check for association from one time step to the next
TT = irf.TimeTable
% i=1;
% while i < ndata,
%     if waveFrequencies(i,1) ~= NaN, %| waveFrequencies(i,4) ~= NaN,
%         eventStartTime = t(i);
%         j=1;
%         while (i+j+4<ndata & waveFrequencies(i+j,1) < waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j,3) > waveFrequencies(i+j-1,1))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+1,1) < waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j+1,3) > waveFrequencies(i+j-1,1))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+2,1) < waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j+2,3) > waveFrequencies(i+j-1,1))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+3,1) < waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j+3,3) > waveFrequencies(i+j-1,1))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+4,1) < waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j+4,3) > waveFrequencies(i+j-1,1)), % | ...
% %                 ...
% %                 (i+j+4<ndata & waveFrequencies(i+j,1) < waveFrequencies(i+j-1,6)...
% %                 & waveFrequencies(i+j,3) > waveFrequencies(i+j-1,4))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+1,1) < waveFrequencies(i+j-1,6)...
% %                 & waveFrequencies(i+j+1,3) > waveFrequencies(i+j-1,4))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+2,1) < waveFrequencies(i+j-1,6)...
% %                 & waveFrequencies(i+j+2,3) > waveFrequencies(i+j-1,4))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+3,1) < waveFrequencies(i+j-1,6)...
% %                 & waveFrequencies(i+j+3,3) > waveFrequencies(i+j-1,4))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+4,1) < waveFrequencies(i+j-1,6)...
% %                 & waveFrequencies(i+j+4,3) > waveFrequencies(i+j-1,4)),
%             if (waveFrequencies(i+j,1) < waveFrequencies(i+j-1,3)) | ...
%                     (waveFrequencies(i+j,1) < waveFrequencies(i+j-1,6)),
%                 j=j+1;
%             elseif (waveFrequencies(i+j+1,1) < waveFrequencies(i+j-1,3)) | ...
%                     (waveFrequencies(i+j+1,1) < waveFrequencies(i+j-1,6))
%                 j=j+2;
%             elseif (waveFrequencies(i+j+2,1) < waveFrequencies(i+j-1,3)) | ...
%                     (waveFrequencies(i+j+2,1) < waveFrequencies(i+j-1,6))
%                 j=j+3;
%             elseif (waveFrequencies(i+j+3,1) < waveFrequencies(i+j-1,3)) | ...
%                     (waveFrequencies(i+j+3,1) < waveFrequencies(i+j-1,6))
%                 j=j+4;
%             else
%                 j=j+5;
%             end
%         end
% %         while (i+j<ndata & waveFrequencies(i+j,1) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j,3) < waveFrequencies(i+j-1,1)),
% %             j=j+1;
% %         end
% % 
%         eventEndTime = t(i+j);
%         if j > 3,  %NOTE: CHANGE THIS IN TWO PLACES
%             tint = [eventStartTime eventEndTime];
%             TT=add(TT,tint,{'EMIC events','for Maarble'},'test');
%         end
%         i=i+j;
%     end
%     i=i+1;
% end
% i=1
% while i < ndata,
%     if waveFrequencies(i,4) ~= NaN, %| waveFrequencies(i,4) ~= NaN,
%         eventStartTime = t(i);
%         j=1;
%         while (i+j+4<ndata & waveFrequencies(i+j,4) > waveFrequencies(i+j-1,6)...
%                 & waveFrequencies(i+j,6) < waveFrequencies(i+j-1,4))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+1,4) > waveFrequencies(i+j-1,6)...
%                 & waveFrequencies(i+j+1,6) < waveFrequencies(i+j-1,4))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+2,4) > waveFrequencies(i+j-1,6)...
%                 & waveFrequencies(i+j+2,6) < waveFrequencies(i+j-1,4))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+3,4) > waveFrequencies(i+j-1,6)...
%                 & waveFrequencies(i+j+3,6) < waveFrequencies(i+j-1,4))...
%                 | (i+j+4<ndata & waveFrequencies(i+j+4,4) > waveFrequencies(i+j-1,6)...
%                 & waveFrequencies(i+j+4,6) < waveFrequencies(i+j-1,4)),
% %                 ...
% %                 | (i+j+4<ndata & waveFrequencies(i+j,4) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j,6) < waveFrequencies(i+j-1,1))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+1,4) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j+1,6) < waveFrequencies(i+j-1,1))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+2,4) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j+2,6) < waveFrequencies(i+j-1,1))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+3,4) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j+3,6) < waveFrequencies(i+j-1,1))...
% %                 | (i+j+4<ndata & waveFrequencies(i+j+4,4) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j+4,6) < waveFrequencies(i+j-1,1)),
%             if (waveFrequencies(i+j,4) > waveFrequencies(i+j-1,6)) | ...
%                     (waveFrequencies(i+j,4) > waveFrequencies(i+j-1,3)),
%                 j=j+1;
%             elseif (waveFrequencies(i+j+1,4) > waveFrequencies(i+j-1,6)) | ...
%                     (waveFrequencies(i+j+1,4) > waveFrequencies(i+j-1,3))
%                 j=j+2;
%             elseif (waveFrequencies(i+j+2,4) > waveFrequencies(i+j-1,6)) | ...
%                     (waveFrequencies(i+j+2,4) > waveFrequencies(i+j-1,3))
%                 j=j+3;
%             elseif (waveFrequencies(i+j+3,4) > waveFrequencies(i+j-1,6)) | ...
%                     (waveFrequencies(i+j+3,4) > waveFrequencies(i+j-1,3))
%                 j=j+4;
%             else
%                 j=j+5;
%             end
%         end
% %         while (i+j<ndata & waveFrequencies(i+j,1) > waveFrequencies(i+j-1,3)...
% %                 & waveFrequencies(i+j,3) < waveFrequencies(i+j-1,1)),
% %             j=j+1;
% %         end
% % 
%         eventEndTime = t(i+j);
%         if j > 3, %NOTE: CHANGE THIS IN TWO PLACES
%             tint = [eventStartTime eventEndTime];
%             TT=add(TT,tint,{'EMIC events','for Maarble'},'test');
%         end
%         i=i+j;
%     end
%     i=i+1;
% end
   
i=1;
while i < ndata,
    if waveFrequencies(i,1) ~= NaN, %| waveFrequencies(i,4) ~= NaN,
        eventStartTime = t(i);
        j=1;
        while (i+j+4<ndata & waveFrequencies(i+j,1) < waveFrequencies(i,3)...
                & waveFrequencies(i+j,3) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+1,1) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+1,3) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+2,1) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+2,3) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+3,1) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+3,3) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+4,1) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+4,3) > waveFrequencies(i,1))  | ...
                ...
                (i+j+4<ndata & waveFrequencies(i+j,4) < waveFrequencies(i,3)...
                & waveFrequencies(i+j,6) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+1,4) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+1,6) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+2,4) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+2,6) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+3,4) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+3,6) > waveFrequencies(i,1))...
                | (i+j+4<ndata & waveFrequencies(i+j+4,4) < waveFrequencies(i,3)...
                & waveFrequencies(i+j+4,6) > waveFrequencies(i,1)),
            if (waveFrequencies(i+j,1) < waveFrequencies(i,3)) | ...
                    (waveFrequencies(i+j,4) < waveFrequencies(i,3)),
                j=j+1;
            elseif (waveFrequencies(i+j+1,1) < waveFrequencies(i,3)) | ...
                    (waveFrequencies(i+j+1,4) < waveFrequencies(i,3))
                j=j+2;
            elseif (waveFrequencies(i+j+2,1) < waveFrequencies(i,3)) | ...
                    (waveFrequencies(i+j+2,4) < waveFrequencies(i,3))
                j=j+3;
            elseif (waveFrequencies(i+j+3,1) < waveFrequencies(i,3)) | ...
                    (waveFrequencies(i+j+3,4) < waveFrequencies(i,3))
                j=j+4;
            else
                j=j+5;
            end
        end
%         while (i+j<ndata & waveFrequencies(i+j,1) > waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j,3) < waveFrequencies(i+j-1,1)),
%             j=j+1;
%         end
% 
        eventEndTime = t(i+j);
        if j > 4,  %NOTE: CHANGE THIS IN TWO PLACES
            tint = [eventStartTime eventEndTime];
            TT=add(TT,tint,{'EMIC events','for Maarble'},'test');
        end
        i=i+j;
    end
    i=i+1;
end
i=1;
while i < ndata,
    if waveFrequencies(i,4) ~= NaN, %| waveFrequencies(i,4) ~= NaN,
        eventStartTime = t(i);
        j=1;
        while (i+j+4<ndata & waveFrequencies(i+j,4) > waveFrequencies(i,6)...
                & waveFrequencies(i+j,6) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+1,4) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+1,6) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+2,4) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+2,6) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+3,4) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+3,6) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+4,4) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+4,6) < waveFrequencies(i,4)) ...
                ...
                | (i+j+4<ndata & waveFrequencies(i+j,1) > waveFrequencies(i,6)...
                & waveFrequencies(i+j,3) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+1,1) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+1,3) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+2,1) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+2,3) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+3,1) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+3,3) < waveFrequencies(i,4))...
                | (i+j+4<ndata & waveFrequencies(i+j+4,1) > waveFrequencies(i,6)...
                & waveFrequencies(i+j+4,3) < waveFrequencies(i,4)),
            if (waveFrequencies(i+j,4) > waveFrequencies(i+j-1,6)) | ...
                    (waveFrequencies(i+j,1) > waveFrequencies(i,6)),
                j=j+1;
            elseif (waveFrequencies(i+j+1,4) > waveFrequencies(i,6)) | ...
                    (waveFrequencies(i+j+1,1) > waveFrequencies(i,6))
                j=j+2;
            elseif (waveFrequencies(i+j+2,4) > waveFrequencies(i,6)) | ...
                    (waveFrequencies(i+j+2,1) > waveFrequencies(i,6))
                j=j+3;
            elseif (waveFrequencies(i+j+3,4) > waveFrequencies(i,6)) | ...
                    (waveFrequencies(i+j+3,1) > waveFrequencies(i,6))
                j=j+4;
            else
                j=j+5;
            end
        end
%         while (i+j<ndata & waveFrequencies(i+j,1) > waveFrequencies(i+j-1,3)...
%                 & waveFrequencies(i+j,3) < waveFrequencies(i+j-1,1)),
%             j=j+1;
%         end
% 
        eventEndTime = t(i+j);
        if j > 4, %NOTE: CHANGE THIS IN TWO PLACES
            tint = [eventStartTime eventEndTime];
            TT=add(TT,tint,{'EMIC events','for Maarble'},'test');
        end
        i=i+j;
    end
    i=i+1;
end

numel(TT)
ascii(TT)
%TT.TimeInterval(4,:)
export_ascii(TT,'/Users/meghanmella/Documents/MATLAB/emic.txt')
             
% length(t)
% size(powerCrossCov_SM_plot)
% length(newfreq)
  
%% plot

figure(318), clf 
clear h;
npl=2;clf;
ipl=1;
   it2 = 0:.018:.9; it2=it2';
   it3 = ones(size(it2)); 
   it3(end)=.9;
   it=0:.02:1;it=it'; 
   xcm=[ [0*it flipud(it) it];[it2 it2 it3];...
       [flipud(it3) flipud(it2) flipud(it2)]; [flipud(it) 0*it 0*it]];
   xcm2=[ [it2 it2 it3];[flipud(it3) flipud(it2) flipud(it2)]];
   clear it; clear it2; clear it3;

   %t=b(:,1);
   %t_start_epoch=get_t_start_epoch(t(1,1));
   
       %%%%% B spectra from spectral matrix
   
    h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
    %pcolor(t-t_start_epoch,newfreq,log10(abs(powerCrossCov_SM_plot.'))) 
    pcolor(t-t_start_epoch,newfreq,log10(abs(powerCrossCov_SM_plot.')))
    %pcolor(log10(abs(powerCrossCov_SM_plot.'))) 
    shading flat
    ylabel('f [Hz]')
    set(gca,'yscale','log','tickdir','out');
  %  set(gca,'tickdir','out');
    cmean=nanmean(nanmean(log10(abs(powerCrossCov_SM_plot))));
    axis(axis)
    hold on
    %plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1.6e-19./1.67e-27);
    plot(t-t_start_epoch,hCyclFreq,'color','k');
    plot(t-t_start_epoch,heCyclFreq,'--','color','k');
    plot(t-t_start_epoch,oCyclFreq,'-.','color','k');
    %axis([min(t-t_start_epoch) max(t-t_start_epoch) min(newfreq) max(newfreq)]);
    hold off
    caxis(floor(cmean)+[-3.5 3.5]);
    %caxis([-3.5 3.5]);
    %caxis([-6.5 0.5]);
    colormap(xcm);
    hca2=colorbar('peer',h(1));
    ylabel(hca2,{'log(B)'; '[nT^2/Hz]'});

           %%%%% B spectra from spectral matrix
   
    h(ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
    %pcolor(t-t_start_epoch,newfreq,log10(abs(powerCrossCov_SM_plot.'))) 
    pcolor(t-t_start_epoch,newfreq,log10(abs(power_median_removed.'))) 
    %pcolor(log10(abs(power_median_removed.'))) 
    %pcolor(t-t_start_epoch,newfreq,log10(abs(medianPower.'))) 
    shading flat
    ylabel('f [Hz]')
    set(gca,'yscale','log','tickdir','out');
  %  set(gca,'tickdir','out');
    cmean=nanmean(nanmean(log10(abs(power_median_removed))));
%     hold on
%     %plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1.6e-19./1.67e-27);
%     plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi,'color','k');
%     plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./4,'--','color','k');
%     plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1e-9.*1.6e-19./1.67e-27./2./pi./16,'-.','color','k');
%     axis([min(t-t_start_epoch) max(t-t_start_epoch) min(newfreq) max(newfreq)]);
%     hold off
    %caxis(floor(cmean)+[-3.5 3.5]);
    axis(axis)
    hold on
    plot(t-t_start_epoch,waveFrequencies(:,1),'o','markerfacecolor','k','markeredgecolor','k','markersize',4)    
    plot(t-t_start_epoch,waveFrequencies(:,2),'d','markerfacecolor','k','markeredgecolor','k','markersize',6)    
    plot(t-t_start_epoch,waveFrequencies(:,3),'o','markerfacecolor','k','markeredgecolor','k','markersize',4)
    plot(t-t_start_epoch,waveFrequencies(:,4),'o','markerfacecolor','k','markeredgecolor','k','markersize',4)    
    plot(t-t_start_epoch,waveFrequencies(:,5),'d','markerfacecolor','k','markeredgecolor','k','markersize',6)    
    plot(t-t_start_epoch,waveFrequencies(:,6),'o','markerfacecolor','k','markeredgecolor','k','markersize',4)

    %plot(BVector(:,1)-t_start_epoch,BVector(:,2).*1.6e-19./1.67e-27);
    plot(t-t_start_epoch,hCyclFreq,'color','k');
    plot(t-t_start_epoch,heCyclFreq,'--','color','k');
    plot(t-t_start_epoch,oCyclFreq,'-.','color','k');
    %axis([min(t-t_start_epoch) max(t-t_start_epoch) min(newfreq) max(newfreq)]);

    hold off
    caxis([-2 2]);
    colormap(xcm);
    hca2=colorbar('peer',h(2));
    ylabel(hca2,{'log(B)'; '[nT^2/Hz]'});

   
    %% Add figure menu
    irf_figmenu;


      axes(h(1));
      %title(['Width Morlet wavelet = ' num2str(Morlet_width)]);
      %ht=irf_pl_info([mfilename '  ' datestr(now)]); set(ht,'interpreter','none'); % add information to the plot
      irf_zoom(h,'x',[min(t) max(t)]);
      %irf_legend(h(1),'Figure reference',[0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);
      %following three lines from irf_timeaxis, because the date kept
      %appearing under the position labels
      xlimlast=get(h(end),'xlim');
      start_time = irf_time(xlimlast(1) + t_start_epoch,'vector');
      time_label = datestr( datenum(start_time),1 );
      
      irf_legend(h(1),['C1           ' time_label],[0 1.05],'fontsize',10,'color','cluster');
      irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
    %irf_legend('C4', 'color','cluster');
    %irf_figmenu;
    
%         figure(320), clf
%         %%%%% Frequency spectrum at a given time
% 
%     %plot(log10(abs(power_median_removed(10330,:))),'--')
%     plot(power_median_removed(110106,:),'--')
%     hold on
%     plot(ones(size(power_median_removed(110106,:))))
%     plot(log10(abs(medianPower(110106,:))),'-.')
%     plot(log10(abs(powerCrossCov_SM_plot(110106,:))))
%     hold off
%     ylabel('f [Hz]')
%   %  set(gca,'yscale','log','tickdir','out');
%     set(gca,'tickdir','out');
%     
%     
%     figure(321), clf
%     plot(waveFrequencies(:,1),'.','markerfacecolor','k','markersize',8)    
%     ylabel('f [Hz]')
%   %  set(gca,'yscale','log','tickdir','out');
%     set(gca,'tickdir','out');



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_start_epoch=get_t_start_epoch(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives back the value of t_start_epoch of the figure
% if not  set, sets t_start_epoch of the figure
ud=get(gcf,'userdata');
ii = find(~isnan(t));
if ii,
  valid_time_stamp=t(ii(1));
else
  valid_time_stamp=[];
end

if isfield(ud,'t_start_epoch'),
  t_start_epoch=ud.t_start_epoch;
elseif valid_time_stamp,
  if valid_time_stamp > 1e8, % set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
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


function y = nanmedian(x)
%NANMEDIAN NaN protected median value.
%   NANMEDIAN(X) returns the median treating NaNs as missing values.
%   For vectors, NANMEDIAN(X) is the median value of the non-NaN
%   elements in X.  For matrices, NANMEDIAN(X) is a row vector
%   containing the median value of each column, ignoring NaNs.
%
%   See also NANMEAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision$  $Date$

[m,n] = size(x);
x = sort(x); % NaNs are forced to the bottom of each column

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));
if min(size(x))==1,
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