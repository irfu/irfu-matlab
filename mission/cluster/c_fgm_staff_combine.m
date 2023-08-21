function [bout]=c_fgm_staff_combine(bFgm,bStaff,varargin)
% C_FGM_STAFF_COMBINE  combine FGM and STAFF data series
%
% b = C_FGM_STAFF_COMBINE(bFgm,bStaff,[OPTIONS])
%    Combines FGM and STAFF time series. Combination is done using
%    FIR filters with 127 points (normal mode) or 1023 points (burst mode).
%    The returned time series have a timeline of bStaff.
%
% b = C_FGM_STAFF_COMBINE(bFgm,bStaff,'fCut',FCUT)
%    Set merging frequency to FCUT (default 1.3 Hz)
%
% b = C_FGM_STAFF_COMBINE(bFgm,bStaff,'cl_id',CL_ID)
%    Specify Cluster CS number (used for sampling frequency calculation)
%
% b = C_FGM_STAFF_COMBINE(bFgm,bStaff,'plot')
%    Plot the spectra of FGM, STAFF and merged time series.
%    This is useful to verify the merging. In cases the STAFF spectra have
%    lower amplitude than the FGM spectra at the FCUT frequency of 1.3Hz,
%    one should use a lower FCUT where the two spcetra agree.

%% check input
if nargout==0 && nargin == 0, help c_fgm_staff_combine;return;end

args = varargin;
if ~isempty(args), have_options = 1; else, have_options = 0; end

%% Default values that can be override by options
fCut = 1.3; % cut Frequency (Hz)
flag_plot = false;
cl_id = [];
while have_options
  l = 1;
  switch(lower(args{1}))
    case 'plot'
      flag_plot = true;
    case 'fcut'
      if length(args)>1 && isnumeric(args{2}) && args{2}>0
        fCut = args{2};
        l = 2;
      else % TODO implement string tint
        errS = 'fCut must be a positive number';
        irf.log('critical',errS), error(errS)
      end
    case 'cl_id'
      if length(args)>1 && isnumeric(args{2}) && ismember(args{2},1:4)
        cl_id = args{2};
        l = 2;
      else % TODO implement string tint
        errS = 'fCut must be a positive number';
        irf.log('critical',errS), error(errS)
      end
    otherwise
      irf.log('notice','unrecognized input')
  end
  args = args(l+1:end);
  if isempty(args), break, end
end

%% define common time interval in which to combine data sets
tint=[max(bFgm(1,1),bStaff(1,1)) min(bFgm(end,1),bStaff(end,1))]; % define common time interval
bFgm=irf_tlim(bFgm,tint); % limit time series to common time interval
bStaff=irf_tlim(bStaff,tint); % limit time series to common time interval

%% check if FGM and STAFF are in normal mode
if isempty(cl_id), fsStaff = c_efw_fsample(bStaff,'hx');
else, fsStaff = c_efw_fsample(bStaff,'hx',cl_id);
end
if fsStaff < 400 % STAFF in normal mode, not implemented
  N = 127;  % FIR filter order (have to be odd number!!!)
  irf.log('notice',sprintf('STAFF data in normal mode [%.2f Hz]',fsStaff));
else
  N = 1023;  % FIR filter order (have to be odd number!!!)
  irf.log('notice',sprintf('STAFF data in burst mode [%.2f Hz]',fsStaff));
end

if 1/mean(diff(bFgm(1:min(end,100),1))) < 60 % FGM in normal mode
  disp('FGM data in normal mode!');
  %Fs_fgm = 22; % Sampling Frequency
else
  if fsStaff == 25
    errS = 'FGM sample rate is higher than STAFF''!';
    irf.log('critical',errS), error(errS)
  end
  irf.log('notice','FGM data in burst mode!');
end
irf.log('warning',sprintf('fCut = %.2f Hz, adjust cut-off frequency if needed',fCut))
%% design filters
% low pass filter for FGM
% higpass filter for STAFF
Sband = 60; %80; % approx. stop band attenuation in dB.

LoP = fir1(N-1,fCut/(fsStaff/2),'low',kaiser(N,Sband/10));

%resample FGM to match STAFF
b_fgm_res = irf_resamp(bFgm, bStaff);

b_fgm_filt=b_fgm_res;
for j=2:4
  b_fgm_filt(:,j) = filtfilt(LoP, 1.0, b_fgm_res(:,j));
end
b_staff_filt=bStaff;
for j=2:4
  b_staff_filt(:,j) = bStaff(:,j) - filtfilt(LoP, 1.0, bStaff(:,j));
end

%merge the two vectors
b=[b_fgm_res(:,1) b_fgm_filt(:,2:end) + b_staff_filt(:,2:end)];

if size(b,2)>4 % the 5th column should be |B|
  b(:,5)=irf_abs(b,1);
end


%% plot PSD
if flag_plot
  wd=2^14; % window size
  h=spectrum.welch('Hann',wd,75);

  compS = 'xyz';
  figure(13), hpl = gobjects(3);
  for k=2:4
    hpl(k-1)=subplot(3,1,k-1);
    psdestFGM = psd(h,b_fgm_res(:,k),'Fs',fsStaff);
    plot(psdestFGM.Frequencies, psdestFGM.Data,'k');
    set(gca, 'XScale', 'log','YScale', 'log');
    ylabel(['PSD(B' compS(k-1) ') [nT^2/Hz]']); xlabel('Frequency [Hz]');
    if k==2, title(sprintf('fCut = %.2f Hz',fCut)), end
    hold on
    psdestSTA = psd(h,bStaff(:,k),'Fs',fsStaff);
    plot(psdestSTA.Frequencies, psdestSTA.Data,'g');
    psdestMER = psd(h,b(:,k),'Fs',fsStaff);
    plot(psdestMER.Frequencies, psdestMER.Data,'r:');
    legend(hpl(k-1),'FGM','STAFF','MERGED')
    ylim = get(gca,'YLim');
    plot([fCut fCut],ylim,'b-.')
    hold off
  end
end


%% check the return
if nargout ~= 0, bout=b; end




