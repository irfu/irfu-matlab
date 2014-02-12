function [b]=c_fgm_staff_combine(bFgm,bStaff,varargin)
% C_FGM_STAFF_COMBINE combine FGM and STAFF data series
% b=C_FGM_STAFF_COMBINE(bFgm,bStaff) combines FGM and STAFF time series
% for the case when both are in burst or normal mode. For filtering are
% used FIR filters with 127 points for the normal mode and 1023 points for
% the burst mode. The returned time series have a timeline of bStaff.
%
% C_FGM_STAFF_COMBINE(bFgm,bStaff,'plot') plot the spectra of FGM, STAFF
% and merged time series. This is useful to see how merging has worked and
% whether merging can be applied on the two time series. In cases when
% STAFF spectra have lower amplitude than FGM spectra just below the cut
% frequency of 1.3Hz, there maybe problems in the merged time series.
%

%% check input
if nargout==0 && nargin == 0, help c_fgm_staff_combine;return;end

%% define common time interval in which to combine data sets
tint=[max(bFgm(1,1),bStaff(1,1)) min(bFgm(end,1),bStaff(end,1))]; % define common time interval
bFgm=irf_tlim(bFgm,tint); % limit time series to common time interval
bStaff=irf_tlim(bStaff,tint); % limit time series to common time interval

%% check if FGM and STAFF are in normal mode
if 1/mean(diff(bStaff(1:min(end,100),1))) < 400 % STAFF in normal mode, not implemented
    fsStaff = 25; % Sampling Frequency
    fCut = 1.3; % cut Frequency (Hz)
    N = 127;  % FIR filter order (have to be odd number!!!)
    fprintf('\n STAFF data in normal mode [%d Hz]! Cut-off at %f Hz\n',fsStaff, fCut);
    %         fprintf('Adjust cut-off frequency if coincides with a kink in the power spectrum');
else
    fsStaff = 450; % Sampling Frequency
    fCut = 1.3; % cut Frequency (Hz)
    N = 1023;  % FIR filter order (have to be odd number!!!)
    fprintf('\n STAFF data in burst mode [%d Hz]! Cut-off at %f Hz\n', fsStaff, fCut);
    
end

if 1/mean(diff(bFgm(1:min(end,100),1))) < 60 % FGM in normal mode
    disp('FGM data in normal mode!');
    %Fs_fgm = 22; % Sampling Frequency
else
    if fsStaff == 25
        disp('FGM sample rate is higher than STAFF''!');
        b=[];
        return;
    end
    disp('FGM data in burst mode!');
    %Fs_fgm = 66.7; % Sampling Frequency
end
fprintf('\n Adjust cut-off frequency if needed\n')
%% design filters
% low pass filter for FGM
% higpass filter for STAFF
Sband = 60; %80; % approx. stop band attenuation in dB.

LoP = fir1(N-1,fCut/(fsStaff/2),'low',kaiser(N,Sband/10));

%resample FGM to match STAFF
b_fgm_res = irf_resamp(bFgm, bStaff);

b_fgm_filt=b_fgm_res;
for j=2:4,
    b_fgm_filt(:,j) = filtfilt(LoP, 1.0, b_fgm_res(:,j));
end
b_staff_filt=bStaff;
for j=2:4,
    b_staff_filt(:,j) = bStaff(:,j) - filtfilt(LoP, 1.0, bStaff(:,j));
end

%merge the two vectors
b=[b_fgm_res(:,1) b_fgm_filt(:,2:end) + b_staff_filt(:,2:end)];

if size(b,2)>4, % the 5th column should be |B|
    b(:,5)=irf_abs(b,1);
end


%% plot PSD

if nargin > 2 && ischar(varargin{1}) && strcmpi(varargin{1},'plot'),
    
    wd=2^14; % window size
    
    lw=1;
    h=spectrum.welch('Hann',wd,75);
    
    
    figure(13)
    for k=2:4;
        
        hpl(k-1)=subplot(3,1,k-1);
        psdestFGM=psd(h,b_fgm_res(:,k),'Fs',fsStaff);
        plot(psdestFGM.Frequencies, psdestFGM.Data);
        set(gca, 'XScale', 'log','YScale', 'log');
        hline = findobj(gca,'Type','line');
        set(hline,'color','r');
        ylabel('nT^2/Hz');   xlabel('Hz');
        
        hold on
        psdestSTA=psd(h,bStaff(:,k),'Fs',fsStaff);
        plot(psdestSTA.Frequencies, psdestSTA.Data);
        hline = findobj(gca,'Type','line','color','b');
        set(hline,'color','g');
        psdestMER=psd(h,b(:,k),'Fs',fsStaff);
        plot(psdestMER.Frequencies, psdestMER.Data);
        hline = findobj(gca,'Type','line','color','b');
        set(hline,'linestyle',':');
        legend(hpl(k-1),'FGM','STAFF','MERGED')
        hold off
        
    end
end


%% check the return
if nargout == 0, b=[]; end




