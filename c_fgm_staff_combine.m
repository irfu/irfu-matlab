function [b]=c_fgm_staff_combine(b_fgm,b_staff)
% C_FGM_STAFF_COMBINE combine FGM and STAFF data series
%    b=C_FGM_STAFF_COMBINE(b_fgm,b_staff)      
%           currently implemented for burst mode data
% 

%% check input 
if nargout==0, help c_fgm_staff_combine;return;end

%% define common time interval in which to combine data sets
tint=[max(b_fgm(1,1),b_staff(1,1)) min(b_fgm(end,1),b_staff(end,1))]; % define common time interval
b_fgm=irf_tlim(b_fgm,tint);     % limit time series to common time interval
b_staff=irf_tlim(b_staff,tint); % limit time series to common time interval

%% check if FGM and STAFF are in normal mode
if 1/mean(diff(b_staff(1:min(end,100),1))) < 400 % STAFF in normal mode, not implemented
    disp('STAFF data in normal mode! FGM and STAFF combination not yet implemented!');
    b=[];
    return;
elseif 1/mean(diff(b_fgm(1:min(end,100),1))) < 60 % FGM in normal  mode
    disp('FGM data in normal mode! FGM and STAFF combination not yet implemented!');
    b=[];
    return;
end    
    
%% design filters
% low pass filter for FGM
    Fs = 66.7;  % Sampling Frequency
    
    Fpass = 2;           % Passband Frequency
    Fstop = 4;           % Stopband Frequency
    Apass = 1;           % Passband Ripple (dB)
    Astop = 20;          % Stopband Attenuation (dB)
    match = 'stopband';  % Band to match exactly
    
    % Construct an FDESIGN object and call its BUTTER method.
    ff  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
    d1_fgm = design(ff, 'butter', 'MatchExactly', match);
% higpass filter for STAFF
    Fs = 450;  % Sampling Frequency
    
    Fstop = 1;           % Stopband Frequency
    Fpass = 3;           % Passband Frequency
    Astop = 80;          % Stopband Attenuation (dB)
    Apass = 1;           % Passband Ripple (dB)
    match = 'stopband';  % Band to match exactly
    
    % Construct an FDESIGN object and call its BUTTER method.
    ff  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
    d1_staff = design(ff, 'butter', 'MatchExactly', match);

b_fgm_filt=b_fgm;
for j=2:4,
    b_fgm_filt(:,j)=filtfilt(d1_fgm.sosMatrix,d1_fgm.ScaleValues,b_fgm(:,j));
end
b_staff_filt=b_staff;
for j=2:4,
    b_staff_filt(:,j)=filtfilt(d1_staff.sosMatrix,d1_staff.ScaleValues,b_staff(:,j));
end
b=irf_add(1,b_staff_filt,1,b_fgm_filt); % FGM data are resampled at STAFF time line
if size(b,2)>4, % the 5th column should be |B|
    b(:,5)=irf_abs(b,1);
end
