function [DCE_SRF_out,PSP_out,ScPot_out,codeVerStr,matVerStr] = vdccal(VDC_inp,EDC_inp,calfile_name)
%SOLO.VDCCAL  Calibrate VDC to get DC E and PSP
%
%    [DCE_SRF,PSP,ScPot,codeVerStr,matVerStr] = solo.vdccal(VDC,EDC,calfilename)
%
% Inputs: VDC,EDC from L2 CWF files, 
%         calfile_name is the name of the calibration file (.mat) one wishes to
%         use. If calfile_name is empty, the code used the calibration file
%         used by BICAS for producing official datasets
%
% Outputs:
%   DCE_SRF    - DC electric field in SRF (Ex=0)
%   PSP        - probe-to-spacecraft potential
%   ScPot      - spacecraft potential (PRELIMINARY PROXY)
%   codeVerStr - Date format version string for function itself. Used by BICAS.
%   matVerStr  - Version string representing .mat file. Currently filename.
%                Used by BICAS.
%
% Loads .mat file produced by solo.correlate_probes_batch (script) if
% calfile_name is empty.
% NOTE: .mat needs to be updated before processing the new data.
%
% NOTE: This function is used by BICAS for producing official L3 datasets.



% Act depending on whether a calibration file is specified or not.
if isempty(calfile_name)
    % Caller did not specify calibration file.
    % IMPORTANT: USES CALIBRATION FILE THAT IS USED BY BICAS FOR PRODUCING
    % OFFICIAL DATASETS.
    calfile_name = 'd23K123_20210521.mat';
else
    % Caller specified calibration file. Useful for debugging/testing new
    % calibrations.
    
    % (Do nothing.)
end
a = load(calfile_name);



%===========================================================================
% Date strings that represent the version of calibration. These strings are
% used by BICAS to set a CDF global attribute in official datasets for
% traceability.
% --
% Version of the function (not .mat file).
% NOTE: This value is meant to be be updated by hand, not by an automatic
% timestamp, so that a constant value represents the same function.
%===========================================================================
codeVerStr = '2021-05-25T12:00:00';
% Version of the .mat file. Using filename, or at least for now.
[~, basename, suffix] = fileparts(calfile_name);
matVerStr  = [basename, suffix];   % Only use filename, not entire path.

%=============================================================================
% Find data points/CDF records for which only V1_DC is available
%
% NOTE: solo.vdccal() only uses DC data (not AC). For mux=0, there are two
% cases:
% (1) all single probes (DC) available
% (2) only probe 1 (DC) available.
% NOTE: Only works for mux=0,2,3,4 (not mux=1).
% NOTE: Currently ignores EDC_inp argument.
%==============================================================================

%load Probe-potential discontinuities
discontTimes=solo.ProbePotDiscontinuities;
mainTint = irf.tint(VDC_inp.time(1),VDC_inp.time(end));
sub_int_times = EpochTT(solo.split_tint(mainTint,discontTimes));

%Predefine output variables:
DCE_SRF_out=irf.ts_vec_xyz(EpochTT([]),double.empty(0,3));
PSP_out=irf.ts_scalar(EpochTT([]),[]);
ScPot_out=irf.ts_scalar(EpochTT([]),[]);

% Perform calibration on each subinterval (if any probe-to-spacecraft
% potential discontinuities are present)

for isub=1:length(sub_int_times)-1
    
    tempTint=sub_int_times(isub:isub+1);
    % Find the closest discontinuities.
    last_discont = EpochTT(max(discontTimes.epoch(tempTint(1).epoch>=discontTimes.epoch)));
    next_discont = EpochTT(min(discontTimes.epoch(tempTint(end).epoch<=discontTimes.epoch)));
 
    % Extend the time interval to the closest discontinuities (this is to
    % avoid problems that occur if VDC contains very little data).
    if ~isempty(last_discont) && ~isempty(next_discont)
        subTint=irf.tint(last_discont,next_discont);
    elseif isempty(last_discont)
        % If there are no discontinuities before the specified time,
        % increase interval by 10 hours before.
        subTint = irf.tint(tempTint(1)+(-10*60*60),next_discont);
    elseif isempty(next_discont)
        % If there are no discontinuities after the specified time,
        % increase interval by 10 hours after.
        subTint = irf.tint(last_discont,tempTint(end)+(10*60*60));
    end

    %%
    VDC = VDC_inp.tlim(subTint);
    
    bSingleProbe = isnan(VDC.y.data) & isnan(VDC.z.data);
    
    % Resample calibration parameters
    d23R  = a.d23.tlim(subTint).resample(VDC);
    K123R = a.K123.tlim(subTint).resample(VDC);
    Gamma0R = a.Gamma0.tlim(subTint).resample(VDC);
    Gamma1R = a.Gamma1.tlim(subTint).resample(VDC);
    ccR = a.CC.tlim(subTint).resample(VDC);
    
    % Start calibration
    V2corr = double(VDC.y.data) -double(d23R.data); %Remove potential offset between 2,3
    V23_corr = (V2corr+double(VDC.z.data))/2; %(V2corr+V3)/2
    V2cmr = double(V2corr)-(Gamma0R.data+V23_corr.*Gamma1R.data)/2; %Remove common mode from V2.
    V3cmr = double(VDC.z.data)+(Gamma0R.data+V23_corr.*Gamma1R.data)/2;
    
    V23 = (V2cmr + V3cmr)/2; % (V2cmr + V3) /2
    
    V23corr = (V23.*K123R.data(:,1) + K123R.data(:,2)); %Correcting V23 to V1
    
    PSP = irf.ts_scalar(VDC.time,(V23corr + double(VDC.x.data))/2); %Compute PSP from corrected quantities.
    
    % Use alternate, simpler "calculation" for single-probe data.
    PSP.data(bSingleProbe) = VDC.x.data(bSingleProbe);
    PSP.units = 'V';
    PSP_out=PSP_out.combine(PSP);
    
    
    PLASMA_POT = 1.5; SHORT_FACTOR = 2.5; % XXX: these are just ad hoc numbers.
    
    ScPot = irf.ts_scalar(VDC.time, -(PSP.data-PLASMA_POT)*SHORT_FACTOR);
    ScPot.units = PSP.units;
    ScPot_out=ScPot_out.combine(ScPot);
    
    % Ey_SRF = V3 - V2, 6.99 - 1/2 of distance between the antennas
    V_delta23_corr = V2cmr-V3cmr; %Fixed V2-V3.
    Ey_SRF = -V_delta23_corr*1e3/6.99;
    
    % Ez_SRF = V23corr - V1
    % Here we use the antenna length of 11.2 m, which correponds to
    % the distance between the center of ANT1 and a symmetric antenna on the
    % other side having voltage V23 corr.
    Ez_SRF = V23corr - double(VDC.x.data);
    Ez_SRF = Ez_SRF*1e3/11.2;
    
    DCE_SRF = irf.ts_vec_xyz(VDC.time,[Ey_SRF*0 Ey_SRF Ez_SRF]);
    DCE_SRF.units = 'mV/m';
    DCE_SRF.coordinateSystem = 'SRF';
    
    DCE_SRF_out=DCE_SRF_out.combine(DCE_SRF);
    
end %for
    % Specify units and coordinate system
    DCE_SRF_out.units = 'mV/m';
    DCE_SRF_out.coordinateSystem = 'SRF';
    PSP_out.units = 'V';
    ScPot_out.units=PSP_out.units;
    
end %function
