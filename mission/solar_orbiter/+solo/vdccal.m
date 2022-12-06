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
%   PSP        - Probe-to-spacecraft potential
%   ScPot      - Spacecraft potential (PRELIMINARY PROXY)
%   codeVerStr - Date format version string for function itself.
%                This is used by BICAS for setting the relevant CDF global
%                attribute to indicate the version of the algorithm used to
%                produce a particular dataset (CDF file).
%   matVerStr  - Version string representing .mat file. Currently filename.
%                This is used by BICAS for setting the relevant CDF global
%                attribute to indicate the calibration data file used to produce
%                a particular dataset (CDF file).
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
    calfile_name = 'd23K123_20220124.mat';
else
    % Caller specified calibration file. Useful for debugging/testing new
    % calibrations.
    
    % (Do nothing.)
end
a = load(calfile_name);



%=============================================================
% Set return values that represent the version of calibration
%=============================================================
% Version of the function (not .mat file).
% NOTE: This value is meant to be be UPDATED BY HAND, not by an automatic
% timestamp, so that a constant value represents the same function/algorithm.
codeVerStr = '2021-10-21T12:00:00';
% Version of the .mat file. Using filename, or at least for now.
% This string is used by BICAS to set a CDF global attribute in official
% datasets for traceability.
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
allDiscontTimes=solo.ProbePotDiscontinuities;
%Time interval defined by the calibration file
calTint = irf.tint(a.d23.time(1),a.d23.time(end));
%We only care about discontinuities inside the calibration interval
discontTimes = EpochTT(allDiscontTimes.epoch(allDiscontTimes.epoch<=calTint.epoch(2)));

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
        % increase interval by 2 days before.
        subTint = irf.tint(tempTint(1)+(-2*24*60*60),next_discont);
    elseif isempty(next_discont)
        % If there are no discontinuities after the specified time,
        % increase interval by 2 days after.
        subTint = irf.tint(last_discont,tempTint(end)+(2*24*60*60));
    end

    %%
    VDC = VDC_inp.tlim(subTint);
    
    bSingleProbe = isnan(VDC.y.data) & isnan(VDC.z.data);
    
    % Resample calibration parameters
    d23R  = a.d23.tlim(subTint).resample(VDC);
    k23R = a.k23.tlim(subTint).resample(VDC);
    K123R = a.K123.tlim(subTint).resample(VDC);
    
    % Start calibration
    
    V1 = double(VDC.x.data);
    V2_scaled = double(VDC.y.data).*k23R.data +double(d23R.data); %Remove potential offset between 2,3
    V3 = double(VDC.z.data);
    V23 = (V2_scaled+V3)/2; % Corresponding to a measurement point between the two antennas.

    V23_scaled = (V23.*K123R.data(:,1) + K123R.data(:,2)); % Correcting V23 to V1
    
    PSP = irf.ts_scalar(VDC.time,(V23_scaled + V1)/2); % Compute PSP from corrected quantities.
    
    % Use alternate, simpler "calculation" for single-probe data.
    PSP.data(bSingleProbe) = VDC.x.data(bSingleProbe);
    PSP.units = 'V';
    PSP_out=PSP_out.combine(PSP);
    
    
    PLASMA_POT = 1.5; SHORT_FACTOR = 2.5; % XXX: these are just ad hoc numbers.
    
    ScPot = irf.ts_scalar(VDC.time, -(PSP.data-PLASMA_POT)*SHORT_FACTOR);
    ScPot.units = PSP.units;
    ScPot_out=ScPot_out.combine(ScPot);
    
    % Ey_SRF = V3 - V2, 6.99 - 1/2 of distance between the antennas
    V_d23 = V2_scaled-V3; %Fixed V2-V3.
    Ey_SRF = -V_d23*1e3/6.99;
    
    % Here we use the effective antenna length of 11.2 m, which correponds to
    % the distance between the center of ANT1 and a symmetric antenna on the
    % other side having voltage V23_scaled (see above).
    
    % For non scaled V23, the effective length would be L_123 = 6.97m, as shown in
    % Steinvall et al., 2021
    Ez_SRF = (V23_scaled - V1)*1e3/11.2;
    
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
