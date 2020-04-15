%
% Quicklook for the content of one BIAS LFR CWF dataset (CDF file), ie. any of the following DATASET_IDs:
% SOLO_L2_RPW-LFR-SBM1-CWF-E
% SOLO_L2_RPW-LFR-SBM2-CWF-E
% SOLO_L2_RPW-LFR-SURV-CWF-E
%
%
% INCOMPLETE
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-28.
%
function hAxesArray = plot_LFR_CWF(filePath)
    % SOLO_L2_RPW-LFR-SBM1-CWF-E_V05.cdf zVariables:
    %
    % Variable Information (0 rVariable, 17 zVariables)
    % ===========================================================
    % Epoch                  CDF_TT2000/1      0:[]    T/
    % ACQUISITION_TIME       CDF_UINT4/1       1:[2]   T/T
    % ACQUISITION_TIME_UNITS CDF_CHAR/16       1:[2]   F/T
    % ACQUISITION_TIME_LABEL CDF_CHAR/32       1:[2]   F/T
    % QUALITY_FLAG           CDF_UINT1/1       0:[]    T/
    % QUALITY_BITMASK        CDF_UINT1/1       0:[]    T/
    % V_LABEL                CDF_CHAR/2        1:[3]   F/T
    % E_LABEL                CDF_CHAR/3        1:[3]   F/T
    % EAC_LABEL              CDF_CHAR/5        1:[3]   F/T
    % V                      CDF_REAL4/1       1:[3]   T/T
    % E                      CDF_REAL4/1       1:[3]   T/T
    % EAC                    CDF_REAL4/1       1:[3]   T/T
    % IBIAS1                 CDF_REAL4/1       0:[]    T/
    % IBIAS2                 CDF_REAL4/1       0:[]    T/
    % IBIAS3                 CDF_REAL4/1       0:[]    T/
    % DELTA_PLUS_MINUS       CDF_INT8/1        0:[]    T/
    % SYNCHRO_FLAG           CDF_UINT1/1       0:[]    T/
    %
    % NOTE: E_LABEL = ["V12","V13","V23"]
    
    % PROPOSAL: Define, "hardcode" array of function calls for each subplot (incl. arguments). For loop to execute calls.
    %   PRO: irf_plot knows how many subplots (before for loop).
    %   PRO: Collects axes handles.
    % PROPOSAL: Internal variable naming convention
    %   PROPOSAL: Follow zVars: VDC, EDC, EAC (not V*_DC etc)

    % TODO-DECISION: Content of figure title
    %   PROPOSAL: Time range
    %
    % TODO/BUG: Functioning ~color map
    % TODO: DC/AC: Detection, label.
    % TODO: Correct sampling frequency for irf_powerfft.
    % BUG: X axes differ between spectra and time series (irfu-matlab bug?).
    %
    % PROPOSAL: Function naming convention: Functions that basically create and fill an irf_panel should be named
    % ~panel.

    FILL_VALUE = single(-1e31);
    
    D = dataobj(filePath);
    
    Epoch = D.data.Epoch.data;
    V1_DC    = D.data.V.data(:,1);
    E12_DC   = D.data.E.data(:,1);
    %V13_DC   = D.data.E.data(:,2);
    E23_DC   = D.data.E.data(:,3);
    E12_AC   = D.data.EAC.data(:,1);
    E23_AC   = D.data.EAC.data(:,3);
    
    V1_DC  = changem(V1_DC,  NaN, FILL_VALUE);
    E12_DC = changem(E12_DC, NaN, FILL_VALUE);
    E23_DC = changem(E23_DC, NaN, FILL_VALUE);
    E12_AC = changem(E12_AC, NaN, FILL_VALUE);
    E23_AC = changem(E23_AC, NaN, FILL_VALUE);
    
    TsV1_DC  = irf.ts_scalar(Epoch, V1_DC);
    TsE12_DC = irf.ts_scalar(Epoch, E12_DC);
    TsE23_DC = irf.ts_scalar(Epoch, E23_DC);
    TsE12_AC = irf.ts_scalar(Epoch, E12_AC);
    TsE23_AC = irf.ts_scalar(Epoch, E23_AC);
    
    irf_plot(3+5,'newfigure');
    
    hAxesArray = [];
    hAxesArray(end+1) = plot_spectrum( 'V1 DC spectrogram', TsV1_DC,   'V1');
    hAxesArray(end+1) = plot_spectrum('V12 DC spectrogram', TsE12_DC, 'V12');
    hAxesArray(end+1) = plot_spectrum('V23 DC spectrogram', TsE23_DC, 'V23');

    hAxesArray(end+1) = plot_time_series( 'V1 DC time series', TsV1_DC,   'V1_DC [V]');
    hAxesArray(end+1) = plot_time_series('V12 DC time series', TsE12_DC, 'V12_DC [V]');
    hAxesArray(end+1) = plot_time_series('V23 DC time series', TsE23_DC, 'V23_DC [V]');
    hAxesArray(end+1) = plot_time_series('V12 AC time series', TsE12_AC, 'V12_AC [V]');
    hAxesArray(end+1) = plot_time_series('V23 AC time series', TsE23_AC, 'V23_AC [V]');

    solo.ql.set_std_title('LFR CWF L2', filePath, hAxesArray(1))
    
    irf_plot_axis_align(hAxesArray)                      % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(hAxesArray, 'x', irf.tint(TsV1_DC.time))    % For aligning the content of the MATLAB axes.
end



% panelTag      : 
% Ts            : irfumatlab TSeries
% yLabelNonUnit : y label with unit.
function hAxes = plot_time_series(panelTag, Ts, yLabel)
    hAxes = irf_panel(panelTag);
    irf_plot(hAxes, Ts)
    ylabel(hAxes, yLabel)
end



% ARGUMENTS
% =========
% panelTag      : 
% Ts            : irfumatlab TSeries (volt).
% yLabelNonUnit : y label without unit (unit is at the color bar; Assumes "Ts" uses volt).
%
function h = plot_spectrum(panelTag, Ts, yLabelNonUnit)
    
    nSamples = numel(Ts.time.epoch);
    fprintf('nSamples = %g\n', nSamples)    
    t = tic;

    N_SAMPLES_PER_SPECTRUM = 2048;
    %N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.
    
    % TEMPORARY FIX?!
    % Assuming there is one sampling frequency.
    % (There is no zVar for sampling frequency.)
    samplingFreqHz = 1/mode(diff(Ts.time.tts));
    
    h = irf_panel(panelTag);
    %irf_plot(hE12Spectra, 'colorbarlabel')
    specrec = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz);
    %irf_spectrogram(hE12Spectra,specrec);
    %irf_colormap(hE12Spectra,'default');

    specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
    %irf_plot(h, specrec);
    irf_spectrogram(h, specrec);   % Replaces irf_plot
    set(h, 'yscale','log')
    ylabel(h, {yLabelNonUnit;'f [Hz]'})   % NOTE: Adding frequency unit on separate row.
    
    tSec = toc(t);
    fprintf('tSec/nSamples = %g [s/sample]\n', tSec/nSamples)
end
