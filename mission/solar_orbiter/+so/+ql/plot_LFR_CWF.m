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
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2020-01-28.
%
function [hAxesArray] = plot_LFR_CWF(filePath)
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

    % TODO-DECISION: Content of figure title
    %   PROPOSAL: ~File path, ~filename
    %   PROPOSAL: Time range
    %   PROPOSAL: DATASET_ID
    %       NOTE: Varies depending on file since code covers three DATASET_IDs.
    %
    % TODO/BUG: Functioning ~color map
    % TODO: DC/AC: Detection, label.
    % TODO: Correct sampling frequency for irf_powerfft.
    % BUG: X axes differ between spectra and time series (irfu-matlab bug?).
    
    warning('Incomplete quicklook code')
    
    FILL_VALUE = single(-1e31);
    
    D = dataobj(filePath);
    
    Epoch = D.data.Epoch.data;
    V1_DC    = D.data.V.data(:,1);
    V12_DC   = D.data.E.data(:,1);
    %E13   = D.data.E.data(:,2);
    V23_DC   = D.data.E.data(:,3);
    V12_AC   = D.data.EAC.data(:,1);
    V23_AC   = D.data.EAC.data(:,3);
    
    V1_DC  = changem(V1_DC,  NaN, FILL_VALUE);
    V12_DC = changem(V12_DC, NaN, FILL_VALUE);
    V23_DC = changem(V23_DC, NaN, FILL_VALUE);
    V12_AC = changem(V12_AC, NaN, FILL_VALUE);
    V23_AC = changem(V23_AC, NaN, FILL_VALUE);
    
    irf_plot(5,'newfigure');
    
    TsV1_DC  = irf.ts_scalar(Epoch, V1_DC);
    TsV12_DC = irf.ts_scalar(Epoch, V12_DC);
    TsV23_DC = irf.ts_scalar(Epoch, V23_DC);
    TsV12_AC = irf.ts_scalar(Epoch, V12_AC);
    TsV23_AC = irf.ts_scalar(Epoch, V23_AC);
    
    hAxesArray = [];
    %h(end+1) = plot_spectrum( 'V1 spectrogram', TsV1,   'V1');
    %h(end+1) = plot_spectrum('E12 spectrogram', TsE12, 'V12');
    %h(end+1) = plot_spectrum('E23 spectrogram', TsE23, 'V23');
    
    hAxesArray(end+1) = plot_time_series( 'V1 DC time series', TsV1_DC,   'V1_DC [V]');
    hAxesArray(end+1) = plot_time_series('V12 DC time series', TsV12_DC, 'V12_DC [V]');
    hAxesArray(end+1) = plot_time_series('V23 DC time series', TsV23_DC, 'V23_DC [V]');
    hAxesArray(end+1) = plot_time_series('V12 AC time series', TsV12_AC, 'V12_AC [V]');
    hAxesArray(end+1) = plot_time_series('V23 AC time series', TsV23_AC, 'V23_AC [V]');

    so.ql.set_std_title('LFR CWF L2', filePath, hAxesArray(1))
    
    irf_plot_axis_align(hAxesArray)                     % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(hAxesArray, 'x', irf.tint(TsV1_DC.time))    % For aligning the content of the MATLAB axes.
    %title(hAxesArray(1), 'LFR CWF')
end



% panelTag      : 
% Ts            : irfumatlab TSeries (volt).
% yLabelNonUnit : y label without unit (unit is at the color bar; Assumes "Ts" uses volt).
function h = plot_spectrum(panelTag, Ts, yLabelNonUnit)
    %N_SAMPLES_PER_SPECTRUM = 2048;
    N_SAMPLES_PER_SPECTRUM = 128;
    
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
end



% panelTag      : 
% Ts            : irfumatlab TSeries
% yLabelNonUnit : y label with unit.
function hAxes = plot_time_series(panelTag, Ts, yLabel)
    hAxes = irf_panel(panelTag);
    irf_plot(hAxes, Ts)
    ylabel(hAxes, yLabel)
end
