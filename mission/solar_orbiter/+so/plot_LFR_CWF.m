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
function hAxes = plot_LFR_CWF(filePath)
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
    
    warning('Incomplete quicklook')
    
    D = dataobj(filePath);
    
    Epoch = D.data.Epoch.data;
    V1    = D.data.V.data(:,1);
    E12   = D.data.E.data(:,1);
    E13   = D.data.E.data(:,2);
    
    hAxes = irf_plot(6,'newfigure');
    
    TsV1  = irf.ts_scalar(Epoch, V1);
    TsE12 = irf.ts_scalar(Epoch, E12);
    TsE13 = irf.ts_scalar(Epoch, E13);
    
    h = [];
    h(end+1) = plot_spectrum('E12 spectrogram', TsE12, 'E12');
    h(end+1) = plot_spectrum('E13 spectrogram', TsE13, 'E13');
    h(end+1) = plot_spectrum( 'V1 spectrogram', TsV1,   'V1');
    
    h(end+1) = plot_time_series('E12 time series', TsE12, 'E12 [V]');
    h(end+1) = plot_time_series('E13 time series', TsE13, 'E13 [V]');
    h(end+1) = plot_time_series( 'V1 time series', TsV1,   'V1 [V]');
    
    irf_plot_axis_align(h)                   % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(h, 'x', irf.tint(TsV1.time))    % For aligning the content of the MATLAB axes.
    title(h(1), 'LFR CWF')
end



function h = plot_spectrum(panelTag, Ts, yLabelNonUnit)
    %N_SAMPLES_PER_SPECTRUM = 2048;
    N_SAMPLES_PER_SPECTRUM = 512;
    
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



function h = plot_time_series(panelTag, Ts, yLabel)
    h = irf_panel(panelTag);
    irf_plot(h, Ts)
    ylabel(h, yLabel)
end