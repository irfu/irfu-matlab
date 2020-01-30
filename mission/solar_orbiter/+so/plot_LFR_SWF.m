%
% Quicklook for the content of one BIAS LFR SWF dataset (CDF file), i.e. DATASET_ID = SOLO_L2_RPW-LFR-SURV-SWF-E
%
% NOTE: Uses bicas.proc_utils.* code.
%
%
% INCOMPLETE
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2020-01-28.
%
function hAxes = plot_LFR_SWF(filePath)
    % SOLO_L2_RPW-LFR-SURV-SWF-E_V05.cdf zVariables:
    %
    % Variable Information (0 rVariable, 18 zVariables)
    % ===========================================================
    % Epoch                 CDF_TT2000/1      0:[]    T/
    % ACQUISITION_TIME      CDF_UINT4/1       1:[2]   T/T
    % ACQUISITION_TIME_UNITS CDF_CHAR/16      1:[2]   F/T
    % ACQUISITION_TIME_LABEL CDF_CHAR/32      1:[2]   F/T
    % QUALITY_FLAG          CDF_UINT1/1       0:[]    T/
    % QUALITY_BITMASK       CDF_UINT1/1       0:[]    T/
    % V_LABEL               CDF_CHAR/2        1:[3]   F/T
    % E_LABEL               CDF_CHAR/3        1:[3]   F/T
    % EAC_LABEL             CDF_CHAR/5        1:[3]   F/T
    % V                     CDF_REAL4/1       2:[2048,3]      T/TT
    % E                     CDF_REAL4/1       2:[2048,3]      T/TT
    % EAC                   CDF_REAL4/1       2:[2048,3]      T/TT
    % IBIAS1                CDF_REAL4/1       1:[2048]        T/T
    % IBIAS2                CDF_REAL4/1       1:[2048]        T/T
    % IBIAS3                CDF_REAL4/1       1:[2048]        T/T
    % DELTA_PLUS_MINUS      CDF_INT8/1        1:[2048]        T/T
    % F_SAMPLE              CDF_REAL4/1       0:[]    T/
    % SYNCHRO_FLAG          CDF_UINT1/1       0:[]    T/
    
    % BUG: Snapshot spectrograms consist of very narrow stripes.
    % BUG: Snapshot spectrogram for F2 misses some snapshots.
    % TODO/BUG: Functioning ~color map
    % TODO: DC/AC: Detection, label.
    % TODO: Correct sampling frequency for irf_powerfft.
    % TODO: Legend connecting line color to V1, E12, V13.
    
    warning('Incomplete quicklook')
    
    N_SAMPLES_PER_SNAPSHOT = 2048;    % Number of Samples Per (LFR) Snapshot (SPS).
    
    D = dataobj(filePath);
    Epoch = D.data.Epoch.data;
    
    hAxes = irf_plot(3*3 + 3,'newfigure');
    
    
    V1  = D.data.V.data(:,:,1);
    E12 = D.data.E.data(:,:,1);
    E13 = D.data.E.data(:,:,2);
    F_SAMPLE = D.data.F_SAMPLE.data;
    
    % LFR sampling frequencies (F0-F3 is LFR's terminology).
    F0Hz = 24576;
    F1Hz =  4096;
    F2Hz =   256;
    %F3Hz =    16;
    % LFR SWF only uses F0-F2.

    bF0 = (F_SAMPLE == F0Hz);
    bF1 = (F_SAMPLE == F1Hz);
    bF2 = (F_SAMPLE == F2Hz);
    
    % NOTE: Using bicas.* code.
    EpochF0 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF0), N_SAMPLES_PER_SNAPSHOT, F_SAMPLE(bF0));
    V1F0    = bicas.proc_utils.convert_N_to_1_SPR_redistribute(V1( bF0, :));
    E12F0   = bicas.proc_utils.convert_N_to_1_SPR_redistribute(E12(bF0, :));
    E13F0   = bicas.proc_utils.convert_N_to_1_SPR_redistribute(E13(bF0, :));
    
    EpochF1 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF1), N_SAMPLES_PER_SNAPSHOT, F_SAMPLE(bF1));
    V1F1    = bicas.proc_utils.convert_N_to_1_SPR_redistribute(V1( bF1, :));
    E12F1   = bicas.proc_utils.convert_N_to_1_SPR_redistribute(E12(bF1, :));
    E13F1   = bicas.proc_utils.convert_N_to_1_SPR_redistribute(E13(bF1, :));
    
    EpochF2 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF2), N_SAMPLES_PER_SNAPSHOT, F_SAMPLE(bF2));
    V1F2    = bicas.proc_utils.convert_N_to_1_SPR_redistribute(V1( bF2, :));
    E12F2   = bicas.proc_utils.convert_N_to_1_SPR_redistribute(E12(bF2, :));
    E13F2   = bicas.proc_utils.convert_N_to_1_SPR_redistribute(E13(bF2, :));
    
    
    
    TsV1F0  = irf.ts_scalar(EpochF0, V1F0);
    TsE12F0 = irf.ts_scalar(EpochF0, E12F0);
    TsE13F0 = irf.ts_scalar(EpochF0, E13F0);
    
    TsV1F1  = irf.ts_scalar(EpochF1, V1F1);
    TsE12F1 = irf.ts_scalar(EpochF1, E12F1);
    TsE13F1 = irf.ts_scalar(EpochF1, E13F1);

    TsV1F2  = irf.ts_scalar(EpochF2, V1F2);
    TsE12F2 = irf.ts_scalar(EpochF2, E12F2);
    TsE13F2 = irf.ts_scalar(EpochF2, E13F2);
    
    TsF0 = irf.ts_scalar(EpochF0, [V1F0, E12F0, E13F0]);
    TsF1 = irf.ts_scalar(EpochF1, [V1F1, E12F1, E13F1]);
    TsF2 = irf.ts_scalar(EpochF2, [V1F2, E12F2, E13F2]);

    
    
    h = [];
    h(end+1) = plot_spectrum( 'V1 F0 spectrogram', TsV1F0,  F0Hz, 'H', 'V1');
    h(end+1) = plot_spectrum('E12 F0 spectrogram', TsE12F0, F0Hz, 'H', 'E12');
    h(end+1) = plot_spectrum('E13 F0 spectrogram', TsE13F0, F0Hz, 'H', 'E13');
    
    h(end+1) = plot_spectrum( 'V1 F1 spectrogram', TsV1F1,  F1Hz, 'M', 'V1');
    h(end+1) = plot_spectrum('E12 F1 spectrogram', TsE12F1, F1Hz, 'M', 'E12');
    h(end+1) = plot_spectrum('E13 F1 spectrogram', TsE13F1, F1Hz, 'M', 'E13');
    
    h(end+1) = plot_spectrum( 'V1 F2 spectrogram', TsV1F2,  F2Hz, 'L', 'V1');
    h(end+1) = plot_spectrum('E12 F2 spectrogram', TsE12F2, F2Hz, 'L', 'E12');
    h(end+1) = plot_spectrum('E13 F2 spectrogram', TsE13F2, F2Hz, 'L', 'E13');
    
    h(end+1) = plot_time_series('V1,E12,E13 F0 time series', TsF0, 'H');
    h(end+1) = plot_time_series('V1,E12,E13 F1 time series', TsF1, 'M');
    h(end+1) = plot_time_series('V1,E12,E13 F2 time series', TsF2, 'L');

    irf_plot_axis_align(h)                     % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(h, 'x', irf.tint(TsV1F2.time))    % For aligning the content of the MATLAB axes.
    title(h(1), 'LFR SWF')
end



function h = plot_spectrum(panelTag, Ts, samplingFreqHz, tlLegend, trLegend)
    N_SAMPLES_PER_SNAPSHOT = 2048;    % Number of Samples Per (LFR) Snapshot (SPS).
    N_SAMPLES_PER_SPECTRUM = N_SAMPLES_PER_SNAPSHOT/4;

    h = irf_panel(panelTag);
    specrec = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz);
    specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
    %irf_plot(h, specrec, 'colorbarlabel', {'',''}, 'fitcolorbarlabel');
    irf_spectrogram(h, specrec);   % Replaces irf_plot    
    set(h, 'yscale','log')

    % Two rows of label causes trouble for the time series ylabels.
    %ylabel(h, {yLabelNonUnit;'f [Hz]'})   % NOTE: Adding frequency unit on separate row.
    %ylabel(h, {yLabelNonUnit})   % NOTE: Adding frequency unit on separate row.

    irf_legend(h, {tlLegend}, [0.02 0.98], 'color', 'k')
    irf_legend(h, trLegend , [0.98 0.98])
end



function h = plot_time_series(panelTag, Ts, tlLegend)
    SIGNALS_LEGEND = {'V1','E12','E13'};
    
    h = irf_panel(panelTag);
    irf_plot(h, Ts)
    ylabel(h, '[V]')
    irf_legend(h, {tlLegend}, [0.02 0.98], 'color', 'k')
    irf_legend(h, SIGNALS_LEGEND, [0.98 0.98])
end