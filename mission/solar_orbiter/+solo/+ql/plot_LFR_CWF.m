%
% Quicklook for the content of one BIAS LFR CWF dataset (CDF file), ie. any of the following DATASET_IDs:
% SOLO_L2_RPW-LFR-SBM1-CWF-E
% SOLO_L2_RPW-LFR-SBM2-CWF-E
% SOLO_L2_RPW-LFR-SURV-CWF-E
%
%
% NOTE: Does not properly derive sampling frequency. (Awaiting new dataset skeletons.)
% NOTE: Color scale is log, therefore negative values (probably).
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
    
    % TODO-DECISION: Content of figure title
    %   PROPOSAL: Time range
    %   PROPOSAL: DOY.
    %       ~PROBLEM/NOTE: Could span multiple days, or part of day.
    %
    % TODO-NI: Detect DC/AC and only plot one of them? Analogous to LFR SWF.
    %   TODO-NI: CWF never uses AC? SBM1/2?
    % TODO-NI: Spectrum overlap 50% also for CWF?
    % TODO: Correct sampling frequency needed for irf_powerfft.
    % TODO: Use zVar SAMPLING_RATE (when introduced) for spectrums.
    % PROPOSAL: Settings arguments to disable/enable hidden functionality
    %   Ex: SWF: Permit/force DC+AC diffs
    %   Ex: Disable spectrograms.
    %
    % PROPOSAL: Change package name to sp (summary plots).
    %
    % ~BUG: Can probably not handle new zVariable names VDC.
    %
    % PROPOSAL: Settings for disabling spectrum etc.
    
    %PERMIT_SIMULTANEOUS_DC_AC_DIFFS = 0;

    D = dataobj(filePath);
    
    epoch = D.data.Epoch.data;
    vDc1  = get_CDF_zv_data(D, 'V',   1);
    vDc12 = get_CDF_zv_data(D, 'E',   1);
    vDc23 = get_CDF_zv_data(D, 'E',   3);
    vAc12 = get_CDF_zv_data(D, 'EAC', 1);
    vAc23 = get_CDF_zv_data(D, 'EAC', 3);

    TsVdc1  = irf.ts_scalar(epoch, vDc1);
    TsVdc12 = irf.ts_scalar(epoch, vDc12);
    TsVdc23 = irf.ts_scalar(epoch, vDc23);
    TsVac12 = irf.ts_scalar(epoch, vAc12);
    TsVac23 = irf.ts_scalar(epoch, vAc23);
    
    % TEMPORARY FIX
    % Assume there is one sampling frequency + datagaps.
    % There is no zVar for sampling frequency in old datasets.
    samplingFreqHz = 1/mode(diff(TsVdc1.time.tts));

    irf_plot(3+5,'newfigure');
    
    hAxesArray = [];
    hAxesArray(end+1) = spectrogram_panel( 'V1 DC spectrogram', TsVdc1,  samplingFreqHz, 'V1\_DC');
    hAxesArray(end+1) = spectrogram_panel('V12 DC spectrogram', TsVdc12, samplingFreqHz, 'V12\_DC');
    hAxesArray(end+1) = spectrogram_panel('V23 DC spectrogram', TsVdc23, samplingFreqHz, 'V23\_DC');

    hAxesArray(end+1) = time_series_panel( 'V1 DC time series', TsVdc1,   'V1_DC [V]');
    hAxesArray(end+1) = time_series_panel('V12 DC time series', TsVdc12, 'V12_DC [V]');
    hAxesArray(end+1) = time_series_panel('V23 DC time series', TsVdc23, 'V23_DC [V]');
    hAxesArray(end+1) = time_series_panel('V12 AC time series', TsVac12, 'V12_AC [V]');
    hAxesArray(end+1) = time_series_panel('V23 AC time series', TsVac23, 'V23_AC [V]');

    % TEST
%     TILE_MRG_SETTINGS = struct(...
%         'mrgLeft',   10, ...
%         'mrgRight',  10, ...
%         'mrgTop',    5, ...
%         'mrgBottom', 5);
%     EJ_library.graph.tile_by_InnerPosition_OuterPosition(hAxesArray(:), TILE_MRG_SETTINGS)
            
    solo.ql.set_std_title('LFR CWF L2', filePath, hAxesArray(1))
    
    irf_plot_axis_align(hAxesArray)                     % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(hAxesArray, 'x', irf.tint(TsVdc1.time))    % For aligning the content of the MATLAB axes.
end



% ARGUMENTS
% =========
% panelTag      : 
% Ts            : irfumatlab TSeries (volt).
% yLabelNonUnit : y label without unit (unit is at the color bar; Assumes "Ts" uses volt).
%
function h = spectrogram_panel(panelTag, Ts, samplingFreqHz, yLabelNonUnit)
    
    %N_SAMPLES_PER_SPECTRUM = 2048;
    N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.
    
    h = irf_panel(panelTag);
    %irf_plot(hE12Spectra, 'colorbarlabel')
    Specrec = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz);
    %irf_spectrogram(hE12Spectra,specrec);
    %irf_colormap(hE12Spectra,'default');

    Specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
    %irf_plot(h, specrec);
    irf_spectrogram(h, Specrec);   % Replaces irf_plot
    set(h, 'yscale','log')
    ylabel(h, {yLabelNonUnit; 'f [Hz]'})   % NOTE: Adding frequency unit on separate row.
end



% panelTag      : 
% Ts            : irfumatlab TSeries
% yLabelNonUnit : y label with unit.
function hAxes = time_series_panel(panelTag, Ts, yLabel)
    hAxes = irf_panel(panelTag);
    irf_plot(hAxes, Ts)
    ylabel(hAxes, yLabel)
end



function data = get_CDF_zv_data(D, zvName, i2)
    fillValue = getfillval(D, zvName);
    data = D.data.(zvName).data(:, i2);
    data = changem(data, NaN, fillValue);
end
