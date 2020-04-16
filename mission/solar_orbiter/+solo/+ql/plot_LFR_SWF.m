%
% Quicklook for the content of one BIAS LFR SWF dataset (CDF file), i.e. DATASET_ID = SOLO_L2_RPW-LFR-SURV-SWF-E
%
%
% NOTE: Only capable (default) of only showing either DC diffs or AC diffs. Hardcoded setting for permitting both.
% NOTE: Uses bicas.proc_utils.* code.
%
% NOTE: Does not yet support spectrogram overlap.
%
% INCOMPLETE
%
%
% ARGUMENTS
% =========
% timeIntervUtc : Optional. Cell array, length 2. {1},{2} = UTC string, begin & end.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-28.
%
function hAxesArray = plot_LFR_SWF(filePath, timeIntervUtc)
    %
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

    % POLICY: BOGIQ for all quicklook plot code: See plot_LFR_CWF.
    %
    % PROPOSAL: Document "proper speed test" since spectrograms are slow.
    %   PROPOSAL: Separate file for this.
    %   PROPOSAL: Vary nSamplesPerSpectrum.
    %
    % TODO-DECISION: How submit data to spectrum_panel?
    %   NEED: Should also work for TDS's varying-length snapshots?
    %   NOTE: Might be necessary to also split up CWF to speed up spectrograms eventually.
    %   TODO-DECISION: Should caller split up data into snapshots?
    %   PROPOSAL: Submit one big TSeries and have the function split up the data based on gaps.
    %   PROPOSAL: Submit one TSeries per snapshot.
    %       PRO: Works for TDS snapshots.
    %       PRO: Splits up the 
    %   PROPOSAL: Submit ~zVar with one snapshot per row.
    %       CON: Does not work per snapshot.
    %   PROPOSAL: Always let the caller split up the time series in segments, snapshots or otherwise.
    %
    % TODO?: Remove interpolation between time series snapshots?
    % TODO: 50% overlap PSD
    
    PERMIT_SIMULTANEOUS_DC_AC_DIFFS = 0;   % YK 2020-04-16: Officially only either DC or AC diffs.
    
    % LFR sampling frequencies (F0-F3 is LFR's terminology).
    % NOTE: LFR SWF only uses F0-F2 (not F3).
    F0Hz = 24576;
    F1Hz =  4096;
    F2Hz =   256;

    
    
    % Interpret argument for optionally restricting covered time range.
    % TEMPORARY? Useful for debugging when plotting large files (slow plotting).
    if nargin == 1
        DATAOBJ_TIME_INTERVAL_ARGS = {};
    elseif nargin == 2
        assert(iscell(timeIntervUtc))
        DATAOBJ_TIME_INTERVAL_ARGS = {'tint', [iso2epoch(timeIntervUtc{1}), iso2epoch(timeIntervUtc{2})]};
    end
    
    

    D = dataobj(filePath, DATAOBJ_TIME_INTERVAL_ARGS{:});
    
    Epoch    = D.data.Epoch.data;
    F_SAMPLE = D.data.F_SAMPLE.data;
    VDC1     = get_CDF_zv_data(D, 'V',   1);
    VDC12    = get_CDF_zv_data(D, 'E',   1);
    VDC23    = get_CDF_zv_data(D, 'E',   3);
    VAC12    = get_CDF_zv_data(D, 'EAC', 1);
    VAC23    = get_CDF_zv_data(D, 'EAC', 3);
    
    hasDcDiffData = any(~isnan(VDC12(:))) || any(~isnan(VDC23(:)));
    hasAcDiffData = any(~isnan(VAC12(:))) || any(~isnan(VAC23(:)));
    if hasDcDiffData && hasAcDiffData
        error('Dataset (CDF file) contains both DC diff and AC diff data. Can not handle this case.')
    elseif ~hasDcDiffData && ~hasAcDiffData && ~PERMIT_SIMULTANEOUS_DC_AC_DIFFS
        error('Dataset (CDF file) contains neither DC diff nor AC diff data. Can not handle this case.')
    end
    displayDcDiffData = hasDcDiffData;
    displayAcDiffData = hasAcDiffData;
    
    
    
    % B = Boolean/Logical (true/false for every index value).
    bF0 = (F_SAMPLE == F0Hz);
    bF1 = (F_SAMPLE == F1Hz);
    bF2 = (F_SAMPLE == F2Hz);



%     nSps = size(VDC1, 2);   % SPS = Samples Per Snapshot
%     
%     % NOTE: Using bicas.* code.
%     EpochF0  = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF0), nSps, F_SAMPLE(bF0));
%     %VDC1(:,end) = NaN;   % TEST
%     Vdc1F0  = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC1( bF0, :));
%     Vdc12F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC12(bF0, :));
%     Vdc23F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC23(bF0, :));
%     Vac12F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC12(bF0, :));
%     Vac23F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC23(bF0, :));
%     
%     EpochF1 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF1), nSps, F_SAMPLE(bF1));
%     Vdc1F1  = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC1( bF1, :));
%     Vdc12F1 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC12(bF1, :));
%     Vdc23F1 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC23(bF1, :));
%     Vac12F1 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC12(bF1, :));
%     Vac23F1 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC23(bF1, :));
%     
%     EpochF2 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF2), nSps, F_SAMPLE(bF2));
%     Vdc1F2  = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC1( bF2, :));
%     Vdc12F2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC12(bF2, :));
%     Vdc23F2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC23(bF2, :));
%     Vac12F2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC12(bF2, :));
%     Vac23F2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC23(bF2, :));
% 
    % Create TSeries representing 3 scalar time series each.
%     TsVdcF0 = irf.ts_scalar(EpochF0, [Vdc1F0, Vdc12F0, Vdc23F0]);
%     TsVdcF1 = irf.ts_scalar(EpochF1, [Vdc1F1, Vdc12F1, Vdc23F1]);
%     TsVdcF2 = irf.ts_scalar(EpochF2, [Vdc1F2, Vdc12F2, Vdc23F2]);
    % Create TSeries representing 3 scalar time series each.
%     TsVdcacF0 = irf.ts_scalar(EpochF0, [Vdc1F0, Vac12F0, Vac23F0]);
%     TsVdcacF1 = irf.ts_scalar(EpochF1, [Vdc1F1, Vac12F1, Vac23F1]);
%     TsVdcacF2 = irf.ts_scalar(EpochF2, [Vdc1F2, Vac12F2, Vac23F2]);
%     % Create TSeries representing 2 scalar time series each.
%     TsVacF0 = irf.ts_scalar(EpochF0, [        Vac12F0, Vac23F0]);
%     TsVacF1 = irf.ts_scalar(EpochF1, [        Vac12F1, Vac23F1]);
%     TsVacF2 = irf.ts_scalar(EpochF2, [        Vac12F2, Vac23F2]);



    pcfcList = {};    % PCFC = Panel Creation Function Call
    %=================
    % F0 spectrograms
    %=================
    pcfcList{end+1}     = @() (spectrum_panel( 'V1 DC F0 spectrogram', Epoch(bF0), VDC1( bF0, :), F0Hz, 'F0', 'V1\_DC'));
    if displayDcDiffData
        pcfcList{end+1} = @() (spectrum_panel('V12 DC F0 spectrogram', Epoch(bF0), VDC12(bF0, :), F0Hz, 'F0', 'V12\_DC'));
        pcfcList{end+1} = @() (spectrum_panel('V23 DC F0 spectrogram', Epoch(bF0), VDC23(bF0, :), F0Hz, 'F0', 'V13\_DC'));
    end
    if displayAcDiffData
        pcfcList{end+1} = @() (spectrum_panel('V12 AC F0 spectrogram', Epoch(bF0), VAC12(bF0, :), F0Hz, 'F0', 'V12\_AC'));
        pcfcList{end+1} = @() (spectrum_panel('V23 AC F0 spectrogram', Epoch(bF0), VAC23(bF0, :), F0Hz, 'F0', 'V23\_AC'));
    end
    if 1   % DEBUG
    %=================
    % F1 spectrograms
    %=================
    pcfcList{end+1} =     @() (spectrum_panel( 'V1 DC F1 spectrogram', Epoch(bF1), VDC1( bF1, :), F1Hz, 'F1', 'V1\_DC'));
    if displayDcDiffData
        pcfcList{end+1} = @() (spectrum_panel('V12 DC F1 spectrogram', Epoch(bF1), VDC12(bF1, :), F1Hz, 'F1', 'V12\_DC'));
        pcfcList{end+1} = @() (spectrum_panel('V23 DC F1 spectrogram', Epoch(bF1), VDC23(bF1, :), F1Hz, 'F1', 'V23\_DC'));
    end
    if displayAcDiffData
        pcfcList{end+1} = @() (spectrum_panel('V12 AC F1 spectrogram', Epoch(bF1), VAC12(bF1, :), F1Hz, 'F1', 'V12\_AC'));
        pcfcList{end+1} = @() (spectrum_panel('V23 AC F1 spectrogram', Epoch(bF1), VAC23(bF1, :), F1Hz, 'F1', 'V23\_AC'));
    end
    %=================
    % F2 spectrograms
    %=================
    pcfcList{end+1}     = @() (spectrum_panel( 'V1 DC F2 spectrogram', Epoch(bF2), VDC1( bF2, :), F2Hz, 'F2', 'V1\_DC'));
    if displayDcDiffData
        pcfcList{end+1} = @() (spectrum_panel('V12 DC F2 spectrogram', Epoch(bF2), VDC12(bF2, :), F2Hz, 'F2', 'V12\_DC'));
        pcfcList{end+1} = @() (spectrum_panel('V23 DC F2 spectrogram', Epoch(bF2), VDC23(bF2, :), F2Hz, 'F2', 'V23\_DC'));
    end
    if displayAcDiffData
        pcfcList{end+1} = @() (spectrum_panel('V12 AC F2 spectrogram', Epoch(bF2), VAC12(bF2, :), F2Hz, 'F2', 'V12\_AC'));
        pcfcList{end+1} = @() (spectrum_panel('V23 AC F2 spectrogram', Epoch(bF2), VAC23(bF2, :), F2Hz, 'F2', 'V23\_AC'));
    end
    end    % DEBUG
    %===================
    % F0-F2 time series 
    % IMPLEMENTATION NOTE: Panel tags have to be unique, or otherwise they will be reused.
    %===================
    if displayDcDiffData

        % DC single + DC diffs
        SIGNALS_LEGEND_DC = EJ_library.graph.escape_str({'V1_DC','V12_DC','V23_DC'});
        tempFuncPtr = @(bFx, freqStr, samplingFreqHz) (@() (time_series_panel2(...
            sprintf('V1,V12,V23 DC %s time series', freqStr), ...
            Epoch(bFx), {VDC1(bFx, :), VDC12(bFx, :), VDC23(bFx, :)}, ...
            samplingFreqHz, freqStr, SIGNALS_LEGEND_DC)));% 
        pcfcList{end+1} = tempFuncPtr(bF0, 'F0', F0Hz);
        pcfcList{end+1} = tempFuncPtr(bF1, 'F1', F1Hz);
        pcfcList{end+1} = tempFuncPtr(bF2, 'F2', F2Hz);
%         pcfcList{end+1} = @() (time_series_panel2('V1,V12,V23 DC F0 time series', EpochF0, {VDC1(bF0, :),VDC12(bF0, :),VDC23(bF0, :)}, 'F0', SIGNALS_LEGEND_DC));
%         pcfcList{end+1} = @() (time_series_panel2('V1,V12,V23 DC F1 time series', EpochF1, {VDC1(bF1, :),VDC12(bF1, :),VDC23(bF1, :)}, 'F1', SIGNALS_LEGEND_DC));
%         pcfcList{end+1} = @() (time_series_panel2('V1,V12,V23 DC F2 time series', EpochF2, {VDC1(bF2, :),VDC12(bF2, :),VDC23(bF2, :)}, 'F2', SIGNALS_LEGEND_DC));

    end    
    if ~displayDcDiffData && displayAcDiffData

        % DC single + AC diffs
        SIGNALS_LEGEND_DC_AC = EJ_library.graph.escape_str({'V1_DC','V12_AC','V23_AC'});
        tempFuncPtr = @(bFx, freqStr, samplingFreqHz) (@() (time_series_panel2(...
            sprintf('V1,V12,V23 DC/AC %s time series', freqStr), ...
            Epoch(bFx), {VDC1(bFx, :), VAC12(bFx, :), VAC23(bFx, :)}, ...
            samplingFreqHz, freqStr, SIGNALS_LEGEND_DC_AC)));

        pcfcList{end+1} = tempFuncPtr(bF0, 'F0', F0Hz);
        pcfcList{end+1} = tempFuncPtr(bF1, 'F1', F1Hz);
        pcfcList{end+1} = tempFuncPtr(bF2, 'F2', F2Hz);
%         pcfcList{end+1} = @() (time_series_panel('V1,V12,V23 DC/AC F0 time series', TsVdcacF0, 'F0', SIGNALS_LEGEND_DC_AC));
%         pcfcList{end+1} = @() (time_series_panel('V1,V12,V23 DC/AC F1 time series', TsVdcacF1, 'F1', SIGNALS_LEGEND_DC_AC));
%         pcfcList{end+1} = @() (time_series_panel('V1,V12,V23 DC/AC F2 time series', TsVdcacF2, 'F2', SIGNALS_LEGEND_DC_AC));

    end    
    if displayDcDiffData && displayAcDiffData

        % AC diffs
        SIGNALS_LEGEND_AC = EJ_library.graph.escape_str({'V12_AC','V23_AC'});
        tempFuncPtr = @(bFx, freqStr, samplingFreqHz) (@() (time_series_panel2(...
            sprintf('V12,V23 AC %s time series', freqStr), ...
            Epoch(bFx), {VAC12(bFx, :), VAC23(bFx, :)}, ...
            samplingFreqHz, freqStr, SIGNALS_LEGEND_AC)));

        pcfcList{end+1} = tempFuncPtr(bF0, 'F0', F0Hz);
        pcfcList{end+1} = tempFuncPtr(bF1, 'F1', F1Hz);
        pcfcList{end+1} = tempFuncPtr(bF2, 'F2', F2Hz);
%         pcfcList{end+1} = @() (time_series_panel2('V12,V23 AC F0 time series', TsVacF0, 'F0', SIGNALS_LEGEND_AC));
%         pcfcList{end+1} = @() (time_series_panel2('V12,V23 AC F1 time series', TsVacF1, 'F1', SIGNALS_LEGEND_AC));
%         pcfcList{end+1} = @() (time_series_panel2('V12,V23 AC F2 time series', TsVacF2, 'F2', SIGNALS_LEGEND_AC));

    end

    %======
    % Plot
    %======
    irf_plot(numel(pcfcList), 'newfigure');
    hAxesArray = [];
    for i = 1:numel(pcfcList)
        funcPtr = pcfcList{i};
        hAxesArray(end+1) = funcPtr();
    end

    solo.ql.set_std_title('LFR SWF L2', filePath, hAxesArray(1))

    irf_plot_axis_align(hAxesArray)                      % For aligning MATLAB axes (taking color legends into account).
    %irf_zoom(hAxesArray, 'x', irf.tint(TsVdcF0.time))    % For aligning the content of the MATLAB axes.    
    irf_zoom(hAxesArray, 'x', irf.tint(Epoch(1), Epoch(end)))    % For aligning the content of the MATLAB axes.    
end



function data = get_CDF_zv_data(D, zvName, i3)
    fillValue = getfillval(D, zvName);
    data = D.data.(zvName).data(:, :, i3);
    data = changem(data, NaN, fillValue);
end



% ARGUMENTS
% =========
% TsCa     : Cell array of TSeries. All TSeries separately describe different time segments of one single scalar time
%            series. In practice, each TSeries should represent one snapshot. This speeds up the spectrograms by a lot.
% tlLegend : Top-left  (TL) legend string.
% trLegend : Top-right (TR) legend string.
%
function h = spectrum_panel(panelTag, zvEpoch, zvData, samplingFreqHz, tlLegend, trLegend)
    % NOTE: Multiple-row labels causes trouble for the time series ylabels.
    % IMPLEMENTATION NOTE: Implemented to be able to handle TDS snapshots that vary in length (in theory; untested).



    % Fraction of the (minimum) time distance between snapshots (centers) that will be used for displaying the spectra.
    % =1 : Spectras are adjacent between snapshot (for minimum snapshot distance).
    SNAPSHOT_WIDTH_FRACTION  = 0.85;
    SPECTRUM_OVERLAP_PERCENT = 0;  % Percent, not fraction. 50 does not work yet.
    %SPECTRUM_OVERLAP_PERCENT = 50;
    
    % NOTE: More samples per spectrum is faster.
    %N_SAMPLES_PER_SPECTRUM = 2048 / 4;   % TEST
    N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.

    
    TsCa  = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz);

%     t = tic;
%     nSamplesTotal = 0;

    h = irf_panel(panelTag);    
    
    %====================
    % Calculate spectras
    %====================
    SpecrecCa = {};
    ssCenterEpochUnixArray = [];
    for i = 1:numel(TsCa)    
        Ts = TsCa{i};
        
%         nSamples      = numel(Ts.time.epoch);
%         fprintf('nSamples = %g\n', nSamples)
%         nSamplesTotal = nSamplesTotal + nSamples;
        
        Specrec = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz, SPECTRUM_OVERLAP_PERCENT);
        SpecrecCa{end+1} = Specrec;
        
        % IMPLEMENTATION NOTE: Later needs the snapshot centers in the same time system as Specrec.t (epoch Unix).
        ssCenterEpochUnixArray(end+1) = (Ts.time.start.epochUnix + Ts.time.stop.epochUnix)/2;
    end
    ssPeriodSec = min(diff(ssCenterEpochUnixArray));
    
    %==================================================================================================================
    % Set the display locations of individual spectras (override defaults). Separately stretch out the collection of
    % spectras that stems from every snapshot.
    % IMPLEMENTATION NOTE: This can not be done in the first loop in order to derive the (minimum) snapshot time
    % distance.
    %==================================================================================================================
    for i = 1:numel(TsCa)
        ssLengthSec = TsCa{i}.time.stop.epochUnix - TsCa{i}.time.start.epochUnix;
        
        % Stretch out spectra (for given snapshot) in time to be ALMOST adjacent between snapshots.
        % NOTE: Specrec.dt is not set by irf_powerfft so there is no default value that can be scaled up.
        % NOTE: Uses original spectrum positions and re-positions them relative to snapshot center.
        scaleFactor = ssPeriodSec/ssLengthSec * SNAPSHOT_WIDTH_FRACTION;
        SpecrecCa{i}.t  = ssCenterEpochUnixArray(i) + (SpecrecCa{i}.t - ssCenterEpochUnixArray(i)) * scaleFactor;
        SpecrecCa{i}.dt = ones(size(SpecrecCa{i}.t)) * min(diff(SpecrecCa{i}.t)) * 0.5;
    end
    
    Specrec = merge_specrec(SpecrecCa);
    
    Specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
    
    irf_spectrogram(h, Specrec);   % Replaces irf_plot
    
    set(h, 'yscale','log')

    irf_legend(h, tlLegend, [0.02 0.98], 'color', 'k')
    irf_legend(h, trLegend, [0.98 0.98])

%     tSec = toc(t);
%     fprintf('tSec/nSamples = %g [s/sample]\n', tSec/nSamplesTotal)
end



% Wrapper for converting data.
function h = time_series_panel2(panelTag, zvEpoch, zvDataList, samplingFreqHz, tlLegend, trLegend)
    
    % Convert from N samples/record --> 1 samples/record (row).
    nSps    = size(zvDataList{1}, 2);   % SPS = Samples Per Snapshot
    zvEpoch = bicas.proc_utils.convert_N_to_1_SPR_Epoch(zvEpoch, nSps, ones(size(zvEpoch))*samplingFreqHz);
    for i = 1:numel(zvDataList)
        zvDataList{i} = bicas.proc_utils.convert_N_to_1_SPR_redistribute(zvDataList{i});
    end
    
    % NOTE: Effectively serves as an assertion on zv sizes.
    Ts = irf.ts_scalar(zvEpoch, [zvDataList{:}]);
    
    h = time_series_panel(panelTag, Ts, tlLegend, trLegend);
end



% ARGUMENTS
% =========
% tlLegend : Top-left  (TL) legend.
% trLegend : Top-right (TR) legend. Cell array of strings, one per scalar time series.
function h = time_series_panel(panelTag, Ts, tlLegend, trLegend)
    h = irf_panel(panelTag);
    irf_plot(h, Ts)
    ylabel(h, '[V]')
    irf_legend(h, tlLegend, [0.02 0.98], 'color', 'k')
    irf_legend(h, trLegend, [0.98 0.98])
end



% Convert zVar-like variables for snapshots to cell array of TSeries.
%
% ARGUMENTS
% =========
% zvEpoch : Nx1 array.
% zvData  : NxM array. (iRecord, iSampleWithinSnapshot). 1 record=1 snapshot.
% TsCa    : (iSnapshot) 1D cell array of TSeries.
%           IMPLEMENTATION NOTE: Can not(?) be struct array since MATLAB confuses indexing a TSeries array (with
%           brackets) with some special TSeries functionality for calling its code with brackets (calling TSeries'
%           method "subsref").
%
function TsCa = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz)
    % IMPLEMENTATION NOTE: Function is written to some day be easily extended to be used for use with TDS's
    % length-varying snapshots.
    %
    % NOTE: No special treatment of snapshots with only NaN.
    
    assert(isscalar(samplingFreqHz))
    EJ_library.assert.size(zvData, [NaN, NaN])
    assert(size(zvEpoch, 1) == size(zvData, 1))   % Same number of records
    bicas.proc_utils.assert_Epoch(zvEpoch)
    
    nRecords = size(zvData, 1);
    nSps     = size(zvData, 2);
    assert(nSps >= 2)
    
    epochRelArray = int64([0:(nSps-1)] * 1/samplingFreqHz * 1e9);    % Relative timestamps inside CDF record/snapshot.
    
    TsCa = {};
    for i = 1:nRecords
        epochRecord = zvEpoch(i) + epochRelArray;
        TsCa{i}  = irf.ts_scalar(epochRecord, zvData(i, :));
    end
    
end



% Merge multiple instances of "specrec" structs as returned by irf_powerfft 
% NOTE: Optionally added fields must be added after merging.
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% SpecrecCa : Cell array of "specrec" as returned by irf_powerfft, but with .dt (column array) added to it.
%             Unsure if irf_powerfft can return more cases than can be handled here.
%             NOTE: Includes dt (column array of scalars).
%             NOTE: Assumes that all specrec use the same frequencies.
%             IMPLEMENTATION NOTE: Uses cell array instead of struct array to be able to handle (and ignore) specrec =
%             [] which can be returned by irf_powerfft.
% Specrec   : Struct array that can be used by irf_spectrogram.
%
function Specrec = merge_specrec(SpecrecCa)
    % PROPOSAL: Assertion for frequencies.
    
    Specrec.f = [];
    Specrec.p = {[]};
    Specrec.t = [];
    Specrec.dt = [];
    
    for i = 1:numel(SpecrecCa)
        
        S = SpecrecCa{i};
        if ~isempty(S)
            EJ_library.assert.struct(S, {'f', 'p', 't', 'dt'}, {});
            assert(iscolumn(S.dt))
            assert(numel(S.dt) == numel(S.t), 'Badly formatted SpecrecCa{%i}.', i)
            
            Specrec.f    = S.f;
            Specrec.p{1} = [Specrec.p{1}; S.p{1}];
            Specrec.t    = [Specrec.t;    S.t(:)];
            Specrec.dt   = [Specrec.dt;   S.dt(:)];
        end
    end
end
