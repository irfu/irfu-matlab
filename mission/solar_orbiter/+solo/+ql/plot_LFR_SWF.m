%
% Quicklook for the content of one BIAS LFR SWF dataset (CDF file), i.e. DATASET_ID = SOLO_L2_RPW-LFR-SURV-SWF-E
%
%
% NOTE: Only capable (default) of only showing either DC diffs or AC diffs. There are hardcoded settings for permitting
% or forcing both.
% NOTE: Uses bicas.proc_utils.* code.
% NOTE: Does not yet support spectrogram overlap.
% NOTE: Time series panels interpolate between snapshots.
% NOTE: Color scale is log, therefore negative values (probably).
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
    % PROPOSAL: Argument for time interval should use some more irfu-matlab-way of specifying a time interval.
    % PROPOSAL: Some clear way of distinguishing AC & DC visually?
    %
    % TODO?: Remove interpolation between time series snapshots?
    % TODO: 50% overlap PSD
    
    % YK 2020-04-16: Officially only either DC or AC diffs.
    % NOTE: solo_L2_rpw-lfr-surv-cwf-e-cdag_20200228_V01.cdf contains both DC & AC diffs.
    ALWAYS_SIMULTANEOUS_DC_AC_DIFFS_PLOTS = 0;   % DEFAULT 0. Useful for debugging (runs through all code).
    PERMIT_SIMULTANEOUS_DC_AC_DIFFS       = 1;   % DEFAULT 0.
    ENABLE_SPECTROGRAMS                   = 1;   % DEFAULT 1.
    
    % Info associated with LFR sampling rates (F0-F3 is LFR's terminology).
    % NOTE: LFR SWF only uses F0-F2 (not F3).
    % RATIONALE: Useful to be able to submit (to a function) all info associated with one sampling rate at once.
    F0.str    = 'F0';
    F1.str    = 'F1';
    F2.str    = 'F2';
    F0.freqHz = 24576;
    F1.freqHz =  4096;
    F2.freqHz =   256;

    
    
    % Interpret argument for optionally restricting covered time range.
    % TEMPORARY? Useful for debugging when plotting large files (slow plotting).
    if nargin == 1
        DATAOBJ_TIME_INTERVAL_ARGS = {};
    elseif nargin == 2
        assert(iscell(timeIntervUtc))
        DATAOBJ_TIME_INTERVAL_ARGS = {'tint', [iso2epoch(timeIntervUtc{1}), iso2epoch(timeIntervUtc{2})]};
    end
    
    

    D = dataobj(filePath, DATAOBJ_TIME_INTERVAL_ARGS{:});
    
    epoch    = D.data.Epoch.data;
    F_SAMPLE = D.data.F_SAMPLE.data;
    vDc1     = get_CDF_zv_data(D, 'V',   1);
    vDc12    = get_CDF_zv_data(D, 'E',   1);
    vDc23    = get_CDF_zv_data(D, 'E',   3);
    vAc12    = get_CDF_zv_data(D, 'EAC', 1);
    vAc23    = get_CDF_zv_data(D, 'EAC', 3);
    
    % B = Boolean/Logical (true/false for every index value).
    F0.bRecords = (F_SAMPLE == F0.freqHz);
    F1.bRecords = (F_SAMPLE == F1.freqHz);
    F2.bRecords = (F_SAMPLE == F2.freqHz);
    assert(all(F0.bRecords | F1.bRecords | F2.bRecords))
    

    
    %=================================================================
    % Determine whether DC diffs, AC diffs, or both should be plotted
    %=================================================================
    hasDcDiffs = any(~isnan(vDc12(:))) || any(~isnan(vDc23(:)));
    hasAcDiffs = any(~isnan(vAc12(:))) || any(~isnan(vAc23(:)));
    if ALWAYS_SIMULTANEOUS_DC_AC_DIFFS_PLOTS
        displayDcDiffs = 1;
        displayAcDiffs = 1;
    else
        % ASSERTIONS
        if      hasDcDiffs && hasAcDiffs && ~PERMIT_SIMULTANEOUS_DC_AC_DIFFS
            error('Dataset (CDF file) contains both DC diff and AC diff data. Can not handle this case.')
        elseif ~hasDcDiffs && ~hasAcDiffs
            error('Dataset (CDF file) contains neither DC diff nor AC diff data. Can not handle this case.')
        end
        
        displayDcDiffs = hasDcDiffs;
        displayAcDiffs = hasAcDiffs;
    end
    

    
    pcfcList = {};    % PCFC = Panel Creation Function Call
    if ENABLE_SPECTROGRAMS
        %=================
        % F0 spectrograms
        %=================
        pcfcList{end+1}     = @() (spectrogram_panel2( 'V1 DC', epoch, vDc1,  F0, 'V1\_DC'));
        if displayDcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 DC', epoch, vDc12, F0, 'V12\_DC'));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 DC', epoch, vDc23, F0, 'V13\_DC'));
        end
        if displayAcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 AC', epoch, vAc12, F0, 'V12\_AC'));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 AC', epoch, vAc23, F0, 'V23\_AC'));
        end
        %=================
        % F1 spectrograms
        %=================
        pcfcList{end+1} =     @() (spectrogram_panel2( 'V1 DC', epoch, vDc1,  F1, 'V1\_DC'));
        if displayDcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 DC', epoch, vDc12, F1, 'V12\_DC'));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 DC', epoch, vDc23, F1, 'V23\_DC'));
        end
        if displayAcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 AC', epoch, vAc12, F1, 'V12\_AC'));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 AC', epoch, vAc23, F1, 'V23\_AC'));
        end
        %=================
        % F2 spectrograms
        %=================
        pcfcList{end+1}     = @() (spectrogram_panel2( 'V1 DC', epoch, vDc1,  F2, 'V1\_DC'));
        if displayDcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 DC', epoch, vDc12, F2, 'V12\_DC'));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 DC', epoch, vDc23, F2, 'V23\_DC'));
        end
        if displayAcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 AC', epoch, vAc12, F2, 'V12\_AC'));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 AC', epoch, vAc23, F2, 'V23\_AC'));
        end
    end
    %==========================================================================================
    % F0-F2 time series 
    % IMPLEMENTATION NOTE: Panel tags have to be unique, or otherwise the axes will be reused.
    %==========================================================================================
    if displayDcDiffs
        %======================
        % DC single + DC diffs
        %======================
        SIGNALS_LEGEND_DC = EJ_library.graph.escape_str({'V1_DC','V12_DC','V23_DC'});
        tempFuncPtr = @(Fx) (@() (time_series_panel2('V1,V12,V23 DC', epoch, {vDc1, vDc12, vDc23}, Fx, SIGNALS_LEGEND_DC)));
        
        pcfcList{end+1} = tempFuncPtr(F0);
        pcfcList{end+1} = tempFuncPtr(F1);
        pcfcList{end+1} = tempFuncPtr(F2);

    end
    if ~displayDcDiffs && displayAcDiffs
        %======================
        % DC single + AC diffs
        %======================
        SIGNALS_LEGEND_DC_AC = EJ_library.graph.escape_str({'V1_DC','V12_AC','V23_AC'});
        tempFuncPtr = @(Fx) (@() (time_series_panel2('V1,V12,V23 DC/AC', epoch, {vDc1, vAc12, vAc23}, Fx, SIGNALS_LEGEND_DC_AC)));

        pcfcList{end+1} = tempFuncPtr(F0);
        pcfcList{end+1} = tempFuncPtr(F1);
        pcfcList{end+1} = tempFuncPtr(F2);

    end    
    if displayDcDiffs && displayAcDiffs
        %======================
        % AC diffs (no single)
        %======================
        SIGNALS_LEGEND_AC = EJ_library.graph.escape_str({'V12_AC','V23_AC'});
        tempFuncPtr = @(Fx) (@() (time_series_panel2('V12,V23 AC', epoch, {vAc12, vAc23}, Fx, SIGNALS_LEGEND_AC)));

        pcfcList{end+1} = tempFuncPtr(F0);
        pcfcList{end+1} = tempFuncPtr(F1);
        pcfcList{end+1} = tempFuncPtr(F2);

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
    irf_zoom(hAxesArray, 'x', irf.tint(epoch(1), epoch(end)))    % For aligning the content of the MATLAB axes.    
end



% Convenient wrapper around spectrum_panel.
% Converts from zVar-like variables to what is actually used for plotting.
function h = spectrogram_panel2(panelTagSignalsStr, zvEpoch, zvData, SamplingRateInfo, trLegend)
    h = spectrogram_panel(...
        sprintf('%s %s spectrogram', panelTagSignalsStr, SamplingRateInfo.str), ...
        zvEpoch(SamplingRateInfo.bRecords, :), ...
        zvData(SamplingRateInfo.bRecords, :), ...
        SamplingRateInfo.freqHz, ...
        SamplingRateInfo.str, ...
        trLegend);
end



% ARGUMENTS
% =========
% TsCa     : Cell array of TSeries. All TSeries separately describe different time segments of one single scalar time
%            series. In practice, each TSeries should represent one snapshot. This speeds up the spectrograms by a lot.
% tlLegend : Top-left  (TL) legend string.
% trLegend : Top-right (TR) legend string.
%
function h = spectrogram_panel(panelTag, zvEpoch, zvData, samplingFreqHz, tlLegend, trLegend)
    % NOTE: Multiple-row labels causes trouble for the time series ylabels.
    % IMPLEMENTATION NOTE: Implemented to potentially be modified to handle TDS snapshots that vary in length.

    % Fraction of the (minimum) time distance between snapshots (centers) that will be used for displaying the spectra.
    % Value 1 : Spectras are adjacent between snapshot (for minimum snapshot distance).
    SNAPSHOT_WIDTH_FRACTION  = 0.90;
    %SPECTRUM_OVERLAP_PERCENT = 0;    % Percent, not fraction.
    SPECTRUM_OVERLAP_PERCENT = 50;    % Percent, not fraction.
    
    % NOTE: More samples per spectrum is faster (sic!).
    %N_SAMPLES_PER_SPECTRUM = 2048 / 2;   % TEST
    N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.
    

    
    TsCa  = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz);

    h = irf_panel(panelTag);    
    
    %====================
    % Calculate spectras
    %====================
    SpecrecCa = {};
    ssCenterEpochUnixArray = [];
    for i = 1:numel(TsCa)    
        Ts = TsCa{i};
        
        SpecrecCa{end+1} = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz, SPECTRUM_OVERLAP_PERCENT);
        
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

end



% Convenient wrapper around time_series_panel.
% Converts from zVar-like variables (N samples/record; all records) to what is actually used for plotting.
function h = time_series_panel2(panelTagSignalsStr, zvEpoch, zvDataList, SamplingRateInfo, trLegend)
    
    nSps     = size(zvDataList{1}, 2);   % SPS = Samples Per Snapshot
    
    zvEpoch  = zvEpoch(SamplingRateInfo.bRecords);
    nRecords = size(zvEpoch, 1);   % NOTE: After selecting records.
    zvEpoch  = bicas.proc_utils.convert_N_to_1_SPR_Epoch(zvEpoch, nSps, ones(nRecords, 1)*SamplingRateInfo.freqHz);
    
    for i = 1:numel(zvDataList)
        zvData = zvDataList{i}(SamplingRateInfo.bRecords, :);
        zvDataList{i} = bicas.proc_utils.convert_N_to_1_SPR_redistribute(zvData);
    end
    
    % NOTE: Effectively serves as an assertion on zv sizes.
    Ts = irf.ts_scalar(zvEpoch, [zvDataList{:}]);
    
    % PLOT
    h = time_series_panel(...
        sprintf('%s %s time series', panelTagSignalsStr, SamplingRateInfo.str), ...
        Ts, SamplingRateInfo.str, trLegend);
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
    bicas.proc_utils.assert_zv_Epoch(zvEpoch)
    
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



% Merge multiple instances of "specrec" structs as returned by irf_powerfft.
% NOTE: Optionally added fields must be added after merging.
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% SpecrecCa : Cell array of "specrec" as returned by irf_powerfft, but with .dt (column array) added to it.
%             Unsure if irf_powerfft can return more cases than can be handled here.
%             NOTE: Requires dt (column array of scalars).
%             NOTE: Assumes that all specrec use the same frequencies.
%             IMPLEMENTATION NOTE: Uses cell array instead of struct array to be able to handle (and ignore) the case
%             specrec = [] which can be returned by irf_powerfft.
% Specrec   : Struct array that can be used by irf_spectrogram.
%
function Specrec = merge_specrec(SpecrecCa)
    % PROPOSAL: Assertion for frequencies.
    
    Specrec.f  = [];
    Specrec.p  = {[]};   % NOTE: Must 1x1 cell array. The array INSIDE the cell array is added to.
    Specrec.t  = [];
    Specrec.dt = [];
    
    for i = 1:numel(SpecrecCa)
        
        S = SpecrecCa{i};
        if ~isempty(S)
            EJ_library.assert.struct(S, {'f', 'p', 't', 'dt'}, {});
            assert(iscolumn(S.dt))
            assert(numel(S.dt) == numel(S.t), 'Badly formatted SpecrecCa{%i}.', i)
            
            Specrec.f    = S.f;                       % NOTE: Not adding to array, but setting it in its entirety.
            Specrec.p{1} = [Specrec.p{1}; S.p{1}];    % NOTE: Add to array inside cell array.
            Specrec.t    = [Specrec.t;    S.t(:)];    % NOTE: Has to be column vector.
            Specrec.dt   = [Specrec.dt;   S.dt(:)];   % NOTE: Has to be column vector.
        end
    end
end



function data = get_CDF_zv_data(D, zvName, i3)
    fillValue = getfillval(D, zvName);
    data = D.data.(zvName).data(:, :, i3);
    data = changem(data, NaN, fillValue);
end
