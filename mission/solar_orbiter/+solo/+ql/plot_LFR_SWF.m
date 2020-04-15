%
% Quicklook for the content of one BIAS LFR SWF dataset (CDF file), i.e. DATASET_ID = SOLO_L2_RPW-LFR-SURV-SWF-E
%
% NOTE: Uses bicas.proc_utils.* code.
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
    % BUG: Snapshot spectrograms consist of very narrow stripes.
    % BUG: Snapshot spectrogram for F2 misses some snapshots.
    % TODO/BUG: Functioning ~color map
    % TODO: DC/AC: Detection, label.
    % TODO: Correct sampling frequency for irf_powerfft.
    % TODO: Legend connecting line color to V1, V12, V13.
    % TODO: Read fill values from file.
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
    % TODO: Add AC spectrograms+AC time series.
    
    
    FILL_VALUE = single(-1e31);
    N_SAMPLES_PER_SNAPSHOT = 2048;    % Number of Samples Per (LFR) Snapshot (SPS).
    
    % Interpret argument for optionally restricting covered time range.
    % TEMPORARY? Useful for debugging when plotting large files (slow plotting).
    if nargin == 1
        DATAOBJ_TIME_INTERVAL_ARGS = {};
    elseif nargin == 2
        assert(iscell(timeIntervUtc))
        DATAOBJ_TIME_INTERVAL_ARGS = {'tint', [iso2epoch(timeIntervUtc{1}), iso2epoch(timeIntervUtc{2})]};
    end
    
    

    D = dataobj(filePath, DATAOBJ_TIME_INTERVAL_ARGS{:});
    Epoch = D.data.Epoch.data;
    %fprintf('nRecords = %i\n', numel(Epoch));

    VDC1  = D.data.V.data(:,:,1);
    VDC12 = D.data.E.data(:,:,1);
    VDC23 = D.data.E.data(:,:,3);
    VAC12 = D.data.EAC.data(:,:,1);
    VAC23 = D.data.EAC.data(:,:,3);
    F_SAMPLE = D.data.F_SAMPLE.data;
    
    VDC1  = changem(VDC1,  NaN, FILL_VALUE);
    VDC12 = changem(VDC12, NaN, FILL_VALUE);
    VDC23 = changem(VDC23, NaN, FILL_VALUE);
    VAC12 = changem(VAC12, NaN, FILL_VALUE);
    VAC23 = changem(VAC23, NaN, FILL_VALUE);

    
    % LFR sampling frequencies (F0-F3 is LFR's terminology).
    F0Hz = 24576;
    F1Hz =  4096;
    F2Hz =   256;
    %F3Hz =    16;
    % LFR SWF only uses F0-F2.

    % B = Boolean/Logical (true/false for every index value).
    bF0 = (F_SAMPLE == F0Hz);
    bF1 = (F_SAMPLE == F1Hz);
    bF2 = (F_SAMPLE == F2Hz);
    
    % DEBUG
%     utc1 = erikpgjohansson.utils.CDF_tt2000_to_UTC_str(Epoch(1));
%     utc2 = erikpgjohansson.utils.CDF_tt2000_to_UTC_str(Epoch(end));
%     fprintf('%s -- %s\n', utc1, utc2);
    
    
    % NOTE: Using bicas.* code.
    EpochF0  = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF0), N_SAMPLES_PER_SNAPSHOT, F_SAMPLE(bF0));
    Vdc1F0  = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC1( bF0, :));
    Vdc12F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC12(bF0, :));
    Vdc23F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC23(bF0, :));
    Vac12F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC12(bF0, :));
    Vac23F0 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VAC23(bF0, :));
    
    EpochF1 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF1), N_SAMPLES_PER_SNAPSHOT, F_SAMPLE(bF1));
    Vdc1F1  = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC1( bF1, :));
    Vdc12F1 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC12(bF1, :));
    Vdc23F1 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC23(bF1, :));
    
    EpochF2 = bicas.proc_utils.convert_N_to_1_SPR_Epoch(Epoch(bF2), N_SAMPLES_PER_SNAPSHOT, F_SAMPLE(bF2));
    Vdc1F2  = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC1( bF2, :));
    Vdc12F2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC12(bF2, :));
    Vdc23F2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(VDC23(bF2, :));

    
    
    % Create cell arrays of TSeries for each snapshot (on each channel).
    TsVdc1F0_array  = snapshot_per_record_2_TSeries(Epoch(bF0), VDC1( bF0, :), F0Hz);
    TsVdc12F0_array = snapshot_per_record_2_TSeries(Epoch(bF0), VDC12(bF0, :), F0Hz);
    TsVdc23F0_array = snapshot_per_record_2_TSeries(Epoch(bF0), VDC23(bF0, :), F0Hz);
    TsVac12F0_array = snapshot_per_record_2_TSeries(Epoch(bF0), VAC12(bF0, :), F0Hz);
    TsVac23F0_array = snapshot_per_record_2_TSeries(Epoch(bF0), VAC23(bF0, :), F0Hz);
    
    TsVdc1F1_array  = snapshot_per_record_2_TSeries(Epoch(bF1), VDC1( bF1, :), F1Hz);
    TsVdc12F1_array = snapshot_per_record_2_TSeries(Epoch(bF1), VDC12(bF1, :), F1Hz);
    TsVdc23F1_array = snapshot_per_record_2_TSeries(Epoch(bF1), VDC23(bF1, :), F1Hz);
    TsVac12F1_array = snapshot_per_record_2_TSeries(Epoch(bF1), VAC12(bF1, :), F1Hz);
    TsVac23F1_array = snapshot_per_record_2_TSeries(Epoch(bF1), VAC23(bF1, :), F1Hz);
    
    TsVdc1F2_array  = snapshot_per_record_2_TSeries(Epoch(bF2), VDC1( bF2, :), F2Hz);
    TsVdc12F2_array = snapshot_per_record_2_TSeries(Epoch(bF2), VDC12(bF2, :), F2Hz);
    TsVdc23F2_array = snapshot_per_record_2_TSeries(Epoch(bF2), VDC23(bF2, :), F2Hz);
    TsVac12F2_array = snapshot_per_record_2_TSeries(Epoch(bF2), VAC12(bF2, :), F2Hz);
    TsVac23F2_array = snapshot_per_record_2_TSeries(Epoch(bF2), VAC23(bF2, :), F2Hz);
    
    % Create TSeries representing 3 scalar time series each.
    TsVdcF0 = irf.ts_scalar(EpochF0, [Vdc1F0, Vdc12F0, Vdc23F0]);
    TsVdcF1 = irf.ts_scalar(EpochF1, [Vdc1F1, Vdc12F1, Vdc23F1]);
    TsVdcF2 = irf.ts_scalar(EpochF2, [Vdc1F2, Vdc12F2, Vdc23F2]);


    % TsV12F0_DC
    % TsVDC12F0
    % TsVdc12F0

    irf_plot(5+3+3 +3, 'newfigure');
    %irf_plot(3 + 3,'newfigure');
    
    hAxesArray = [];
    hAxesArray(end+1) = spectrum_panel( 'V1 DC F0 spectrogram', TsVdc1F0_array,  F0Hz, 'F0', 'V1\_DC');
    hAxesArray(end+1) = spectrum_panel('V12 DC F0 spectrogram', TsVdc12F0_array, F0Hz, 'F0', 'V12\_DC');
    hAxesArray(end+1) = spectrum_panel('V23 DC F0 spectrogram', TsVdc23F0_array, F0Hz, 'F0', 'V13\_DC');
    hAxesArray(end+1) = spectrum_panel('V12 AC F0 spectrogram', TsVac12F0_array, F0Hz, 'F0', 'V12\_AC');
    hAxesArray(end+1) = spectrum_panel('V23 AC F0 spectrogram', TsVac23F0_array, F0Hz, 'F0', 'V13\_AC');
    
    hAxesArray(end+1) = spectrum_panel( 'V1 F1 spectrogram', TsVdc1F1_array,  F1Hz, 'F1', 'V1\_DC');
    hAxesArray(end+1) = spectrum_panel('V12 F1 spectrogram', TsVdc12F1_array, F1Hz, 'F1', 'V12\_DC');
    hAxesArray(end+1) = spectrum_panel('V23 F1 spectrogram', TsVdc23F1_array, F1Hz, 'F1', 'V13\_DC');
    
    hAxesArray(end+1) = spectrum_panel( 'V1 F2 spectrogram', TsVdc1F2_array,  F2Hz, 'F2', 'V1\_DC');
    hAxesArray(end+1) = spectrum_panel('V12 F2 spectrogram', TsVdc12F2_array, F2Hz, 'F2', 'V12\_DC');
    hAxesArray(end+1) = spectrum_panel('V23 F2 spectrogram', TsVdc23F2_array, F2Hz, 'F2', 'V13\_DC');

    SIGNALS_LEGEND = EJ_library.graph.escape_str({'V1_DC','V12_DC','V23_DC'});
    hAxesArray(end+1) = time_series_panel('V1,V12,V23 F0 time series', TsVdcF0, 'F0', SIGNALS_LEGEND);
    hAxesArray(end+1) = time_series_panel('V1,V12,V23 F1 time series', TsVdcF1, 'F1', SIGNALS_LEGEND);
    hAxesArray(end+1) = time_series_panel('V1,V12,V23 F2 time series', TsVdcF2, 'F2', SIGNALS_LEGEND);
    
%     hAxesArray = [];
%     nPanels = numel(panelCreationList);
%     for i = 1:nPanels
%         hAxesArray(end+1) = panelCreationList;
%     end

    solo.ql.set_std_title('LFR SWF L2', filePath, hAxesArray(1))
    
    irf_plot_axis_align(hAxesArray)                        % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(hAxesArray, 'x', irf.tint(TsVdcF0.time))    % For aligning the content of the MATLAB axes.    
end



% ARGUMENTS
% =========
% TsCa     : Cell array of TSeries. All TSeries separately describe different time segments of one single scalar time series.
% tlLegend : Top-left  (TL) legend string.
% trLegend : Top-right (TR) legend string.
function h = spectrum_panel(panelTag, TsCa, samplingFreqHz, tlLegend, trLegend)
    % NOTE: Two-rows labels causes trouble for the time series ylabels.
    % IMPLEMENTATION NOTE: Implemented to be able to handle TDS snapshots that vary in length (in theory; untested).
    
    % Fraction of the (minimum) time distance between snapshots (centers) that will be used for displaying the spectra.
    % =1 : Spectras are adjacent between snapshot (for minimum snapshot distance).
    SNAPSHOT_WIDTH_FRACTION = 0.8;
    
    % NOTE: More samples per spectrum is faster.
    %N_SAMPLES_PER_LFR_SNAPSHOT = 2048;    % Number of Samples Per (LFR) Snapshot (SPS).
    %N_SAMPLES_PER_SPECTRUM = N_SAMPLES_PER_LFR_SNAPSHOT / 4;   % TEST
    N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.

    
    %TsCa = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz);
    
    
    t = tic;
    nSamplesTotal = 0;

    h = irf_panel(panelTag);
    
    
    %====================
    % Calculate spectras
    %====================
    SpecrecCa = {};
    snapshotCenterSecArray = [];
    for i = 1:numel(TsCa)    
        Ts = TsCa{i};
        
        nSamples      = numel(Ts.time.epoch);
%         fprintf('nSamples = %g\n', nSamples)
        nSamplesTotal = nSamplesTotal + nSamples;
        
        Specrec = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz);        
        SpecrecCa{end+1} = Specrec;        
        
        % IMPLEMENTATION NOTE: Later needs the snapshot centers in the same time system as Specrec.t. Therefore using
        % Specrec.t to derive snapshot centers.
        snapshotCenterSecArray(end+1) = mean(Specrec.t);
    end
    snapshotPeriodSec = min(diff(snapshotCenterSecArray));
    
    %==================================================================================================================
    % Set the display locations of individual spectras (override defaults). Separately stretch out the collection of
    % spectras that stems from every snapshot.
    %==================================================================================================================
    for i = 1:numel(TsCa)
        nSpectra = numel(SpecrecCa{i}.t);
        
        % Make spectra for snapshot stretch out in time to be almost adjacent 
        relativeArraySec = ([0:(nSpectra-1)]' + 0.5 - nSpectra/2) * snapshotPeriodSec/nSpectra * SNAPSHOT_WIDTH_FRACTION;
        SpecrecCa{i}.t  = snapshotCenterSecArray(i) + relativeArraySec;
        SpecrecCa{i}.dt = ones(nSpectra, 1) * snapshotPeriodSec/nSpectra * 0.5;
    end
    
    Specrec = merge_specrec(SpecrecCa);
    
    Specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
    %Specrec.dt      = nSamples / samplingFreqHz * 0.4 * 1e3;    % TEST
    
    irf_spectrogram(h, Specrec);   % Replaces irf_plot
    
    set(h, 'yscale','log')

    irf_legend(h, {tlLegend}, [0.02 0.98], 'color', 'k')
    irf_legend(h, trLegend ,  [0.98 0.98])

    tSec = toc(t);
    fprintf('tSec/nSamples = %g [s/sample]\n', tSec/nSamplesTotal)
end



% Convert zVar-like data for snapshots to cell array of TSeries.
%
% ARGUMENTS
% =========
% zvEpoch : Nx1 array.
% zvData  : NxM array.
% TsCa    : 1D cell array of TSeries.
%           IMPLEMENTATION NOTE: Can not(?) be struct array since MATLAB confuses indexing a TSeries array (with
%           brackets) with some special TSeries functionality for calling it with brackets (calling TSeries' method
%           "subsref").
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



% ARGUMENTS
% =========
% tlLegend : Top-left  (TL) legend. String.
% trLegend : Top-right (TR) legend. Cell array of strings, one per scalar time series.
function h = time_series_panel(panelTag, Ts, tlLegend, trLegend)
    h = irf_panel(panelTag);
    irf_plot(h, Ts)
    ylabel(h, '[V]')
    irf_legend(h, {tlLegend}, [0.02 0.98], 'color', 'k')
    irf_legend(h, trLegend,   [0.98 0.98])
end



% Merge multiple instances of "specrec" structs as returned by irf_powerfft
% NOTE: Optionally added fields must be added after merging.
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% SpecrecCa : Cell array of "specrec". Unsure if irf_powerfft supports more cases
%             than can be handled here.
%             NOTE: Includes dt (array of scalars).
%             NOTE: Assumes that all use the same frequencies.
%             IMPLEMENTATION NOTE: Uses cell array instead of struct array to be able to handle (and ignore) specrec =
%             [] as can be returned by irf_powerfft.
% Specrec   : Struct array that irf_spectrogram can use.
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
            
            Specrec.f    = S.f;
            Specrec.p{1} = [Specrec.p{1}; S.p{1}];
            Specrec.t    = [Specrec.t;    S.t(:)];
            Specrec.dt   = [Specrec.dt;   S.dt(:)];
        end
    end
end
