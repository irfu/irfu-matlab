%
% Quicklook for the content of one BIAS LFR SWF dataset (CDF file), i.e.
% DATASET_ID = SOLO_L2_RPW-LFR-SURV-SWF-E.
%
%
% NOTE: Only capable (default) of only showing either DC diffs or AC diffs.
% There are hard-coded settings for permitting or forcing both.
% NOTE: Uses bicas.proc_utils.* code.
% NOTE: Does not yet support spectrogram overlap.
% NOTE: Time series panels interpolate between snapshots.
% NOTE: Color scale is log, therefore negative values (probably).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-28.
%
function hAxesArray = plot_LFR_SWF(filePath)
    % SPEED
    % =====
    % Execution time
    % With    irf_spectrogram : 26.592199 s
    % Without irf_spectrogram : 36.683050 s
    % With    irf_spectrogram, downsample_Specrec(nBins=2) for every
    % snapshot                : 31.648818 s
    % ==> irf_spectrogram is only a part of the problem. It is not enough to use
    % downsample_Specrec.
    % --
    % N_SAMPLES_PER_SPECTRUM = 128 ==> 29.063167 s
    % N_SAMPLES_PER_SPECTRUM = 512 ==> 17.913298 s
    % --
    % NOTE: 24 spectrums per snapshot (N_SAMPLES_PER_SPECTRUM = 128)
    % ==> not much use downsampling further).
    %
    %
    % BOGIQ:
    % =====
    % POLICY: Combined BOGIQ for all quicklook plot code: See plot_LFR_CWF.
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
    % PROPOSAL: Some clear way of distinguishing AC & DC visually?
    %
    % TODO?: Remove interpolation (straight lines) between time series snapshots?
    %
    % ~BUG: Algorithm for calculating enlargement of snapshot spectra not always appropriate.
    %   Ex: solo_L2_rpw-lfr-surv-swf-e-cdag_20200228_V01.cdf
    %   PROPOSAL: Enlarge each snapshot spectrum until it reaches neighbour (roughly).
    %       Use min(timeToSnapshotBefore, timeToSnapshotafter) as max radius.
    %
    % Old TODO: YK's "fixup" 2020-10-13
    %   hswf = solo.ql.plot_LFR_SWF([RPWPATH LRFFILE]);%% fixup
    %   load cmap
    %   colormap(cmap)
    %   set(hswf(1:3),'YTick',[0.1 1 10])
    %   caxis(hswf(1),[-11 -7])
    %   caxis(hswf(2),[-13 -10])
    %   caxis(hswf(3),[-13 -10])
    %   caxis(hswf(4),[-9 -5])
    %   set(hswf(4:6),'YTick',[0.1 1 10])
    %   caxis(hswf(5),[-11 -8])
    %   caxis(hswf(6),[-11 -8])
    %   set(hswf(7:9),'YTick',[0.1 1 10 100])
    %   caxis(hswf(7),[-8 -5])
    %   caxis(hswf(8),[-10 -7])
    %   caxis(hswf(9),[-10 -7])
    %   2020-10-13: Should be implemented.
    
    

    % YK 2020-04-16: Officially only either DC or AC diffs.
    % NOTE: solo_L2_rpw-lfr-surv-cwf-e-cdag_20200228_V01.cdf contains both DC &
    % AC diffs.
    ALWAYS_SIMULTANEOUS_DC_AC_DIFFS_PLOTS = 0;   % DEFAULT 0. Useful for debugging (runs through all code).
    PERMIT_SIMULTANEOUS_DC_AC_DIFFS       = 1;   % DEFAULT 0.
    ENABLE_SPECTROGRAMS                   = 1;   % DEFAULT 1.
    
    D = dataobj(filePath);
    
    epoch    = D.data.Epoch.data;
    F_SAMPLE = get_CDF_zv_data(D, 'SAMPLING_RATE', 1);
    vDc1     = get_CDF_zv_data(D, 'VDC', 1);
    vDc12    = get_CDF_zv_data(D, 'EDC', 1);
    vDc23    = get_CDF_zv_data(D, 'EDC', 3);
    vAc12    = get_CDF_zv_data(D, 'EAC', 1);
    vAc23    = get_CDF_zv_data(D, 'EAC', 3);
    clear D
    
    if 0
        % DEBUG: Limit records
        I1 = 1;
        I2 = 100;
        
        epoch = epoch(I1:I2);
        F_SAMPLE = F_SAMPLE(I1:I2);
        vDc1  = vDc1( I1:I2, :);
        vDc12 = vDc12(I1:I2, :);
        vDc23 = vDc23(I1:I2, :);
        vAc12 = vAc12(I1:I2, :);
        vAc23 = vAc23(I1:I2, :);
        
        fprintf('Limiting records to %s -- %s\n', ...
            EJ_library.cdf.tt2000_to_UTC_str(epoch(1)), ...
            EJ_library.cdf.tt2000_to_UTC_str(epoch(end)))
    end
    
    
    
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
        % IMPLEMENTATION NOTE: One could almost make a for loop over LFR
        % sampling frequencies (F0..F2) here. The current structure is however
        % useful for (1) making it possible to manually (and temporarily)
        % disable selected spectrograms, (2) setting frequency- and
        % channel-dependent constants.
        %=================
        % F0 spectrograms
        %=================
        pcfcList{end+1}     = @() (spectrogram_panel2( 'V1 DC', epoch, vDc1,  F_SAMPLE, 1, 'V1\_DC', [-11,-7]));
        if displayDcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 DC', epoch, vDc12, F_SAMPLE, 1, 'V12\_DC', [-13,-10]));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 DC', epoch, vDc23, F_SAMPLE, 1, 'V23\_DC', [-13,-10]));
        end
        if displayAcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 AC', epoch, vAc12, F_SAMPLE, 1, 'V12\_AC', [-13,-10]));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 AC', epoch, vAc23, F_SAMPLE, 1, 'V23\_AC', [-13,-10]));
        end
        %=================
        % F1 spectrograms
        %=================
        pcfcList{end+1} =     @() (spectrogram_panel2( 'V1 DC', epoch, vDc1,  F_SAMPLE, 2, 'V1\_DC', [-9,-5]));
        if displayDcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 DC', epoch, vDc12, F_SAMPLE, 2, 'V12\_DC', [-11,-8]));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 DC', epoch, vDc23, F_SAMPLE, 2, 'V23\_DC', [-11,-8]));
        end
        if displayAcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 AC', epoch, vAc12, F_SAMPLE, 2, 'V12\_AC', [-11,-8]));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 AC', epoch, vAc23, F_SAMPLE, 2, 'V23\_AC', [-11,-8]));
        end
        %=================
        % F2 spectrograms
        %=================
        pcfcList{end+1}     = @() (spectrogram_panel2( 'V1 DC', epoch, vDc1,  F_SAMPLE, 3, 'V1\_DC', [-8,-5]));
        if displayDcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 DC', epoch, vDc12, F_SAMPLE, 3, 'V12\_DC', [-10,-7]));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 DC', epoch, vDc23, F_SAMPLE, 3, 'V23\_DC', [-10,-7]));
        end
        if displayAcDiffs
            pcfcList{end+1} = @() (spectrogram_panel2('V12 AC', epoch, vAc12, F_SAMPLE, 3, 'V12\_AC', [-10,-7]));
            pcfcList{end+1} = @() (spectrogram_panel2('V23 AC', epoch, vAc23, F_SAMPLE, 3, 'V23\_AC', [-10,-7]));
        end
    end
    %===========================================================================
    % F0-F2 time series
    % IMPLEMENTATION NOTE: Panel tags have to be unique, or otherwise the axes
    % will be reused.
    %===========================================================================
    if displayDcDiffs
        %======================
        % DC single + DC diffs
        %======================
        SIGNALS_LEGEND_DC = EJ_library.graph.escape_str({'V1_DC','V12_DC','V23_DC'});
        tempFuncPtr = @(iLsf) (@() (time_series_panel2(...
            'V1,V12,V23 DC', epoch, {vDc1, vDc12, vDc23}, F_SAMPLE, iLsf, SIGNALS_LEGEND_DC)));
        
        pcfcList{end+1} = tempFuncPtr(1);
        pcfcList{end+1} = tempFuncPtr(2);
        pcfcList{end+1} = tempFuncPtr(3);

    end
    if ~displayDcDiffs && displayAcDiffs
        %======================
        % DC single + AC diffs
        %======================
        SIGNALS_LEGEND_DC_AC = EJ_library.graph.escape_str({'V1_DC','V12_AC','V23_AC'});
        tempFuncPtr = @(iLsf) (@() (time_series_panel2(...
            'V1,V12,V23 DC/AC', epoch, {vDc1, vAc12, vAc23}, F_SAMPLE, iLsf, SIGNALS_LEGEND_DC_AC)));

        pcfcList{end+1} = tempFuncPtr(1);
        pcfcList{end+1} = tempFuncPtr(2);
        pcfcList{end+1} = tempFuncPtr(3);

    end
    if displayDcDiffs && displayAcDiffs
        %======================
        % AC diffs (no single)
        %======================
        % NOTE: Assumes that DC single+diffs have already been plotted (in
        % separate panels).
        SIGNALS_LEGEND_AC = EJ_library.graph.escape_str({'V12_AC','V23_AC'});
        tempFuncPtr = @(iLsf) (@() (time_series_panel2(...
            'V12,V23 AC', epoch, {vAc12, vAc23}, F_SAMPLE, iLsf, SIGNALS_LEGEND_AC)));

        pcfcList{end+1} = tempFuncPtr(1);
        pcfcList{end+1} = tempFuncPtr(2);
        pcfcList{end+1} = tempFuncPtr(3);

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

    % For aligning MATLAB axes (taking color legends into account).
    irf_plot_axis_align(hAxesArray)
    % For aligning the content of the MATLAB axes.    
    irf_zoom(hAxesArray, 'x', irf.tint(epoch(1), epoch(end)))
    
end



% Convenient wrapper around spectrum_panel.
% Converts from zVar-like variables to what is actually used for plotting.
%
function h = spectrogram_panel2(panelTagSignalsStr, zvEpoch, zvData, zvSamplFreqHz, iLsf, trLegend, colLimits)
    samplFreqHz = EJ_library.so.constants.LSF_HZ(iLsf);
    lsfName     = EJ_library.so.constants.LSF_NAME_ARRAY{iLsf};
    
    bRecords = (zvSamplFreqHz == samplFreqHz);
    
    h = spectrogram_panel(...
        sprintf('%s %s spectrogram', panelTagSignalsStr, lsfName), ...
        zvEpoch(bRecords, :), ...
        zvData( bRecords, :), ...
        samplFreqHz, ...
        lsfName, ...
        trLegend, ...
        colLimits);
end



% ARGUMENTS
% =========
% tlLegend : Top-left  (TL) legend string.
% trLegend : Top-right (TR) legend string.
%
% SS  : SnapShot
% SSS : SnapShot Spectrogram
%
function h = spectrogram_panel(panelTag, zvEpoch, zvData, samplingFreqHz, tlLegend, trLegend, colLimits)
    % NOTE: Multiple-row labels causes trouble for the time series ylabels.
    % IMPLEMENTATION NOTE: Implemented to potentially be modified to handle TDS
    % snapshots that vary in length.

    % Fraction of the (minimum) time distance between snapshots (centers) that
    % will be used for displaying the spectra.
    % Value 1 : Spectras are adjacent between snapshot (for minimum snapshot
    % distance).
    SNAPSHOT_WIDTH_FRACTION  = 0.90;
    %SPECTRUM_OVERLAP_PERCENT = 0;    % Percent, not fraction.
    SPECTRUM_OVERLAP_PERCENT = 50;    % Percent, not fraction.
    
    % NOTE: More samples per spectrum is faster (sic!).
    N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.
    %N_SAMPLES_PER_SPECTRUM = 512;    % Speed test
    
    erikpgjohansson.assert.sizes(colLimits, [1,2])
    


    TsCa = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz);
    nTs  = numel(TsCa);

    h = irf_panel(panelTag);    
    
    %====================
    % Calculate spectras
    %====================
    % IMPLEMENTATION NOTE: irf_powerfft is the most time-consuming part of this
    % code.
    %
    % NOTE: Using for-->parfor speeds up plot_LFR_SWF by
    % 29.912231 s-->21.303145 s (irony). /2020-09-04
    %
    SpecrecCa = cell(nTs, 1);
    ssCenterEpochUnixArray = zeros(nTs, 1);
    parfor i = 1:nTs    % PARFOR
        Ts = TsCa{i};
        
        SpecrecCa{i} = irf_powerfft(Ts, N_SAMPLES_PER_SPECTRUM, samplingFreqHz, SPECTRUM_OVERLAP_PERCENT);
        
        % IMPLEMENTATION NOTE: Later needs the snapshot centers in the same time
        % system as Specrec.t (epoch Unix).
        ssCenterEpochUnixArray(i) = (Ts.time.start.epochUnix + Ts.time.stop.epochUnix)/2;
    end
    sssMaxWidthSecArray = derive_max_spectrum_width(ssCenterEpochUnixArray);
    
    %===========================================================================
    % Set the display locations of individual spectras (override defaults).
    % Separately stretch out the collection of spectras that stems from every
    % snapshot.
    % IMPLEMENTATION NOTE: This can not be done in the first loop in order to
    % derive the (minimum) snapshot time distance.
    %===========================================================================
    for i = 1:numel(TsCa)
        bKeep(i) = ~isempty(SpecrecCa{i});
        if ~isempty(SpecrecCa{i})
            %ssLengthSec = TsCa{i}.time.stop.epochUnix - TsCa{i}.time.start.epochUnix;
            
            sssWidthSec = sssMaxWidthSecArray(i) * SNAPSHOT_WIDTH_FRACTION;
            
            %SpecrecCa{i} = solo.ql.downsample_Specrec(SpecrecCa{i}, 10);    % TEST
            
            % Stretch out spectra (for given snapshot) in time to be ALMOST
            % adjacent between snapshots.
            % NOTE: Specrec.dt is not set by irf_powerfft so there is no default
            % value that can be scaled up.
            % NOTE: Uses original spectrum positions and re-positions them
            % relative to snapshot center.
            
            % Number of timestamps, but also spectras (within snapshot).
            nTime = numel(SpecrecCa{i}.t);
            % Distance from SS center to center of first/last FFT.
            distToSssEdgeT = sssWidthSec/2 - sssWidthSec/(2*nTime);
            SpecrecCa{i}.t  = ssCenterEpochUnixArray(i) + linspace(-distToSssEdgeT, distToSssEdgeT, nTime);
            SpecrecCa{i}.dt = ones(nTime, 1) * sssWidthSec/(2*nTime);
        end
    end
    
    SpecrecCa(~bKeep) = [];
    Specrec = merge_specrec(SpecrecCa);    
    
    Specrec.p_label = {'log_{10} [V^2/Hz]'};    % Replaces colorbarlabel
    irf_spectrogram(h, Specrec);   % Replaces irf_plot
    
    set(h, 'yscale','log')

    irf_legend(h, tlLegend, [0.02 0.98], 'color', 'k')
    irf_legend(h, trLegend, [0.98 0.98])
    
    % Technically potentially slower to load file every time, but I don't want
    % to make this an argument, and the OS probable caches the small file.
    load cmap
    colormap(cmap)
    
    caxis(h, colLimits)
    set(h, 'YTick', [0.1, 1, 10, 100])
end



% For every snapshot, return the available width (in time; centered on snapshot
% center) for displaying the snapshot spectrogram. Time offset and unit
% unimportant. Argument and return values have same unit.
%
function sssMaxWidthArray = derive_max_spectrum_width(ssCenterArray)
    % Use distance to nearest snapshot for each snapshot separately.
    % NOTE: Should NOT be multiplied by two, since using entire distance.
    sssMaxWidthArray = min([Inf; diff(ssCenterArray(:))], [diff(ssCenterArray(:)); Inf]);
    
    % Use smallest distance between any two consecutive snapshots (one global
    % value for all snapshots). Sometimes yields too narrow spectrograms.
    % Ex: solo_L2_rpw-lfr-surv-swf-e-cdag_20200228_V01.cdf
    %sssMaxWidthArray = min(diff(ssCenterArray)) * ones(size(ssCenterArray));
    
    % NOTE: Can not assume that both input and output have same size, only same
    % length.
    assert(numel(ssCenterArray) == numel(sssMaxWidthArray))
end



% Convenient wrapper around time_series_panel.
% Converts from zVar-like variables (N samples/record; all records) to what is
% actually used for plotting.
function h = time_series_panel2(panelTagSignalsStr, zvEpoch, zvDataList, zvSamplFreqHz, iLsf, trLegend)
    
    bRecords = (zvSamplFreqHz == EJ_library.so.constants.LSF_HZ(iLsf));
    samplFreqHz = EJ_library.so.constants.LSF_HZ(iLsf);
    lsfName     = EJ_library.so.constants.LSF_NAME_ARRAY{iLsf};
    
    nSps       = size(zvDataList{1}, 2);   % SPS = Samples Per Snapshot
    
    zvEpoch  = zvEpoch(bRecords);
    nRecords = size(zvEpoch, 1);   % NOTE: After selecting records.
    zvEpoch  = EJ_library.so.convert_N_to_1_SPR_Epoch(zvEpoch, nSps, ones(nRecords, 1)*samplFreqHz);
    
    for i = 1:numel(zvDataList)
        zvData        = zvDataList{i}(bRecords, :);
        zvDataList{i} = EJ_library.so.convert_N_to_1_SPR_redistribute(zvData);
    end
    
    % NOTE: Effectively serves as an assertion on zv sizes.
    Ts = irf.ts_scalar(zvEpoch, [zvDataList{:}]);
    
    % PLOT
    h = time_series_panel(...
        sprintf('%s %s time series', panelTagSignalsStr, lsfName), ...
        Ts, lsfName, trLegend);
end



% ARGUMENTS
% =========
% tlLegend : Top-left  (TL) legend.
% trLegend : Top-right (TR) legend. Cell array of strings, one per scalar time
%            series.
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
%           IMPLEMENTATION NOTE: Can not(?) be struct array since MATLAB
%           confuses indexing a TSeries array (with brackets) with some special
%           TSeries functionality for calling its code with brackets (calling
%           TSeries' method "subsref").
%
% IMPLEMENTATION NOTE: Function is written to some day be easily extended to be
% used for use with TDS's length-varying snapshots.
%
function TsCa = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz)
    %
    % NOTE: No special treatment of snapshots with only NaN.
    
    assert(isscalar(samplingFreqHz))
    EJ_library.assert.sizes(zvData, [NaN, NaN])
    assert(size(zvEpoch, 1) == size(zvData, 1))   % Same number of records
    bicas.proc_utils.assert_zv_Epoch(zvEpoch)
    
    nRecords = size(zvData, 1);
    nSps     = size(zvData, 2);
    assert(nSps >= 2)
    
    % Relative timestamps inside CDF record/snapshot.
    epochRelArray = int64([0:(nSps-1)] * 1/samplingFreqHz * 1e9);
    
    TsCa = {};
    for i = 1:nRecords
        epochRecord = zvEpoch(i) + epochRelArray;
        TsCa{i}  = irf.ts_scalar(epochRecord, zvData(i, :));
    end
    
end



% Merge multiple instances of "specrec" structs as returned by irf_powerfft,
% with identical frequencies.
% NOTE: Optionally added fields must be added after merging.
% NOTE: Cf EJ_library.utils.merge_Specrec which is more powerful but which is
% unnecessary here since only merging spectras with the same frequencies
% (underlying data uses the same sampling frequency).
%
% ARGUMENTS
% =========
% SpecrecCa : Cell array of "Specrec" structs as returned by irf_powerfft, but
%             with .dt (column array) added to it.
%             NOTE: Requires dt (column array of scalars).
%             NOTE: Assumes that all specrec use the same frequencies.
%             IMPLEMENTATION NOTE: Uses cell array instead of struct array to be
%             able to handle (and ignore) the case specrec = [] which can be
%             returned by irf_powerfft.
%
% RETURN VALUE
% ============
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
            assert(iscolumn(S.dt), 'S.dt is not a column.')
            assert(numel(S.dt) == numel(S.t), 'Badly formatted SpecrecCa{%i}.', i)
            
            Specrec.f    = S.f;                       % NOTE: Not adding to array, but setting it in its entirety.
            Specrec.p{1} = [Specrec.p{1}; S.p{1}];    % NOTE: Add to array inside cell array.
            Specrec.t    = [Specrec.t;    S.t(:)];    % NOTE: Has to be column vector.
            Specrec.dt   = [Specrec.dt;   S.dt(:)];   % NOTE: Has to be column vector.
        end
    end
    
    assert(issorted(Specrec.t))   % Not sure if good assertion.
end



function data = get_CDF_zv_data(D, zvName, i3)
    
    % TEMPORARY: For backward compatibility.
    if strcmp(zvName, 'SAMPLING_RATE') && isfield(D.data, 'F_SAMPLE')
        zvName = 'F_SAMPLE';
    elseif strcmp(zvName, 'VDC') && isfield(D.data, 'V')
        zvName = 'V';
    elseif strcmp(zvName, 'EDC') && isfield(D.data, 'E')
        zvName = 'E';
    end
    
    fillValue = getfillval(D, zvName);
    data = D.data.(zvName).data(:, :, i3);
    data = changem(data, NaN, fillValue);
end
