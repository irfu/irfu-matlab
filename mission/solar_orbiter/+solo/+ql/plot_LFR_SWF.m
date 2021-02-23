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
    % POLICY: Combined BOGIQ for ~all quicklook plot code: See solo.ql.summary_plot.
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
    % Old TODO: YK's "fixup" 2020-10-13 -- DONE
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
    %   2020-10-13: Should already have been implemented.
    %
    % PROPOSAL: Rename. Not (MATLAB) plot, but (MATLAB) figure.
    
    

    % YK 2020-04-16: Officially only either DC or AC diffs.
    % NOTE: solo_L2_rpw-lfr-surv-cwf-e-cdag_20200228_V01.cdf contains both DC &
    % AC diffs.
    ALWAYS_SIMULTANEOUS_DC_AC_DIFFS_PLOTS = 0;   % DEFAULT 0. Useful for debugging (runs through all code).
    PERMIT_SIMULTANEOUS_DC_AC_DIFFS       = 1;
    ENABLE_SPECTROGRAMS                   = 1;   % DEFAULT 1. Useful for debugging non-spectrogram code.
    
    D = dataobj(filePath);
    
    epoch    = D.data.Epoch.data;
    F_SAMPLE = get_CDF_zv_data(D, 'SAMPLING_RATE', 1);
    vDc1     = get_CDF_zv_data(D, 'VDC', 1);
    vDc12    = get_CDF_zv_data(D, 'EDC', 1);
    vDc23    = get_CDF_zv_data(D, 'EDC', 3);
    vAc12    = get_CDF_zv_data(D, 'EAC', 1);
    vAc23    = get_CDF_zv_data(D, 'EAC', 3);
    clear D
    
%     if 0
%         % DEBUG: Limit records
%         I1 = 1;
%         I2 = 100;
%         
%         epoch = epoch(I1:I2);
%         F_SAMPLE = F_SAMPLE(I1:I2);
%         vDc1  = vDc1( I1:I2, :);
%         vDc12 = vDc12(I1:I2, :);
%         vDc23 = vDc23(I1:I2, :);
%         vAc12 = vAc12(I1:I2, :);
%         vAc23 = vAc23(I1:I2, :);
%         
%         fprintf('Limiting records to %s -- %s\n', ...
%             EJ_library.cdf.TT2000_to_UTC_str(epoch(1)), ...
%             EJ_library.cdf.TT2000_to_UTC_str(epoch(end)))
%     end
    
    
    
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
    
    
    
    Sp = solo.ql.summary_plot();
    
    
    
    if ENABLE_SPECTROGRAMS
        % IMPLEMENTATION NOTE: One could almost make a for loop over LFR
        % sampling frequencies (F0..F2) here. The current structure is however
        % useful for (1) making it possible to manually (and temporarily)
        % disable selected spectrograms, (2) setting frequency- and
        % channel-dependent constants.
        %=================
        % F0 spectrograms
        %=================
        Sp.add_panel_spectrogram_SWF_LSF( 'V1 DC', epoch, vDc1,  F_SAMPLE, 1, 'V1\_DC', [-11,-7]);
        if displayDcDiffs
            Sp.add_panel_spectrogram_SWF_LSF('V12 DC', epoch, vDc12, F_SAMPLE, 1, 'V12\_DC', [-13,-10]);
            Sp.add_panel_spectrogram_SWF_LSF('V23 DC', epoch, vDc23, F_SAMPLE, 1, 'V23\_DC', [-13,-10]);
        end
        if displayAcDiffs
            Sp.add_panel_spectrogram_SWF_LSF('V12 AC', epoch, vAc12, F_SAMPLE, 1, 'V12\_AC', [-13,-10]);
            Sp.add_panel_spectrogram_SWF_LSF('V23 AC', epoch, vAc23, F_SAMPLE, 1, 'V23\_AC', [-13,-10]);
        end
        %=================
        % F1 spectrograms
        %=================
        Sp.add_panel_spectrogram_SWF_LSF( 'V1 DC', epoch, vDc1,  F_SAMPLE, 2, 'V1\_DC', [-9,-5]);
        if displayDcDiffs
            Sp.add_panel_spectrogram_SWF_LSF('V12 DC', epoch, vDc12, F_SAMPLE, 2, 'V12\_DC', [-11,-8]);
            Sp.add_panel_spectrogram_SWF_LSF('V23 DC', epoch, vDc23, F_SAMPLE, 2, 'V23\_DC', [-11,-8]);
        end
        if displayAcDiffs
            Sp.add_panel_spectrogram_SWF_LSF('V12 AC', epoch, vAc12, F_SAMPLE, 2, 'V12\_AC', [-11,-8]);
            Sp.add_panel_spectrogram_SWF_LSF('V23 AC', epoch, vAc23, F_SAMPLE, 2, 'V23\_AC', [-11,-8]);
        end
        %=================
        % F2 spectrograms
        %=================
        Sp.add_panel_spectrogram_SWF_LSF( 'V1 DC', epoch, vDc1,  F_SAMPLE, 3, 'V1\_DC', [-8,-5]);
        if displayDcDiffs
            Sp.add_panel_spectrogram_SWF_LSF('V12 DC', epoch, vDc12, F_SAMPLE, 3, 'V12\_DC', [-10,-7]);
            Sp.add_panel_spectrogram_SWF_LSF('V23 DC', epoch, vDc23, F_SAMPLE, 3, 'V23\_DC', [-10,-7]);
        end
        if displayAcDiffs
            Sp.add_panel_spectrogram_SWF_LSF('V12 AC', epoch, vAc12, F_SAMPLE, 3, 'V12\_AC', [-10,-7]);
            Sp.add_panel_spectrogram_SWF_LSF('V23 AC', epoch, vAc23, F_SAMPLE, 3, 'V23\_AC', [-10,-7]);
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
        for iLsf = 1:3
            Sp.add_panel_time_series_SWF_LSF(...
                'V1,V12,V23 DC', epoch, {vDc1, vDc12, vDc23}, F_SAMPLE, iLsf, ...
                EJ_library.graph.escape_str({'V1_DC','V12_DC','V23_DC'}), [0,0,0]);
        end

    end
    if ~displayDcDiffs && displayAcDiffs
        %======================
        % DC single + AC diffs
        %======================
        for iLsf = 1:3
            Sp.add_panel_time_series_SWF_LSF(...
                'V1,V12,V23 DC/AC', epoch, {vDc1, vAc12, vAc23}, F_SAMPLE, iLsf, ...
                EJ_library.graph.escape_str({'V1_DC','V12_AC','V23_AC'}), [0,1,1]);
        end
    end
    if displayDcDiffs && displayAcDiffs
        %======================
        % AC diffs (no single)
        %======================
        % NOTE: Assumes that DC single+diffs have already been plotted (in
        % separate panels).
        for iLsf = 1:3
            Sp.add_panel_time_series_SWF_LSF(...
                'V12,V23 AC', epoch, {vAc12, vAc23}, F_SAMPLE, iLsf, ...
                EJ_library.graph.escape_str({'V12_AC','V23_AC'}), [1,1]);
        end
    end

    hAxesArray = Sp.finalize('LFR SWF L2', filePath);
    
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
