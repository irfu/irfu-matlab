%
% Generate official quicklook for the content of one BIAS LFR SWF dataset (CDF
% file), i.e. DATASET_ID = SOLO_L2_RPW-LFR-SURV-SWF-E.
%
%
% NOTES
% =====
% * CAN TAKE A SIGNIFICANT AMOUNT OF TIME TO COMPLETE, DUE TO SPECTROGRAMS(?).
% * Only capable (default) of only showing either DC diffs or AC diffs.
%   There are hard-coded settings for permitting or forcing both.
% * Does not yet support spectrogram overlap.
% * Time series panels interpolate between snapshots.
% * Color scale is log. Therefore colors represent negative values (probably).
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
    % POLICY: Combined BOGIQ for ~all quicklook plot code: See solo.sp.summary_plot.
    %
    % TODO-DEC: How submit data to spectrum_panel?
    %   NEED: Should also work for TDS's varying-length snapshots?
    %   NOTE: Might be necessary to also split up CWF to speed up spectrograms eventually.
    %   TODO-DEC: Should caller split up data into snapshots?
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
    % PROPOSAL: If no EDC or EAC (or both), print that somehow on plot.
    %
    % ~BUG: Algorithm for calculating enlargement of snapshot spectra not always appropriate.
    %   Ex: solo_L2_rpw-lfr-surv-swf-e-cdag_20200228_V01.cdf
    %   PROPOSAL: Enlarge each snapshot spectrum until it reaches neighbour (roughly).
    %       Use min(timeToSnapshotBefore, timeToSnapshotafter) as max radius.
    %
    % Old TODO: YK's "fixup" 2020-10-13 -- DONE
    %   hswf = solo.sp.plot_LFR_SWF([RPWPATH LRFFILE]);%% fixup
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
    
    
    % Constants for how to handle certain special cases (unusual datasets)
    % --------------------------------------------------------------------
    % YK 2020-04-16: Officially only either DC or AC diffs.
    % NOTE: solo_L2_rpw-lfr-surv-cwf-e-cdag_20200228_V01.cdf contains both DC &
    % AC diffs.
    ALWAYS_SIMULTANEOUS_DC_AC_DIFFS_PLOTS = 0;   % DEFAULT 0. Useful for debugging (runs through all code).
    PERMIT_SIMULTANEOUS_DC_AC_DIFFS       = 1;
    ENABLE_SPECTROGRAMS                   = 1;   % DEFAULT 1. Useful for debugging non-spectrogram code.
    PERMIT_NEITHER_DC_AC_DIFFS            = 1;   % DEFAULT 1?
    
    D = dataobj(filePath);
    
    zvEpoch       = D.data.Epoch.data;
    zvSamplRateHz = get_CDF_zv_data(D, 'SAMPLING_RATE', 1);
    zvDc1         = get_CDF_zv_data(D, 'VDC', 1);
    zvDc12        = get_CDF_zv_data(D, 'EDC', 1);
    zvDc23        = get_CDF_zv_data(D, 'EDC', 3);
    zvAc12        = get_CDF_zv_data(D, 'EAC', 1);
    zvAc23        = get_CDF_zv_data(D, 'EAC', 3);
    clear D
    
    hasDcDiffs = any(~isnan(zvDc12(:))) || any(~isnan(zvDc23(:)));
    hasAcDiffs = any(~isnan(zvAc12(:))) || any(~isnan(zvAc23(:)));
    
    % ASSERTIONS
    if      hasDcDiffs && hasAcDiffs && ~PERMIT_SIMULTANEOUS_DC_AC_DIFFS
        error('Dataset (CDF file) contains both DC diff and AC diff data. Can not handle this case.')
    elseif ~hasDcDiffs && ~hasAcDiffs
        MSG = 'Dataset (CDF file) contains neither DC diff nor AC diff data.';
        if PERMIT_NEITHER_DC_AC_DIFFS
            warning('Dataset (CDF file) contains neither DC diff nor AC diff data.')
        else
            error('%s. Can not handle this case.', MSG)
        end
    end

    zvDcAc12 = solo.sp.utils.merge_zvs(zvDc12, zvAc12);
    zvDcAc23 = solo.sp.utils.merge_zvs(zvDc23, zvAc23);
    
    %=================================================================
    % Determine whether DC diffs, AC diffs, or both should be plotted
    %=================================================================
    if ALWAYS_SIMULTANEOUS_DC_AC_DIFFS_PLOTS
        displayDcDiffs = 1;
        displayAcDiffs = 1;
    else
        
        displayDcDiffs = hasDcDiffs;
        displayAcDiffs = hasAcDiffs;
    end
    
    
    
    Sp = solo.sp.summary_plot();



    if ENABLE_SPECTROGRAMS
        % IMPLEMENTATION NOTE: One could almost make a for loop over LFR
        % sampling frequencies (F0..F2) here. The current structure is however
        % useful for (1) making it possible to manually (and temporarily)
        % disable selected spectrograms, (2) setting frequency- and
        % channel-dependent constants.
        %=================
        % F0 spectrograms
        %=================
        Sp.add_panel_spectrogram_SWF_LSF('V1 DC',     zvEpoch, zvDc1,    zvSamplRateHz, 1, 'V1\_DC',     [-11, -7]);
        Sp.add_panel_spectrogram_SWF_LSF('V12 DC/AC', zvEpoch, zvDcAc12, zvSamplRateHz, 1, 'V12\_DC/AC', [-13,-10]);
        Sp.add_panel_spectrogram_SWF_LSF('V23 DC/AC', zvEpoch, zvDcAc23, zvSamplRateHz, 1, 'V23\_DC/AC', [-13,-10]);
        %=================
        % F1 spectrograms
        %=================
        Sp.add_panel_spectrogram_SWF_LSF('V1 DC',     zvEpoch, zvDc1,    zvSamplRateHz, 2, 'V1\_DC',     [ -9,-5]);
        Sp.add_panel_spectrogram_SWF_LSF('V12 DC/AC', zvEpoch, zvDcAc12, zvSamplRateHz, 2, 'V12\_DC/AC', [-11,-8]);
        Sp.add_panel_spectrogram_SWF_LSF('V23 DC/AC', zvEpoch, zvDcAc23, zvSamplRateHz, 2, 'V23\_DC/AC', [-11,-8]);
        %=================
        % F2 spectrograms
        %=================
        Sp.add_panel_spectrogram_SWF_LSF('V1 DC',     zvEpoch, zvDc1,    zvSamplRateHz, 3, 'V1\_DC',     [ -8,-5]);
        Sp.add_panel_spectrogram_SWF_LSF('V12 DC/AC', zvEpoch, zvDcAc12, zvSamplRateHz, 3, 'V12\_DC/AC', [-10,-7]);
        Sp.add_panel_spectrogram_SWF_LSF('V23 DC/AC', zvEpoch, zvDcAc23, zvSamplRateHz, 3, 'V23\_DC/AC', [-10,-7]);
    end
    
    
    
    %===========================================================================
    % F0-F2 time series
    % -----------------
    % IMPLEMENTATION NOTE: Panel tags have to be unique, or otherwise the axes
    % will be reused.
    %===========================================================================
    if displayDcDiffs
        %======================
        % DC single + DC diffs
        %======================
        for iLsf = 1:3
            Sp.add_panel_time_series_SWF_LSF(...
                'V1,V12,V23 DC', zvEpoch, {zvDc1, zvDc12, zvDc23}, zvSamplRateHz, iLsf, ...
                irf.graph.escape_str({'V1_DC','V12_DC','V23_DC'}), [0,0,0]);
        end

    end
    if ~displayDcDiffs && displayAcDiffs
        %======================
        % DC single + AC diffs
        %======================
        for iLsf = 1:3
            Sp.add_panel_time_series_SWF_LSF(...
                'V1,V12,V23 DC/AC', zvEpoch, {zvDc1, zvAc12, zvAc23}, zvSamplRateHz, iLsf, ...
                irf.graph.escape_str({'V1_DC','V12_AC','V23_AC'}), [0,1,1]);
        end
    end
    if displayDcDiffs && displayAcDiffs
        %======================
        % AC diffs (no single)
        %======================
        % NOTE: Assumes that DC single+diffs have already been plotted (in
        % separate panels). There should be another "if" statement for that.
        for iLsf = 1:3
            Sp.add_panel_time_series_SWF_LSF(...
                'V12,V23 AC', zvEpoch, {zvAc12, zvAc23}, zvSamplRateHz, iLsf, ...
                irf.graph.escape_str({'V12_AC','V23_AC'}), [1,1]);
        end
    end

    hAxesArray = Sp.finalize('LFR SWF L2', filePath);
    
end



function zv = get_CDF_zv_data(D, zvName, i3)
    
    % TEMPORARY: For backward compatibility.
    if strcmp(zvName, 'SAMPLING_RATE') && isfield(D.data, 'F_SAMPLE')
        zvName = 'F_SAMPLE';
    elseif strcmp(zvName, 'VDC') && isfield(D.data, 'V')
        zvName = 'V';
    elseif strcmp(zvName, 'EDC') && isfield(D.data, 'E')
        zvName = 'E';
    end
    
    fv = getfillval(D, zvName);
    zv = D.data.(zvName).data(:, :, i3);
    zv = changem(zv, NaN, fv);
end
