%
% Generate official quicklook for the content of one BIAS LFR CWF dataset (CDF
% file), ie. any of the following DATASET_IDs:
%   SOLO_L2_RPW-LFR-SBM1-CWF-E
%   SOLO_L2_RPW-LFR-SBM2-CWF-E
%   SOLO_L2_RPW-LFR-SURV-CWF-E
%
%
% NOTE: CAN BE VERY SLOW DUE TO SPECTROGRAMS.
% NOTE: Color scale is log. Therefore colors represent negative values (probably).
%
%
% SELECTED BUGFIXES
% =================
% NOTE: List is kept for now to keep track of fixed bugs, since not all
% quicklooks have been re-generated, and old quicklooks may display the effects
% of already fixed bugs.
% --
% FIXED BUG: Empirically: When there is any EAC data, for any time interval,
% then summary plot only contains VDC1-spectrum for the time interval (VDC1) for
% which there is EAC data, and plots no EDC spectras.
%   Ex: 2020-04-09: LFR-CWF 107 MiB
%   Ex: 2020-05-07: LFR-CWF 284 MiB
%   Ex: 2020-05-08: LFR-CWF 391 MiB
%   Ex: 2020-07-08: LFR-CWF 243 miB
% Seems to be due to how samplFreqHz is selected for irf_powerfft. Code
% selects the most common one, weighted by records (not time). Data with
% lower sampling frequencies then come to be counted as having many
% data gaps. ==> NaN-valued spectra.
% (Probably the same also for other EDC12/EAC12, EDC23/EAC23.)
%
%
% BUGS
% ====
% CAN BE VERY SLOW DUE TO SPECTROGRAMS.
% YK 2020-09-01: Plotting of spectrograms slows it down and that spectrograms
% should be averaged before plotting.
% --
% EJ 2020-09-01: Run implies/hints that wall time increases much faster than dataset size.
% erjo@brain /home/erjo/temp/sp/2020-08-31> grep cwf.*Wall  so_sp_batch.2020-08-31_12.52.27.log
%     solo_L3_rpw-lfr-surv-cwf-e_20200808_V01.png: Wall time used for plotting: 38127.7 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200809_V01.png: Wall time used for plotting: 37890.2 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200810_V01.png: Wall time used for plotting: 136.988 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200811_V01.png: Wall time used for plotting: 114.704 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200812_V01.png: Wall time used for plotting: 132.951 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200813_V01.png: Wall time used for plotting: 133.784 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200819_V01.png: Wall time used for plotting: 66.5982 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200820_V01.png: Wall time used for plotting: 134.468 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200821_V01.png: Wall time used for plotting: 116.474 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200822_V01.png: Wall time used for plotting: 138.96 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200823_V01.png: Wall time used for plotting: 116.116 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200824_V01.png: Wall time used for plotting: 132.725 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200825_V01.png: Wall time used for plotting: 135.157 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200826_V01.png: Wall time used for plotting: 133.753 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200827_V01.png: Wall time used for plotting: 134.075 [s]
%     solo_L3_rpw-lfr-surv-cwf-e_20200828_V01.png: Wall time used for plotting: 79.0095 [s]
% erjo@brain /data/solo/remote/data/L2/lfr_wf_e/2020> ll -h */*cwf*202008{08..31}*
% /.../
%  807M 2020-08-20 22.19:43 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200808_V01.cdf
%  807M 2020-08-21 07.55:23 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200809_V02.cdf
%   95M 2020-08-20 19.31:17 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200810_V01.cdf
%  115M 2020-08-20 21.27:53 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200811_V01.cdf
%  105M 2020-08-20 21.01:52 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200812_V02.cdf
%  104M 2020-08-20 18.16:54 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200813_V01.cdf
%   67M 2020-08-26 12.31:51 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200819_V03.cdf
%  105M 2020-08-26 14.19:01 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200820_V03.cdf
%  105M 2020-08-26 13.12:55 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200821_V03.cdf
%  105M 2020-08-26 15.03:52 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200822_V02.cdf
%  105M 2020-08-26 18.19:03 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200823_V02.cdf
%  105M 2020-08-28 13.53:03 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200824_V03.cdf
%  105M 2020-08-28 14.01:53 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200825_V02.cdf
%  105M 2020-08-28 14.07:03 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200826_V01.cdf
%  105M 2020-08-29 14.02:19 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200827_V02.cdf
%   62M 2020-08-29 15.35:15 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200828_V01.cdf
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-28.
%
function hAxesArray = plot_LFR_CWF(filePath)
    %
    % TODO-NI: Detect DC/AC and only plot one of them? Analogous to LFR SWF.
    %   TODO-NI: CWF never uses AC? SBM1/2?
    %
    % PROPOSAL: Settings arguments to disable/enable hidden functionality
    %   Ex: SWF: Permit/force DC+AC diffs
    %   Ex: Disable spectrograms.
    %
    % PROPOSAL: Settings for disabling spectrum etc.

    D = dataobj(filePath);
    zvEpoch       = D.data.Epoch.data;
    zvDc1         = get_CDF_ZV_data(D, 'VDC', 1);
    zvDc12        = get_CDF_ZV_data(D, 'EDC', 1);
    zvDc23        = get_CDF_ZV_data(D, 'EDC', 3);
    zvAc12        = get_CDF_ZV_data(D, 'EAC', 1);
    zvAc23        = get_CDF_ZV_data(D, 'EAC', 3);
    zvSamplFreqHz = D.data.SAMPLING_RATE.data;
    clear D    % To prove that D will not be used later.

    zvAc12 = remove_mean(zvAc12);
    zvAc23 = remove_mean(zvAc23);

    zvDcAc12 = solo.sp.utils.merge_zvs(zvDc12, zvAc12);
    zvDcAc23 = solo.sp.utils.merge_zvs(zvDc23, zvAc23);

    Sp = solo.sp.summary_plot();

    Sp.add_panel_spectrogram_CWF( 'V1 DC spectrogram',    zvEpoch, zvDc1,    zvSamplFreqHz, 'V1\_DC',     [-7, -3]);
    Sp.add_panel_spectrogram_CWF('V12 DC/AC spectrogram', zvEpoch, zvDcAc12, zvSamplFreqHz, 'V12\_DC/AC', [-8, -4]);
    % Good color scale? Which use?
    % Was [-9, -5] for DC only till YK suggested adding AC. /2021-03-15
    % Ex: 2020-08-06, 2020-06-29: high intensities ==> max
    Sp.add_panel_spectrogram_CWF('V23 DC/AC spectrogram', zvEpoch, zvDcAc23, zvSamplFreqHz, 'V23\_DC/AC', [-8.5, -6.0]);

    Sp.add_panel_time_series_CWF( 'V1 DC time series',    zvEpoch,  zvDc1,           'V1_DC [V]');
    Sp.add_panel_time_series_CWF('V12 DC/AC time series', zvEpoch, [zvDc12, zvAc12], 'V12_DC/AC [V]', 'trLegend', {'DC', 'AC'});
    Sp.add_panel_time_series_CWF('V23 DC/AC time series', zvEpoch, [zvDc23, zvAc23], 'V23_DC/AC [V]', 'trLegend', {'DC', 'AC'});

    hAxesArray = Sp.finalize('LFR CWF L2', filePath);

end



% IMPLEMENTATION NOTE: Must be able to handle NaN.
% Empirically: Using detrend() returns NaN if any samples are NaN.
function x = remove_mean(x)
    x = x - mean(x, 'omitnan');
end



function zv = get_CDF_ZV_data(D, zvName, i2)

    % TEMPORARY: For backward compatibility.
    if strcmp(zvName, 'VDC') && isfield(D.data, 'V')
        zvName = 'V';
    elseif strcmp(zvName, 'EDC') && isfield(D.data, 'E')
        zvName = 'E';
    end

    fv = getfillval(D, zvName);
    zv = D.data.(zvName).data(:, i2);
    zv = changem(zv, NaN, fv);
end
