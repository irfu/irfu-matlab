%
% Quicklook for the content of one BIAS LFR CWF dataset (CDF file), ie. any of
% the following DATASET_IDs:
%   SOLO_L2_RPW-LFR-SBM1-CWF-E
%   SOLO_L2_RPW-LFR-SBM2-CWF-E
%   SOLO_L2_RPW-LFR-SURV-CWF-E
%
%
% NOTE: Color scale is log, therefore negative values (probably).
%
% 
% SELECTED BUGFIXES
% =================
% NOTE: List is kept for now to keep track of fixed bugs, since not all
% quicklooks have been re-generated, and old quicklooks may display the effects
% of already fixed bugs.
% --
% FIXED BUG: Empirically: When there is any EAC data, for any time interval,
% then summay plot only contains VDC1-spectrum for the time interval (VDC1) for
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
% Can be very slow due to spectrograms.
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
% -rw-r--r-- 1 erjo solarorbiter 807M 2020-08-20 22.19:43 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200808_V01.cdf
% -rw-r--r-- 1 erjo solarorbiter 807M 2020-08-21 07.55:23 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200809_V02.cdf
% -rw-r--r-- 1 erjo solarorbiter  95M 2020-08-20 19.31:17 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200810_V01.cdf
% -rw-r--r-- 1 erjo solarorbiter 115M 2020-08-20 21.27:53 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200811_V01.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-20 21.01:52 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200812_V02.cdf
% -rw-r--r-- 1 erjo solarorbiter 104M 2020-08-20 18.16:54 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200813_V01.cdf
% -rw-r--r-- 1 erjo solarorbiter  67M 2020-08-26 12.31:51 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200819_V03.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-26 14.19:01 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200820_V03.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-26 13.12:55 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200821_V03.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-26 15.03:52 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200822_V02.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-26 18.19:03 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200823_V02.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-28 13.53:03 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200824_V03.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-28 14.01:53 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200825_V02.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-28 14.07:03 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200826_V01.cdf
% -rw-r--r-- 1 erjo solarorbiter 105M 2020-08-29 14.02:19 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200827_V02.cdf
% -rw-r--r-- 1 erjo solarorbiter  62M 2020-08-29 15.35:15 08/solo_L2_rpw-lfr-surv-cwf-e-cdag_20200828_V01.cdf
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-28.
%
function hAxesArray = plot_LFR_CWF(filePath)
    
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
    % PROPOSAL: Settings for disabling spectrum etc.
    
    %PERMIT_SIMULTANEOUS_DC_AC_DIFFS = 0;

    D = dataobj(filePath);
    epoch = D.data.Epoch.data;
    vDc1          = get_CDF_zv_data(D, 'VDC', 1);
    vDc12         = get_CDF_zv_data(D, 'EDC', 1);
    vDc23         = get_CDF_zv_data(D, 'EDC', 3);
    vAc12         = get_CDF_zv_data(D, 'EAC', 1);
    vAc23         = get_CDF_zv_data(D, 'EAC', 3);
    zvSamplFreqHz = D.data.SAMPLING_RATE.data;
    clear D
    
    if 0
        % DEBUG: Limit records
        I1 = 1.12e5;
        I2 = 1.20e5;
        
        epoch = epoch(I1:I2);
        vDc1  = vDc1( I1:I2);
        vDc12 = vDc12(I1:I2);
        vDc23 = vDc23(I1:I2);
        vAc12 = vAc12(I1:I2);
        vAc23 = vAc23(I1:I2);
        zvSamplFreqHz = zvSamplFreqHz(I1:I2);
        
        fprintf('Limiting records to %s -- %s\n', ...
            EJ_library.cdf.tt2000_to_UTC_str(epoch(1)), ...
            EJ_library.cdf.tt2000_to_UTC_str(epoch(end)))
    end
    
    TsVdc1  = irf.ts_scalar(epoch, vDc1);
    TsVdc12 = irf.ts_scalar(epoch, vDc12);
    TsVdc23 = irf.ts_scalar(epoch, vDc23);
    TsVac12 = irf.ts_scalar(epoch, vAc12);
    TsVac23 = irf.ts_scalar(epoch, vAc23);
    
    irf_plot(3+5,'newfigure');
    
    hAxesArray = [];
    hAxesArray(end+1) = spectrogram_panel( 'V1 DC spectrogram', TsVdc1,  zvSamplFreqHz, 'V1\_DC');
    hAxesArray(end+1) = spectrogram_panel('V12 DC spectrogram', TsVdc12, zvSamplFreqHz, 'V12\_DC');
    hAxesArray(end+1) = spectrogram_panel('V23 DC spectrogram', TsVdc23, zvSamplFreqHz, 'V23\_DC');

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
    
    irf_plot_axis_align(hAxesArray)               % For aligning MATLAB axes (taking color legends into account).
    irf_zoom(hAxesArray, 'x', irf.tint(epoch))    % For aligning the content of the MATLAB axes.
end



% ARGUMENTS
% =========
% panelTag      : 
% Ts            : irfumatlab TSeries (volt).
% yLabelNonUnit : y label without unit (unit is at the color bar; Assumes "Ts" uses volt).
%
function hAxes = spectrogram_panel(panelTag, Ts, zvSamplingFreqHz, yLabelNonUnit)
    
    N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.
    
    hAxes = irf_panel(panelTag);
    
    SpecrecCa = {};
    [iSs1Array, iSs2Array, nSs] = EJ_library.utils.split_by_change(zvSamplingFreqHz);
    for jSs = 1:nSs
        
        iSs = iSs1Array(jSs) : iSs2Array(jSs);

        SpecrecCa{jSs} = irf_powerfft(Ts(iSs), N_SAMPLES_PER_SPECTRUM, zvSamplingFreqHz(iSs1Array(jSs)));
        
%         Specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
%         hold on    % Required so that successive calls to irf_spectrogram do not replace preceedings ones.
%         irf_spectrogram(hAxes, Specrec);   % Replaces irf_plot
    end

    Specrec = EJ_library.utils.merge_Specrec(SpecrecCa);
    Specrec.p_label = {'[V^2/Hz]'};    % Replaces colorbarlabel
    irf_spectrogram(hAxes, Specrec);   % Replaces irf_plot    
    
    
    set(hAxes, 'yscale','log')
    ylabel(hAxes, {yLabelNonUnit; 'f [Hz]'})   % NOTE: Adding frequency unit on separate row.

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
    
    % TEMPORARY: For backward compatibility.
    if strcmp(zvName, 'VDC') && isfield(D.data, 'V')
        zvName = 'V';
    elseif strcmp(zvName, 'EDC') && isfield(D.data, 'E')
        zvName = 'E';
    end
    
    fillValue = getfillval(D, zvName);
    data = D.data.(zvName).data(:, i2);
    data = changem(data, NaN, fillValue);
end
