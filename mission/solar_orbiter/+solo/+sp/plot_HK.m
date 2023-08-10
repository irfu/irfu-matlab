%
% Generate official quicklook for the content of one BIAS HK dataset (CDF file),
% i.e. DATASET_ID = SOLO_HK_RPW-BIA.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-27.
%
function hAxesArray = plot_HK(filePath)
    % SOLO_HK_RPW-BIA_V01 zVariables:
    %
    % Variable Information (0 rVariable, 38 zVariables)
    % ===========================================================
    % Epoch                     CDF_TT2000/1      0:[]    T/
    % ACQUISITION_TIME          CDF_UINT4/1       1:[2]   T/T
    % SYNCHRO_FLAG              CDF_UINT1/1       0:[]    T/
    % ACQUISITION_TIME_UNITS    CDF_CHAR/16       1:[2]   F/T
    % ACQUISITION_TIME_LABEL    CDF_CHAR/32       1:[2]   F/T
    % PA_DPU_HK_REPORT_SID      CDF_UINT1/1       0:[]    T/   var "HK Report SID."  Meaningful?
    % HK_BIA_BIAS1              CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_BIAS2              CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_BIAS3              CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_M1                 CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_M2                 CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_M3                 CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_REF_VOLTAGE_H      CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_REF_VOLTAGE_L      CDF_UINT1/1       0:[]    T/   var
    % HK_BIA_TEMP_ANT1_LF_PA    CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_TEMP_ANT2_LF_PA    CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_TEMP_ANT3_LF_PA    CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_TEMP_PCB           CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_NHV                CDF_UINT1/1       0:[]    T/   var
    % HK_BIA_PHV                CDF_UINT1/1       0:[]    T/   var
    % HK_BIA_REF2_VOLT          CDF_UINT2/1       0:[]    T/   var
    % HK_BIA_MODE_VERSION_NR    CDF_UINT1/1       0:[]    T/   var, in principle, but very constant.
    % HK_BIA_MODE_ACTIVE_LINK   CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_SWEEP_BUSY    CDF_UINT1/1       0:[]    T/   bit?  "Counter of parity error."
    % HK_BIA_MODE_MUX_SET       CDF_UINT1/1       0:[]    T/   var
    % HK_BIA_MODE_HV_ENABLED    CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_BIAS3_ENABLED CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_BIAS2_ENABLED CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_BIAS1_ENABLED CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_DIFF_PROBE    CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_BYPASS_PROBE3 CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_BYPASS_PROBE2 CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_MODE_BYPASS_PROBE1 CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_DIFF_GAIN          CDF_UINT1/1       0:[]    T/   bit
    % HK_BIA_CMD_CNT            CDF_UINT1/1       0:[]    T/   var   "Command counter (counts number of writes made to BIAS)." Meaningful?
    % HK_BIA_CUR_SELECTED_PAGE  CDF_UINT1/1       0:[]    T/   var?  "HK_BIA_CUR_SELECTED_PAGE
    % HK_BIA_DUMMY              CDF_UINT2/1       0:[]    T/   var   "Last written value to serial link."
    % PA_RPW_HK_SPARE8_1        CDF_UINT1/1       0:[]    T/   var?  "Unused."
    % ==> 38 ZVs
    %
    % PROPOSAL: New name: plot_BIAS_HK
    %
    % PROBLEM: How handle triplets of scalar/bit ZVs (as list of ZVs)?
    %   PROPOSAL: Put in lists as one string=3 ZVs?
    %   PROPOSAL: Put in lists as three recursively grouped ZVs (cell array in cell array).
    %
    % TODO-DEC: How plot bits?
    %   PROPOSAL: Combine triplets of bit series, somehow.
    %       PROPOSALS: Three graphs together in same panel
    %           PROPOSAL: Slight offset between them so that lines do not
    %                     overlap when horizontal.
    %               CON: Will overlap when vertical, i.e. entirely for short
    %                    flips.
    %           PROPOSAL: Great offset between them so that lines do not
    %                     overlap when horizontal or vertical.
    %   PROPOSAL: Flip between colors. No line.
    %       CON: Flip might be invisible for very short durations
    %           Ex: Old bug: HK_BIA_MODE_BYPASS_PROBE1/2/3 does not show bit=1
    %               (short time interval; calibration sweep) for 2020-07-01 when
    %               saving to file.
    %   PROPOSAL: Graph.
    %       PRO: Can always see flips, no matter how short the time interval.
    %       CON: Hard to distinguish 0 and 1, low & high, especially if many bit plots in a row.
    %           PRO: May blend with boundary between plots (panels).
    %   PROPOSAL: Graph with different colors/symbols for 0 and 1 (technically two
    %             graphs).
    %       PROBLEM/CON: Unclear if this works with using irf_plot + TimeSeries.
    %
    % BUG: Multiple dates below shared x axis, after zooming (automatic global re-setting of x axis?).
    %
    % 
    % PROPOSAL: Return hAxesArray like other plot_* functions.
    %
    % PROPOSAL: ylabel [TM], [TM units].
    

    
    BIT_PANEL_HEIGHT = 0.011;
    %BIT_PANEL_HEIGHT = 0.03;   % DEBUG
    
    ZVAR_EXCLUDED_LIST = {...
        'Epoch', ...
        'HK_BIA_MODE_VERSION_NR', ...
        'ACQUISITION_TIME_UNITS', ...
        'ACQUISITION_LABEL', ...
        'PA_DPU_HK_REPORT_SID', ...
        'HK_BIA_CUR_SELECTED_PAGE', ...
        'HK_BIA_DUMMY', ...
        'PA_RPW_HK_SPARE8_1'};
    
    ZVAR_SINGLE_SCALAR_LIST = {
        'HK_BIA_MODE_MUX_SET', ...
        'HK_BIA_REF_VOLTAGE_H', ...
        'HK_BIA_REF_VOLTAGE_L', ...
        'HK_BIA_TEMP_PCB', ...
        'HK_BIA_NHV', ...
        'HK_BIA_PHV', ...
        'HK_BIA_REF2_VOLT'};
    
    ZVAR_TRIPLET_SCALAR_LIST = {
        zv_name_triplet('HK_BIA_BIAS%i'), ...
        zv_name_triplet('HK_BIA_TEMP_ANT%i_LF_PA'), ...
        zv_name_triplet('HK_BIA_M%i')};
    
    % zVariables that are (presumably) bit flags.
    ZVAR_BIT_FLAG_LIST = {...
        'SYNCHRO_FLAG', ...
        'HK_BIA_MODE_ACTIVE_LINK', ...
        'HK_BIA_MODE_SWEEP_BUSY', ...
        'HK_BIA_MODE_HV_ENABLED', ...
        'HK_BIA_MODE_DIFF_PROBE', ...
        'HK_BIA_MODE_BYPASS_PROBE1', ...
        'HK_BIA_MODE_BYPASS_PROBE2', ...
        'HK_BIA_MODE_BYPASS_PROBE3', ...
        'HK_BIA_MODE_BIAS1_ENABLED', ...
        'HK_BIA_MODE_BIAS2_ENABLED', ...
        'HK_BIA_MODE_BIAS3_ENABLED', ...
        'HK_BIA_DIFF_GAIN' ...
        };
    
    % ASSERTION
    ZVAR_ALL_LIST = irf.utils.union(...
        ZVAR_EXCLUDED_LIST, ...
        ZVAR_SINGLE_SCALAR_LIST, ...
        unpack(ZVAR_TRIPLET_SCALAR_LIST), ...
        ZVAR_BIT_FLAG_LIST);
    
    assert(...
        numel(ZVAR_ALL_LIST) == ...
        numel(ZVAR_EXCLUDED_LIST) + ...
        numel(ZVAR_SINGLE_SCALAR_LIST) + ...
        numel(unpack(ZVAR_TRIPLET_SCALAR_LIST)) + ...
        numel(ZVAR_BIT_FLAG_LIST), ...
        'Separate sets of specified ZVs overlap.')
    
    D = dataobj(filePath);
    Epoch = D.data.Epoch.data;

    

    Sp = solo.sp.summary_plot();
    
    
    
    %=====================================
    % Add panels for scalar numeric ZVs
    %=====================================
    for i = 1:numel(ZVAR_SINGLE_SCALAR_LIST)
        lineWidth = 0.5;  % Default.
        if strcmp(ZVAR_SINGLE_SCALAR_LIST{i}, 'HK_BIA_MODE_MUX_SET')
            lineWidth = 5.0;
        end
        
        axesPropCa = {};
        if strcmp(ZVAR_SINGLE_SCALAR_LIST{i}, 'HK_BIA_MODE_MUX_SET')
            % IMPLEMENTATION NOTE: Skip lowest tick/label value (YTick), so that the label does not overlap with the
            % label of the panel below.
            axesPropCa = {'YLim', [0,7], 'YTick', 1:7};
        end
        Sp.add_panel_time_series1_HK(D, ZVAR_SINGLE_SCALAR_LIST{i}, {'LineWidth', lineWidth}, axesPropCa)
    end
    
    %==================================================================
    % Add panels for triplets of numeric ZVs (one zVar per antenna).
    %==================================================================
    for i = 1:numel(ZVAR_TRIPLET_SCALAR_LIST)
        Sp.add_panel_time_series3_HK(D, ZVAR_TRIPLET_SCALAR_LIST{i});
    end

    %==================================================================
    % Add panels for bit-valued ZVs.
    %==================================================================
%     hBitArray = [];
    for i = 1:numel(ZVAR_BIT_FLAG_LIST)
        zVarName = ZVAR_BIT_FLAG_LIST{i};
        Sp.add_panel_plot_bit_series(zVarName, Epoch, D.data.(zVarName).data, zVarName, BIT_PANEL_HEIGHT);
    end


    
    hAxesArray = Sp.finalize('BIAS HK', filePath);

end



% Utility function
%
% c  : Cell array of (1) non-cells, and/or (2) cell arrays.
% ca : Column cell array.
function ca = unpack(c)
    ca = cell(0,1);
    for i = 1:numel(c)
        if iscell(c{i})
            ca = [ca; c{i}(:)];
        else
            ca = [ca; c(i)];
        end
    end
end



% Utility function
function zVarNameTripletCa = zv_name_triplet(zVarNamePattern)
    zVarNameTripletCa = {};
    for iAnt = 1:3
        zVarNameTripletCa{iAnt} = sprintf(zVarNamePattern, iAnt);
    end
end
