%
% Quicklook for the content of one BIAS HK dataset (CDF file), i.e. DATASET_ID = SOLO_HK_RPW-BIA.
%
%
% INCOMPLETE
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2020-01-27.
%
function plot_HK(filePath)
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
    % ==> 38 zVars
    %
    % PROPOSAL: New name: plot_BIAS_HK
    % PROPOSAL: Combine triplets of bit series.
    % PROPOSAL: Improve bit plotting. The locations of bit flips are technically somewhat misrepresented. E.g. singular
    % ones (surrounded by zeros) should be represented as "infinitely" thin vertical lines, whereas two ones in a row
    % (surrounded by zeros) are not.
    %
    % ~BUG: No grid lines for bits empty zVars.
    % BUG: Multiple dates below shared x axis, after zooming (automatic globa re-setting of x axis?).
    %
    % PROBLEM: How handle triplets of scalar/bit zVars (as list of zVars)?
    %   PROPOSAL: Put in lists as one string=3 zVars?
    %   PROPOSAL: Put in lists as three recursively grouped zVars (cell array in cell array).
    %
    % PROPOSAL: ylabel [TM], [TM units].
    
    %warning('Incomplete quicklook code')
    
    BIT_PANEL_HEIGHT = 0.011;
    
    ZVAR_EXCLUDED_LIST = {...
        'Epoch', ...
        'HK_BIA_MODE_VERSION_NR', ...
        'ACQUISITION_TIME_UNITS', ...
        'ACQUISITION_LABEL', ...
        'PA_DPU_HK_REPORT_SID', ...
        'HK_BIA_CUR_SELECTED_PAGE', ...
        'HK_BIA_DUMMY', ...
        'PA_RPW_HK_SPARE8_1'};
    
    ZVAR_SCALAR_LIST = {
        'HK_BIA_MODE_MUX_SET', ...
        'HK_BIA_REF_VOLTAGE_H', ...
        'HK_BIA_REF_VOLTAGE_L', ...
        'HK_BIA_TEMP_PCB', ...
        'HK_BIA_NHV', ...
        'HK_BIA_PHV', ...
        'HK_BIA_REF2_VOLT'};

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
    ZVAR_ALL_LIST = union(union(ZVAR_EXCLUDED_LIST, ZVAR_SCALAR_LIST), ZVAR_BIT_FLAG_LIST);
    assert(numel(ZVAR_ALL_LIST) == numel(ZVAR_EXCLUDED_LIST)+numel(ZVAR_SCALAR_LIST)+numel(ZVAR_BIT_FLAG_LIST), 'zVars overlap')
    
    D = dataobj(filePath);
    Epoch = D.data.Epoch.data;
    
    hAxesArray = irf_plot(numel(ZVAR_SCALAR_LIST)+3+numel(ZVAR_BIT_FLAG_LIST), 'newfigure');
    % Panel height: fixed size + weight for each panel, used for distributing segments.
    panelHeightFswArray = zeros(0,2);   % FSW = Fixed Size + Weight

    %=====================================
    % Add panels for scalar numeric zVars
    %=====================================
    for i = 1:numel(ZVAR_SCALAR_LIST)
        lineWidth = 0.5;  % Default.
        if strcmp(ZVAR_SCALAR_LIST{i}, 'HK_BIA_MODE_MUX_SET')
            lineWidth = 5.0;
        end
        
        [hAxes, hLines, hLegendText] = plot_time_series1(D, ZVAR_SCALAR_LIST{i}, lineWidth);        
        
        if strcmp(ZVAR_SCALAR_LIST{i}, 'HK_BIA_MODE_MUX_SET')
            % IMPLEMENTATION NOTE: Skip lowest tick/label value (YTick), so that the label does not overlap with the
            % label of the panel below.
            set(hAxes, 'YLim', [0,7])
            set(hAxes, 'YTick', 1:7)
        end
        
        panelHeightFswArray(end+1, :) = [0, 1];
    end
    
    %==================================================================
    % Add panels for triplets of numeric zVars (one zVar per antenna).
    %==================================================================
    plot_time_series3(D, 'HK_BIA_BIAS%s');             panelHeightFswArray(end+1, :) = [0, 1];
    plot_time_series3(D, 'HK_BIA_TEMP_ANT%s_LF_PA');   panelHeightFswArray(end+1, :) = [0, 1];
    plot_time_series3(D, 'HK_BIA_M%s'   );             panelHeightFswArray(end+1, :) = [0, 1];

    %==================================================================
    % Add panels for bit-valued zVars.
    %==================================================================
    hBitArray = [];
    for i = 1:numel(ZVAR_BIT_FLAG_LIST)
        zVarName = ZVAR_BIT_FLAG_LIST{i};
        hBitArray(end+1) = plot_bit_series(zVarName, Epoch, D.data.(zVarName).data, zVarName);
        panelHeightFswArray(end+1, :) = [BIT_PANEL_HEIGHT, 0];
    end



    %=================================
    % Adjust the height of bit panels
    %=================================
    % 'Position' : [left bottom width height]. Size and location, excluding a margin for the labels.
    positionCa = get(hAxesArray, 'Position');    % CA = Cell Array
    yPanelArray1      = cellfun(@(x) ([x(2)]), positionCa);
    % Panel height before distributing height segments. Assumes that panels are adjacent to each other.
    heightPanelArray1 = cellfun(@(x) ([x(4)]), positionCa);
    
%     iPanelsToSetHeightFor = numel(ZVAR_SCALAR_LIST) + 3 + [1:numel(ZVAR_BIT_FLAG_LIST)];
%     heightPanelArray2 = solo.ql.reweight(...
%         heightPanelArray1, ...
%         BIT_PANEL_HEIGHT, ...
%         iPanelsToSetHeightFor);
    heightPanelArray2 = EJ_library.utils.distribute_segments(...
        sum(heightPanelArray1), ...
        panelHeightFswArray(:,1), ...
        panelHeightFswArray(:,2));
    yPanelArray2 = cumsum([heightPanelArray2(2:end); yPanelArray1(end)], 'reverse');
    
    %config = [repmat(-1, 1, numel(ZVAR_SCALAR_LIST)+3), repmat(BIT_PANEL_HEIGHT, 1, numel(ZVAR_BIT_FLAG_LIST))]
    %[heightPanelArray2] = EJ_library.utils.autoscale(sum(heightPanelArray1), config);
    
    for i = 1:numel(hAxesArray)
        position = positionCa{i};
        position(2) = yPanelArray2(i);
        position(4) = heightPanelArray2(i);
        set(hAxesArray(i), 'InnerPosition', position)
    end
    
    
    
    solo.ql.set_std_title('BIAS HK', filePath, hAxesArray(1))
    
    irf_plot_axis_align(hAxesArray)               % For aligning MATLAB axes OuterPosition (taking color legends into account).
    irf_zoom(hAxesArray, 'x', irf.tint(Epoch))    % For aligning the content of the MATLAB axes.

    set(hAxesArray(1:end-1), 'XLabel', [])        % Remove duplicate x labels. Empirically: Must come after irf_zoom.
    %EJ_library.graph.set_shared_dynamic_XYZAxes(hAxesArray, 'X', 'No init')    % Test
end



% Plot one panel for 1 numeric time series.
%
function [hAxes, hLines, hLegendText] = plot_time_series1(D, zvName, lineWidth)
    % NOTE: Automatically derive panel tag.
    panelTag = zvName;
    
    Ts = irf.ts_scalar(D.data.Epoch.data, D.data.(zvName).data);
    [hAxes, hLines, hLegendText] = plot_time_series(panelTag, Ts, '', {EJ_library.graph.escape_str(zvName)});
    set(hLines, 'LineWidth', lineWidth)
end



% Plot one panel for 3 numeric time series (one for each antenna).
%
% ARGUMENTS
% =========
% D               : dataobj
% zVarNamePattern : String used in sprintf, with "%s" (not "%i"!) representing the antenna number 1, 2, or 3.
%
function [hAxes, hLines, hLegendText] = plot_time_series3(D, zVarNamePattern)
    % NOTE: Automatically derive panel tag.
    panelTag = sprintf(zVarNamePattern, 'x');
    
    zvName1 = sprintf(zVarNamePattern, '1');
    zvName2 = sprintf(zVarNamePattern, '2');
    zvName3 = sprintf(zVarNamePattern, '3');
    
    Ts = irf.ts_scalar(D.data.Epoch.data, [...
        D.data.(zvName1).data, ...
        D.data.(zvName2).data, ...
        D.data.(zvName3).data]);

    [hAxes, hLines, hLegendText] = plot_time_series(panelTag, Ts, '', {EJ_library.graph.escape_str(sprintf(zVarNamePattern, '1')), '2', '3'});
end



function [hAxes, hLines, hLegendText] = plot_time_series(panelTag, Ts, yLabel, channelNames)
    LEGEND_POSITION = [0.98, 0.98];
    %LEGEND_COLOR = [0,0,1];
    
    hAxes  = irf_panel(panelTag);
    hLines = irf_plot(hAxes, Ts);   % Array, if multiple scalar time series.
    modify_legend_text_color(hLines)
    
    ylabel(hAxes, yLabel)
    hLegendText = irf_legend(hAxes, channelNames, LEGEND_POSITION);
    %modify_legend_text_color(hLegendText)
end



% Plot one bit-valued array in one irf_panel.
%
% ARGUMENTS
% =========
% channelName : 
%
function [hAxes, hLines, hLegendText] = plot_bit_series(panelTag, Epoch, bitSeries, channelName)
    % IMPLEMENTATION NOTE: Somewhat hackish plotting of "area": y value=0 is replaced by NaN. Plots very wide line which
    % is later clipped/cropped by the axes.
    
    %BIT_COLOR = [0.5, 0.5, 0.5];
    BIT_COLOR = [0, 0, 0];
    %CHANNEL_NAME_POS = [0.02, 0.98];
    CHANNEL_NAME_POS = [0.98, 0.98];
    %LEGEND_COLOR = [0,0,1];
    
    if isempty(bitSeries)
        Epoch(:) = [];
        channelName = ['(', channelName, ': empty zVar)'];
    end
    
    % Create time series object.
    bitSeries = double(bitSeries);      % NOTE: Raw zVariable may be integers which can not be assigned the value NaN.
    assert(all(ismember(bitSeries, [0,1])))
    bitSeries(bitSeries == 0) = NaN;
    Ts = irf.ts_scalar(Epoch, bitSeries);
    
    % Plot
    hAxes = irf_panel(panelTag);
    hLines = irf_plot(hAxes, Ts, 'LineWidth', 1e9, 'Color', BIT_COLOR);
    modify_legend_text_color(hLines)
    
    % NOTE: Command puts the text relative to the specified coordinates in different ways depending on coordinates.
    hLegendText = irf_legend(hAxes, EJ_library.graph.escape_str(channelName), CHANNEL_NAME_POS);    
    %modify_legend_text_color(hLegendText)
    
    set(hAxes, 'Clipping', 'on')   % Should be default, so not necessary.
    set(hAxes, 'ClippingStyle', 'rectangle')
    set(hAxes, 'yTickLabel', {})
    set(hAxes, 'YGrid', 'off')   % Remove horizontal grid lines.
    set(hAxes, 'TickLength', [0, 0])          % Remove ticks by setting their length.
end



function modify_legend_text_color(hTextArray)
    C_FADE = 0.7;
    for i = 1:numel(hTextArray)
        legendColor = get(hTextArray(i), 'Color');
        set(hTextArray(i), 'Color', 1-C_FADE*(1-legendColor));   % Fade color (move toward white).
    end
end

