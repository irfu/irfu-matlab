%
% Quicklook for the content of one BIAS HK dataset (CDF file), i.e. DATASET_ID =
% SOLO_HK_RPW-BIA.
%
%
% INCOMPLETE?
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
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
    %
    % PROBLEM: How handle triplets of scalar/bit zVars (as list of zVars)?
    %   PROPOSAL: Put in lists as one string=3 zVars?
    %   PROPOSAL: Put in lists as three recursively grouped zVars (cell array in cell array).
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
    %
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
    ZVAR_ALL_LIST = EJ_library.utils.union(...
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
        'zVars overlap')
    
    D = dataobj(filePath);
    Epoch = D.data.Epoch.data;
    
    hAxesArray = irf_plot(...
        numel(ZVAR_SINGLE_SCALAR_LIST)+...
        numel(ZVAR_TRIPLET_SCALAR_LIST) + ...
        numel(ZVAR_BIT_FLAG_LIST), ...
        'newfigure');
    % Panel height: fixed size + weight for each panel, used for distributing segments.
    panelHeightFswArray = zeros(0,2);   % FSW = Fixed Size + Weight

    %=====================================
    % Add panels for scalar numeric zVars
    %=====================================
    for i = 1:numel(ZVAR_SINGLE_SCALAR_LIST)
        lineWidth = 0.5;  % Default.
        if strcmp(ZVAR_SINGLE_SCALAR_LIST{i}, 'HK_BIA_MODE_MUX_SET')
            lineWidth = 5.0;
        end
        
        [hAxes, hLines, hLegendText] = plot_time_series1(D, ZVAR_SINGLE_SCALAR_LIST{i}, lineWidth);        
        
        if strcmp(ZVAR_SINGLE_SCALAR_LIST{i}, 'HK_BIA_MODE_MUX_SET')
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
    for i = 1:numel(ZVAR_TRIPLET_SCALAR_LIST)
        
        plot_time_series3(D, ZVAR_TRIPLET_SCALAR_LIST{i});
        panelHeightFswArray(end+1, :) = [0, 1];
    end

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
    
    heightPanelArray2 = EJ_library.utils.distribute_segments(...
        sum(heightPanelArray1), ...
        panelHeightFswArray(:,1), ...
        panelHeightFswArray(:,2));
    yPanelArray2 = cumsum([heightPanelArray2(2:end); yPanelArray1(end)], 'reverse');
    
    for i = 1:numel(hAxesArray)
        position = positionCa{i};
        position(2) = yPanelArray2(i);
        position(4) = heightPanelArray2(i);
        set(hAxesArray(i), 'InnerPosition', position)
    end
    
    
    
    solo.ql.set_std_title('BIAS HK', filePath, hAxesArray(1))
    
    % For aligning MATLAB axes OuterPosition (taking color legends into account).
    irf_plot_axis_align(hAxesArray)               
    % For aligning the content of the MATLAB axes.
    irf_zoom(hAxesArray, 'x', irf.tint(Epoch))

    % Remove duplicate x labels.
    % Empirically: Must come after irf_zoom.
    set(hAxesArray(1:end-1), 'XLabel', [])
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
% D           : dataobj
% zVarNamesCa : Length 3 cell array of zVar names.
%
function [hAxes, hLines, hLegendText] = plot_time_series3(D, zVarNamesCa)
% zVarNamePattern : String used in sprintf, with "%s" (not "%i"!) representing the antenna number 1, 2, or 3.
    assert(numel(zVarNamesCa) == 3)
    
    % NOTE: Automatically derive panel tag. Exact value unimportant.
    panelTag = zVarNamesCa{1};
    
    Ts = irf.ts_scalar(D.data.Epoch.data, [...
        D.data.(zVarNamesCa{1}).data, ...
        D.data.(zVarNamesCa{2}).data, ...
        D.data.(zVarNamesCa{3}).data]);

    [hAxes, hLines, hLegendText] = plot_time_series(panelTag, Ts, '', ...
        {EJ_library.graph.escape_str(zVarNamesCa{1}), '2', '3'});
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
function [hAxes, hLines, hLegendText] = plot_bit_series(panelTag, Epoch, b, channelName)
    % OBSOLETE: IMPLEMENTATION NOTE: Somewhat hackish plotting of "area": y value=0 is
    % replaced by NaN. Plots very wide line which is later clipped/cropped by
    % the axes.

    %BIT_COLOR = [0.5, 0.5, 0.5];
    %BIT_COLOR = [0, 0, 0];
    %BIT_COLOR = [0, 0, 1];
    %CHANNEL_NAME_POS = [0.02, 0.98];
    CHANNEL_NAME_POS = [0.98, 0.98];
    %LEGEND_COLOR = [0,0,1];
    
    % ASSERTIONS
    assert(all(ismember(b, [0,1])))
    
    if isempty(b)
        Epoch(:) = [];
        channelName = ['(', channelName, ': empty zVar)'];
    end
    
   

    % Divide b into b0 and b1, which when plotted as lines, form one line
    % identical to b. This way one can plot separate parts with separate colors.
    % Use NaN for parts that are not plotted by either b0 or b1.
    b = double(b);   % Must be double to be able to assign NaN.
    % Identify indices just after b has flipped.
    iJumpUp = 1 + find((b(1:end-1)==0) & (b(2:end)==1));
    iJumpDn = 1 + find((b(1:end-1)==1) & (b(2:end)==0));
    % Create b0, b1.
    b0 = b;
    b1 = b;
    b0(b0 == 1) = NaN;
    b1(b1 == 0) = NaN;
    b0(iJumpUp) = 1;
    b1(iJumpDn) = 0;
    
    % Create time series object
    % -------------------------
    Ts = irf.ts_scalar(Epoch, [b0, b1]);
    
    % Plot
    hAxes = irf_panel(panelTag);
    hLines = irf_plot(hAxes, Ts, 'LineWidth', 3);
    assert(numel(hLines) == 2)
    hLines(1).Color = [0,0,1];
    hLines(2).Color = [1,0,0];
    
    modify_legend_text_color(hLines)
    
    
    
    % NOTE: Command puts the text relative to the specified coordinates in
    % different ways depending on coordinates.
    hLegendText = irf_legend(hAxes, EJ_library.graph.escape_str(channelName), CHANNEL_NAME_POS);    
    %modify_legend_text_color(hLegendText)
    
    set(hAxes, 'Clipping', 'on')   % Should be default, so not necessary.
    set(hAxes, 'ClippingStyle', 'rectangle')
    set(hAxes, 'yTickLabel', {})
    set(hAxes, 'YGrid', 'off')          % Remove horizontal grid lines.
    %set(hAxes, 'TickLength', [0, 0])    % Remove ticks by setting their length.
    
    set(hAxes, 'YLim', [-0.2, 1.2])
    set(hAxes, 'YTick', [])
    %set(hAxes, 'yTickLabel', {'0', '1'})
end



function modify_legend_text_color(hTextArray)
    C_FADE = 0.7;
    for i = 1:numel(hTextArray)
        legendColor = get(hTextArray(i), 'Color');
        set(hTextArray(i), 'Color', 1-C_FADE*(1-legendColor));   % Fade color (move toward white).
    end
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