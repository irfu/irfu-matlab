%
% Quicklook for the content of one BIAS HK dataset (CDF file), i.e. DATASET_ID = SOLO_HK_RPW-BIA.
%
%
% VERY INCOMPLETE
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
    % PA_DPU_HK_REPORT_SID      CDF_UINT1/1       0:[]    T/
    % HK_BIA_BIAS1              CDF_UINT2/1       0:[]    T/
    % HK_BIA_BIAS2              CDF_UINT2/1       0:[]    T/
    % HK_BIA_BIAS3              CDF_UINT2/1       0:[]    T/
    % HK_BIA_M1                 CDF_UINT2/1       0:[]    T/
    % HK_BIA_M2                 CDF_UINT2/1       0:[]    T/
    % HK_BIA_M3                 CDF_UINT2/1       0:[]    T/
    % HK_BIA_REF_VOLTAGE_H      CDF_UINT2/1       0:[]    T/
    % HK_BIA_REF_VOLTAGE_L      CDF_UINT1/1       0:[]    T/
    % HK_BIA_TEMP_ANT1_LF_PA    CDF_UINT2/1       0:[]    T/
    % HK_BIA_TEMP_ANT2_LF_PA    CDF_UINT2/1       0:[]    T/
    % HK_BIA_TEMP_ANT3_LF_PA    CDF_UINT2/1       0:[]    T/
    % HK_BIA_TEMP_PCB           CDF_UINT2/1       0:[]    T/
    % HK_BIA_NHV                CDF_UINT1/1       0:[]    T/
    % HK_BIA_PHV                CDF_UINT1/1       0:[]    T/
    % HK_BIA_REF2_VOLT          CDF_UINT2/1       0:[]    T/
    % HK_BIA_MODE_VERSION_NR    CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_ACTIVE_LINK   CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_SWEEP_BUSY    CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_MUX_SET       CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_HV_ENABLED    CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_BIAS3_ENABLED CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_BIAS2_ENABLED CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_BIAS1_ENABLED CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_DIFF_PROBE    CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_BYPASS_PROBE3 CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_BYPASS_PROBE2 CDF_UINT1/1       0:[]    T/
    % HK_BIA_MODE_BYPASS_PROBE1 CDF_UINT1/1       0:[]    T/
    % HK_BIA_DIFF_GAIN          CDF_UINT1/1       0:[]    T/
    % HK_BIA_CMD_CNT            CDF_UINT1/1       0:[]    T/
    % HK_BIA_CUR_SELECTED_PAGE  CDF_UINT1/1       0:[]    T/
    % HK_BIA_DUMMY              CDF_UINT2/1       0:[]    T/
    % PA_RPW_HK_SPARE8_1        CDF_UINT1/1       0:[]    T/

    % TODO-DECISION: Content of figure title
    %   PROPOSAL: ~File path, ~filename
    %   PROPOSAL: Time range
    %   PROPOSAL: DATASET_ID
    %
    % PROPOSAL: Improve bit plotting. The locations of bit flips are technically somewhat misrepresented. E.g. singular
    % ones (surrounded by zeros) should be represented as "infinitely" thin vertical lines, whereas two ones in a row
    % (surrounded by zeros) are not.
    %
    % ~BUG: No grid lines for bits empty zVars.
    
    warning('Incomplete quicklook')
    
    D = dataobj(filePath);
    Epoch = D.data.Epoch.data;
    
%     HK_BIA_REF_VOLTAGE_H
%     HK_BIA_REF_VOLTAGE_L

    % zVariables that (presumably) are bit flags.
    ZVAR_BIT_FLAG_LIST = {...
        'TIME_SYNCHRO_FLAG', ...
        'HK_BIA_MODE_ACTIVE_LINK', ...
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

    hAxesArray = irf_plot(3+numel(ZVAR_BIT_FLAG_LIST)+1, 'newfigure');
    
    plot_time_series3(D, 'HK_BIA_BIAS%s');
    plot_time_series3(D, 'HK_BIA_M%s'   );
    plot_time_series3(D, 'HK_BIA_TEMP_ANT%s_LF_PA');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    testBit = round(rand(size(Epoch)));
    %TsTestBit = irf.ts_scalar(Epoch, testBit);
    
    hBitArray = [];
    for i = 1:numel(ZVAR_BIT_FLAG_LIST)
        zVarName = ZVAR_BIT_FLAG_LIST{i};
        hBitArray(end+1) = plot_bit_series(zVarName, Epoch, D.data.(zVarName).data, zVarName);
    end
    
    hBitArray(end+1) = plot_bit_series('BIT TEST', Epoch, testBit, 'BIT_TEST');
    
%     for i = 1:numel(hBitArray)
%         op = get(hBitArray(i), 'OuterPosition');
%         %op(4) = 0.06;
%         set(hBitArray(i), 'OuterPosition', op, 'YGrid', 'off', 'YTick', [])
%     end

    set(hBitArray, 'YGrid', 'off', 'YTick', [])

    %heightArray = cellfun(@(x) (x(4)), get(hAxesArray, 'Position'));
    %heightSum = sum(heightArray)
    %arrayfun(@(h) (set_axes_height(h,0.10)), hAxesArray);


    irf_plot_axis_align(hAxesArray)               % For aligning MATLAB axes OuterPosition (taking color legends into account).

    irf_zoom(hAxesArray, 'x', irf.tint(Epoch))    % For aligning the content of the MATLAB axes.

    %print_axes_height(hAxesArray)
    title(hAxesArray(1), 'BIAS HK')
    %print_axes_height(hAxesArray)
end



% function x = reweight(x, xNew, iSet)
%     xSumNotSet = sum(x) - x(iSet);
%     x(iSet) = xNew;
%     x / sum(x) * xSumNotSet;
% end



function print_axes_height(hAxesArray)
    heightArray = cellfun(@(x) (x(4)), get(hAxesArray, 'Position'));
    heightArray, size(heightArray)
end



function set_axes_height(h, height)
    % NOTE: Must use Position, not OuterPosition.
    pos = get(h, 'Position');
    pos(4) = height;
    set(h, 'Position', pos)
end



% Plot3 time series, one for each antenna.
%
% zVarNamePattern : String used in sprintf, with "%s" (!) representing the antenna number 1, 2, or 3.
function h = plot_time_series3(D, zVarNamePattern)
    panelTag = sprintf(zVarNamePattern, 'x');
    yLabel   = sprintf(zVarNamePattern, 'x');
    
    zvName1 = sprintf(zVarNamePattern, '1');
    zvName2 = sprintf(zVarNamePattern, '2');
    zvName3 = sprintf(zVarNamePattern, '3');
    
    Ts = irf.ts_scalar(D.data.Epoch.data, [...
        D.data.(zvName1).data, ...
        D.data.(zvName2).data, ...
        D.data.(zvName3).data]);
    
    %plot_time_series(panelTag, Ts, yLabel, {'1', '2', '3'})
    h = plot_time_series(panelTag, Ts, '', {escape_str(sprintf(zVarNamePattern, '1')), '2', '3'});
end



function h = plot_bit_series(panelTag, Epoch, bitSeries, channelName)
    % IMPLEMENTATION NOTE: Somewhat hackish plotting of "area": y value=0 is replaced by NaN. Plots very wide line which
    % is later clipped/cropped by the axes.
    
    BIT_COLOR = [0.5,0.5,0.5];
    
    if isempty(bitSeries)
        Epoch(:) = [];
        channelName = ['(', channelName, ': empty zVar)'];
    end
    bitSeries = double(bitSeries);      % NOTE: Raw zVariable may be integers which can not be assigned the value NaN.
    assert(all(ismember(bitSeries, [0,1])))
    bitSeries(bitSeries == 0) = NaN;    
    Ts = irf.ts_scalar(Epoch, bitSeries);
    
    h = irf_panel(panelTag);
    irf_plot(h, Ts, 'linewidth', 1e9, 'color', BIT_COLOR)
    irf_legend(h, escape_str(channelName), [0.02 0.98])
    
    set(h, 'Clipping', 'on')   % Should be default, so not necessary.
    set(h, 'ClippingStyle', 'rectangle')
    set(h, 'yTickLabel', {})
end



function h = plot_time_series(panelTag, Ts, yLabel, channelNames)
    h = irf_panel(panelTag);
    irf_plot(h, Ts)
    ylabel(h, yLabel)
    irf_legend(h, channelNames, [0.98 0.98])
end



% Also works for cell arrays of strings.
function s = escape_str(s)
    s = strrep(s, '_', '\_');
end
