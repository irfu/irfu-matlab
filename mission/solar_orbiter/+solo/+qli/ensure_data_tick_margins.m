%
% Given (1) the min/max ticks, and (2) min/max data value on some axis, derive
% suggested plot min/max limits to use so that ticks are not at the plot min/max
% limits. The algorithm largely assumes that the axis is logarithmic but it
% should work well enough for linear axes too.
%
% This is useful when stacking panels on top of/next to each other without any
% space in between. It is however a crude method that assumes that MATLAB will
% not add ticks in certain ways when changing the plot range.
%
% Intended to simplify/replace sections "Remove overlapping Tics" in
% quicklooks_24_6_2_h() and quicklooks_7days().
%
% NOTE: Might not work for tick limits at zero. Is a special case.
%
%
% ARGUMENTS
% =========
% tickLimits
%       Length 2 vector. Min & max value for ticks in plot on one axis.
% dataLimits
%       Length 2 vector. Min & max value for data in plot on one axis.
%
%
% RETURN VALUES
% =============
% plotLimits
%       Length 2 row vector. Min & max value for plotted range in plot.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function plotLimits = ensure_data_tick_margins(tickLimits, dataLimits)
    % PROPOSAL: Do not assume that data values use a logarithmic axis.
    %   Ex: Nonweekly plots, panel 2 = density is linear.
    %   PRO: Should work better for limit=zero.
    %   --
    %   PROPOSAL: Argument for linear/log.
    %       CON: Has to hardcode lin/log for every panel.
    %           CON: Caller can read it from axis properties.
    %       PROPOSAL: Linear: Add margins which are a multiple of the range of data
    %                 (max minus min).
    %       PROPOSAL: Assume logarithmic, except when one tick=0 and use
    %                 assumption only for padding att tick=0. -- IMPLEMENTED

    % NOTE: 2022-03-22, 24h plot, panel 8/E_SRF has bad y margins for
    % "e723101f Erik P G Johansson (2023-05-10 18:10:25 +0200) SolO QLI:
    % Aesthetics-fix: Panel 2, left y label: Constant position"
    % (some data is outside the boundaries). This was later fixed in
    % "c1d06302 JordiBoldu (2023-05-11 14:09:27 +0200) Plot fixes"
    % on one branch in a code section for which a pre-bugfix version was
    % replaced by this code on another branch in parallel, without the bugfix.
    % Need to check that this code fixes the same bug eventually.
    % /EJ 2023-05-11

    C_LINEAR_MARGIN = 0.1;

    assert(length(tickLimits) == 2)
    assert(length(dataLimits) == 2)

    tickMin = tickLimits(1);
    tickMax = tickLimits(2);
    dataMin = dataLimits(1);
    dataMax = dataLimits(2);

    assert(tickMin <= tickMax)
    assert(dataMin <= dataMax)

    linearMargin = (dataMax - dataMin) * C_LINEAR_MARGIN;

    plotMax =  ensure_high_data_tick_margin( tickMax,  dataMax, linearMargin);
    plotMin = -ensure_high_data_tick_margin(-tickMin, -dataMin, linearMargin);

    plotLimits = [plotMin; plotMax];
end



function plotMax = ensure_high_data_tick_margin(tickMax, dataMax, linearMargin)
    C_DIMINISH = 0.9;
    C_MAGNIFY  = 1.1;

    assert(linearMargin >= 0)

    if tickMax > 0
        minPlotMax = C_MAGNIFY  * tickMax;
    elseif tickMax < 0
        minPlotMax = C_DIMINISH * tickMax;
    else
        % CASE: tickMax == 0

        % minPlotMax = dataMax;    % NOTE: Not tickMax. Bad?
        minPlotMax = linearMargin;
    end
    plotMax = max(dataMax, minPlotMax);
end



% Original code from quicklooks_24_6_2_h.m which this function largely replaces.
% NOTE: The old code should have a bug for the cases
% (1) mintick<=0 and minlim<0, and
% (2) maxtick<=0 and maxlim>0
% due to the use of abs().
% ------------------------------------------------------------------------------
% cax=h(iax);
% mintick = min(cax.YTick);
% maxtick = max(cax.YTick);
% minlim = cax.YLim(1);
% maxlim = cax.YLim(2);
%
% if maxtick>=0
%     if maxlim<1.1*maxtick
%         newmax = 1.1*maxtick;
%     else
%         newmax = maxlim;
%     end
% else
%     if abs(maxlim)>0.9*abs(maxtick)
%         newmax = 0.9*maxtick;
%     else
%         newmax = maxlim;
%     end
% end
%
% if mintick>0
%     if minlim>0.9*mintick
%         newmin = 0.9*mintick;
%     else
%         newmin = minlim;
%     end
% else
%     if abs(minlim)<1.1*abs(mintick)
%         newmin=1.1*mintick;
%     else
%         newmin=minlim;
%     end
% end
% cax.YLim=[newmin,newmax];
% ------------------------------------------------------------------------------
