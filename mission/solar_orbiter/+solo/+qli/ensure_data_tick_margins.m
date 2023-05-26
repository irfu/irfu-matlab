%
% For a given axis (X/Y/Z): given
%   (1) preexisting min/max ticks (min/max values at which there are ticks), and
%   (2) min/max *data* values,
% derive suggested plot min/max limits (e.g. YLim) which will ensure that there
% is a margin between the (old) ticks and the returned plot min/max limits. The
% returned min/max limits may be the same as the submitted ones.
%
% This is useful when stacking panels on top of/next to each other without any
% space in between. Tick labels may overlap in such cases. It is however a crude
% method that does not take font sizes into account and that may fail.
%
% The function assumes that the user ensures that the same ticks are used before
% and after the call to the function. This can be achieved using a command
% h.YTickMode = 'manual' before calling the function.
%
% NOTE: This function (intentionally) does not operate on (read from or write
% to) any graphical objects. It only derives numerical values from other
% numerical values. This is better for e.g. automated tests and modularization.
%
%
% ARGUMENTS
% =========
% tickLimits
%       Length-2 vector. Min & max value for ticks in plot on the relevant axis.
%       These will be inside the final returned plot limits.
% dataLimits
%       Length-2 vector. Min & max value for data in plot on the relevant axis.
% scale
%       String constant. 'linear' or 'log'.
%       NOTE: Same constants as in graphical object properties X/Y/ZScale.
%
%
% RETURN VALUES
% =============
% plotLimits
%       Length-2 row vector. Suggested values for property X/Y/ZLim, i.e. min &
%       max value for the range in plot.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function plotLimits = ensure_data_tick_margins(tickLimits, dataLimits, scale)
    % PROPOSAL: Do not assume that data values use a logarithmic axis.
    %   Ex: Nonweekly plots, panel 2 = density is linear.
    %   PRO: Should work better for limit=zero.
    %   --
    %   PROPOSAL: Argument for linear/log. -- IMPLEMENTED
    %       CON: Has to hardcode lin/log for every panel.
    %           CON: Caller can read it from axis properties.
    %       PROPOSAL: Linear: Add margins which are a multiple of the range of data
    %                 (max minus min).
    %       PROPOSAL: Assume logarithmic, except when one tick=0 and use
    %                 assumption only for padding att tick=0.
    %
    % PROPOSAL: Use terminology/naming similar to property names:
    %           plotLimits : X/Y/ZLim
    %           scale      : X/Y/ZScale
    %
    % PROPOSAL: Arguments for internal constants, C_LINEAR_MARGIN etc.

    % NOTE: 2022-03-22, 24h plot, panel 8/E_SRF has bad y margins for
    % "e723101f Erik P G Johansson (2023-05-10 18:10:25 +0200) SolO QLI:
    % Aesthetics-fix: Panel 2, left y label: Constant position"
    % (some data is outside the boundaries). This was later fixed in
    % "c1d06302 JordiBoldu (2023-05-11 14:09:27 +0200) Plot fixes"
    % on one branch in a code section for which a pre-bugfix version was
    % replaced by this code on another branch in parallel, without the bugfix.
    % Need to check that this code fixes the same bug eventually.
    % /EJ 2023-05-11

    % ~DESIGN BUG: Only has arguments for min/max ticks. However, there might be
    %              multiple ticks outside the data range.
    %   PROPOSAL: Argument for array of all ticks. Use the next larger/smaller
    %             tick and make sure that no tick outside that is used/visible.
    %       PROPOSAL: Return array of new ticks.
    %           PRO: Can remove ticks outside the nearest smaller/larger tick,
    %                but inside the tick+margin.

    C_LINEAR_MARGIN = 0.1;

    assert(length(tickLimits) == 2)
    assert(length(dataLimits) == 2)

    tickMin = tickLimits(1);
    tickMax = tickLimits(2);
    dataMin = dataLimits(1);
    dataMax = dataLimits(2);

    assert(tickMin <= tickMax)
    assert(dataMin <= dataMax)

    %linearMargin = (dataMax - dataMin) * C_LINEAR_MARGIN;
    linearMargin = (tickMax - tickMin) * C_LINEAR_MARGIN;

    if strcmp(scale, 'linear')
        plotMax =  ensure_lin_max_margin( tickMax,  dataMax, linearMargin);
        plotMin = -ensure_lin_max_margin(-tickMin, -dataMin, linearMargin);
    elseif strcmp(scale, 'log')
        plotMax =  ensure_log_max_margin( tickMax,  dataMax, linearMargin);
        plotMin = -ensure_log_max_margin(-tickMin, -dataMin, linearMargin);
    else
        error('Illegal argument scale="%s"', scale)
    end

    plotLimits = [plotMin; plotMax];
end



function plotMax = ensure_lin_max_margin(tickMax, dataMax, linearMargin)
    plotMax = max(dataMax, tickMax + linearMargin);
end



% NOTE: Handles negative values for historical reasons.
function plotMax = ensure_log_max_margin(tickMax, dataMax, linearMargin)
    C_DIMINISH = 0.9;
    C_MAGNIFY  = 1.1;
    %C_DIMINISH = 0.8;
    %C_MAGNIFY  = 1.25;

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
