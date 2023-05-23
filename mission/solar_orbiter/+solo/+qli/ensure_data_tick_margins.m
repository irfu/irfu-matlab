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
    %       PROPOSAL: Linear: Add margins which are a multiple of the range of data
    %                 (max minus min).
    %       PROPOSAL: Assume logarithmic, except when one tick=0 and use
    %                 assumption only for padding att tick=0.
    %
    % PROPOSAL: Refactor to use one function for both min and max.
    %   PRO: Symmetry.
    %   PRO: Avoids duplication of code.
    % PROPOSAL: Constants for 0.9 and 1.1.

    % NOTE: 2022-03-22, 24h plot, panel 8/E_SRF has bad y margins for
    % "e723101f Erik P G Johansson (2023-05-10 18:10:25 +0200) SolO QLI:
    % Aesthetics-fix: Panel 2, left y label: Constant position"
    % (some data is outside the boundaries). This was later fixed in
    % "c1d06302 JordiBoldu (2023-05-11 14:09:27 +0200) Plot fixes"
    % on one branch in a code section for which a pre-bugfix version was
    % replaced by this code on another branch in parallel, without the bugfix.
    % Need to check that this code fixes the same bug eventually.
    % /EJ 2023-05-11

    assert(length(tickLimits) == 2)
    assert(length(dataLimits) == 2)

    tickMin = tickLimits(1);
    tickMax = tickLimits(2);
    dataMin = dataLimits(1);
    dataMax = dataLimits(2);

    assert(tickMin <= tickMax)
    assert(dataMin <= dataMax)

    % linearMargin = (dataMax - dataMin) * 0.1;

    if tickMax>0
        minDataMax = 1.1 * tickMax;
    elseif tickMax<0
        minDataMax = 0.9 * tickMax;
    else
        minDataMax = dataMax;   % NOTE: Not tickMax.
    end
    plotMax = max(dataMax, minDataMax);

    if tickMin>0
        maxDataMin = 0.9 * tickMin;
    elseif tickMin<0
        maxDataMin = 1.1 * tickMin;
    else
        maxDataMin = dataMin;   % NOTE: Not tickMin.
    end
    plotMin = min(dataMin, maxDataMin);

    plotLimits = [plotMin; plotMax];
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
