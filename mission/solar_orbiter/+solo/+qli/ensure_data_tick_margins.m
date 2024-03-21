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
% and after the call to the function. This can be achieved using e.g. command
%   >> hAxes.YTickMode = 'manual'
% before calling the function. It is the caller's responsibility to make sure
% that there are enough ticks in the data range, e.g. by setting
%   >> hAxes.YLimMode  = 'auto'
%   >> hAxes.YTickMode = 'auto'
% first, and afterwards setting as mentioned above.
%   >> hAxes.YTickMode = 'manual'
%
%
% NOTE: This function (intentionally) does not operate on (read from or write
% to) any graphical objects. It only derives numerical values from other
% numerical values. This is better for e.g. automated tests and modularization.
%
%
% ARGUMENTS
% =========
% ticks
%       Vector (or empty) with tick values on the relevant axis.
%       Values may be inside and/or outside dataLimits.
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
%       max value for the displayed range on one axis in a plot.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function plotLimits = ensure_data_tick_margins(ticks, dataLimits, scale)
% PROPOSAL: Use terminology/naming similar to property names:
%           plotLimits : X/Y/ZLim
%           scale      : X/Y/ZScale
%
% PROPOSAL: Arguments for internal constants, C_LINEAR_MARGIN etc.
%
% NOTE: 2022-03-22, 24h plot, panel 8/E_SRF has bad y margins for
% "e723101f Erik P G Johansson (2023-05-10 18:10:25 +0200) SolO QLI:
% Aesthetics-fix: Panel 2, left y label: Constant position"
% (some data is outside the boundaries). This was later fixed in
% "c1d06302 JordiBoldu (2023-05-11 14:09:27 +0200) Plot fixes"
% on one branch in a code section for which a pre-bugfix version was
% replaced by this code on another branch in parallel, without the bugfix.
% Need to check that this code fixes the same bug eventually.
% /EJ 2023-05-11
%

%C_LINEAR_MARGIN = 0.1;
C_LINEAR_MARGIN = 0.05;

assert((isvector(ticks) || isempty(ticks)) && issorted(ticks, 'ascend'))
assert(length(dataLimits) == 2)

dataMin = dataLimits(1);
dataMax = dataLimits(2);

assert(dataMin <= dataMax)

linearMargin = (dataMax - dataMin) * C_LINEAR_MARGIN;
% IMPLEMENTATION NOTE: Does not want to derive linearMargin from ticks
% since:
% (1) there might be zero ticks,
% (2) ticks may be uncorrelated with dataMin/dataMax.

plotMax =  ensure_max_margin( ticks,  dataMax, scale, linearMargin);
plotMin = -ensure_max_margin(-ticks, -dataMin, scale, linearMargin);

plotLimits = [plotMin; plotMax];
end



% Find plotMax (only).
function plotMax = ensure_max_margin(ticks, dataMax, scale, linearMargin)
% PROPOSAL: Better way of deriving value "just below" higherBoundary.
%   PROPOSAL: Different methods depending on "scale".
%   PROPOSAL: higherBoundary-eps(higherBoundary)*C
% PROPOSAL: Special case for zero ticks.
%   PRO: Could simplify algorithm.

EPSILON    = 1e-10;
C_DIMINISH = 0.9;
C_MAGNIFY  = 1.1;
%C_DIMINISH = 0.8;
%C_MAGNIFY  = 1.25;

% Use ticks to define "boundaries" which additionally have values at the
% infinities. This ensures that there are always lower and higher boundaries
% which makes the rest of the algorithm simpler.
boundaries = [-inf; ticks(:); inf];

% Find nearest lower/equal boundary. Should always exist though can be -Inf.
lowerEqBoundary = max(boundaries(boundaries <= dataMax));
% Find nearest higher      boundary. Should always exist though can be +Inf.
higherBoundary  = min(boundaries(boundaries >  dataMax));

%===================================
% Adjust lowerEqBoundary for margin
%===================================
if strcmp(scale, 'linear')

  lowerEqBoundary = lowerEqBoundary + linearMargin;

elseif strcmp(scale, 'log')

  if lowerEqBoundary > 0
    lowerEqBoundary = lowerEqBoundary * C_MAGNIFY;
  elseif lowerEqBoundary < 0
    lowerEqBoundary = lowerEqBoundary * C_DIMINISH;
  else
    error('Illegal tickLowerEq=%d', lowerEqBoundary)
  end

else
  error('Illegal argument scale="%s"', scale)
end

assert(isscalar(lowerEqBoundary))
assert(isscalar(higherBoundary ))

% Use dataMax, but increase value if too low, i.e. to get a margin relative
% to the nearest lower tick (without considering higher ticks).
% NOTE: lowerEqBoundary may be -Inf.
plotMax = max(dataMax, lowerEqBoundary);

% If value is higher than the next higher tick, then decrease value to "just
% below" the next higher tick.
% NOTE: higherBoundary may be +Inf.
plotMax = min(plotMax, higherBoundary-EPSILON);
end



% function modify_by_factor(x, increase, f)
%     assert(f > 1)
%     assert(islogical(increase))
%
%     if ~increase
%         f = 1/f;
%     end
%
%     if x > 0
%         x = x * f;
%     elseif x < 0
%         x = x / f;
%     else
%         error('Illegal =%d', x)
%     end
%
% end



% function plotMax = ensure_lin_max_margin(tickMax, dataMax, linearMargin)
%     plotMax = max(dataMax, tickMax + linearMargin);
% end
%
%
%
% % NOTE: Handles negative values for historical reasons.
% function plotMax = ensure_log_max_margin(tickMax, dataMax, linearMargin)
%     C_DIMINISH = 0.9;
%     C_MAGNIFY  = 1.1;
%     %C_DIMINISH = 0.8;
%     %C_MAGNIFY  = 1.25;
%
%     assert(linearMargin >= 0)
%
%     if tickMax > 0
%         minPlotMax = C_MAGNIFY  * tickMax;
%     elseif tickMax < 0
%         minPlotMax = C_DIMINISH * tickMax;
%     else
%         % CASE: tickMax == 0
%
%         % minPlotMax = dataMax;    % NOTE: Not tickMax. Bad?
%         minPlotMax = linearMargin;
%     end
%     plotMax = max(dataMax, minPlotMax);
% end



% Original code from generate_quicklooks_24h_6h_2h.m which this function largely replaces.
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
