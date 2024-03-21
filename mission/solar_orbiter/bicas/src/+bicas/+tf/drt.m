%
% Class for detrending and retrending a signal e.g. before & after the
% application of a TF to it.
%
% DRT = De/Re-Trending
%
%
% IMPLEMENTATION NOTES
% ====================
% Reasons that de- & re-trending code exists as separate code:
%   ** reusability
%   ** automatically testable
%   ** isolation/encapsulation
%   ** collect (and reuse) the assertions on settings
% --
% Code could also be implemented as a function which takes a function handle to
% the function it effectively wraps. That implementation has not been chosen
% due to:
%   ** general avoidance of function handles
%   ** risk of complexity if combining with other functionality to modify
%      signals before & after application of TF
%   ** the current implementation is an experiment to see if the model is
%      useful.
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% det : detrending
% ret : retrending
% IMPLEMENTATION NOTE: Above is used to make the words harder to confuse with
% each other.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-12
%
classdef drt < handle
  %
  % PROPOSAL: Better name.
  %   PRO: "drt" violates current naming conventions.
  %   Deretrending
  %   Deretr
  %   Detrending
  %   DetrendingRetrending



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Access=private)

    detDegreeOf
    detEnabled   % In principle unnecessary, but clarifies code.
    retEnabled
    yTrend1

    % State to verify that methods are executed in correct order.
    state = 0;
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = drt(detDegreeOf, retEnabled)
      % ASSERTIONS
      assert(obj.state == 0)
      assert(isscalar(detDegreeOf) & isnumeric(detDegreeOf))
      assert(isscalar(retEnabled ) & islogical(retEnabled ), ...
        'retEnabled is not scalar & logical.')

      obj.detDegreeOf = detDegreeOf;
      obj.detEnabled = (detDegreeOf >= 0);
      obj.retEnabled  = retEnabled;
      obj.yTrend1     = [];    % Set later.

      % ASSERTION: No RE-trending if no DE-trending.
      assert(obj.detEnabled || ~obj.retEnabled, ...
        'BICAS:Assertion:IllegalArgument', ...
        ['Illegal combination of', ...
        ' "detDegreeOf" and "retEnabled".'])

      obj.state = 1;
    end



    function y1b = detrend(obj, y1a)
      % ASSERTIONS
      assert(obj.state == 1)
      assert(iscolumn(y1a), 'Argument y1a is not a column vector.')

      if obj.detEnabled
        nSamples         = length(y1a);

        trendFitsCoeffs1 = polyfit((1:nSamples)', y1a, obj.detDegreeOf);
        % NOTE: Set instance variable.
        % NOTE: By necessity hardcodes the vector "type" (row/column).
        obj.yTrend1      = polyval(trendFitsCoeffs1, (1:nSamples)');

        y1b              = y1a - obj.yTrend1;
      else
        y1b = y1a;
      end

      obj.state = 2;
    end



    function y2a = retrend(obj, y2b, scaleFactor)
      % ASSERTIONS
      assert(obj.state == 2)
      assert(iscolumn(y2b), 'Argument y2b is not a column vector.')
      assert(isscalar(scaleFactor))

      if obj.retEnabled
        yTrend2 = obj.yTrend1 * scaleFactor;
        y2a     = y2b + yTrend2;
      else
        y2a     = y2b;
      end

      obj.state = 3;
    end



  end    % methods(Access=public)



end
