%
% Class for miscellaneous GA-related functions
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ga
  % PROPOSAL: Move/rename to bicas.ga.utils.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Given a value, if it is identical to any one in a specified list (column
    % cell array) of values, replace it with a specified value.
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function x = normalize_value(x, beforeValuesCa, afterValue)
      % PROPOSAL: Rename to not clash with other normalization.
      %   Ex: Adding "none" to empty list of strings.
      %   PROPOSAL: normalize_value
      %
      % PROPOSAL: Convert to generic function.
      % NOTE: Compare irf.utils.translate().
      %
      % PROBLEM: Function is used partly for normalizing CDF input. Can not easily
      % distinguish between non-compliant CDFs (and give warning/error) and
      % normalization.

      assert(...
        iscell(beforeValuesCa) & iscolumn(beforeValuesCa), ...
        'beforeValuesCa is not a column cell array.')

      for i = 1:numel(beforeValuesCa)
        beforeValue = beforeValuesCa{i};

        if isequaln(x, beforeValue)
          x = afterValue;
          return
        end
      end

    end



    % If a cell column array is empty, then return a one-element cell array
    % with the specified value. Intended for complying with metadata standards
    % which require replacing empty GAs with a one-entry GA with a string
    % constant.
    function gaCa = normalize_empty_column_array(gaCa, emptyValueStr)
      % TODO-DEC: "normalize" is the wrong term?

      assert(iscolumn(gaCa) & iscell(gaCa))
      assert(ischar(emptyValueStr))

      if isempty(gaCa)
        gaCa = {emptyValueStr};
      end
    end



  end    % methods(Static)



end
