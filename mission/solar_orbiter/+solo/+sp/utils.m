%
% Collection of trivial ~utility functions shared by summary plots.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-03-15
%
classdef utils
  % PROPOSAL: Automatic test code.

  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Merge two ~ZVs, with identical timestamps, but with non-overlapping
    % regions of non-NaN indices.
    function zv = merge_zvs(zv1, zv2)
      irf.assert.sizes(...
        zv1, [-1, -2], ...
        zv2, [-1, -2])

      b1 = ~isnan(zv1);
      b2 = ~isnan(zv2);

      % ASSERTION: Non-NaN components do not overlap, i.e. it is
      % unambiguous which non-NaN values to use.
      assert( all(~(b1 & b2), 'all'))

      zv     = zv1;
      zv(b2) = zv2(b2);

    end



  end

end
