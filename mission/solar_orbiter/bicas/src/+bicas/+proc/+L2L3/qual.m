%
% Collection of code relating to quality ZVs for L2 to L3 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual
  % PROPOSAL: Automatic test code.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function [QUALITY_FLAG, L3_QUALITY_BITMASK] = ...
        get_quality_ZVs_density(badDensityQrcbAr)

      % IMPLEMENTATION NOTE: Function is (as of 2023-12-18) in principle more
      % complicated than necessary w.r.t. L3_QUALITY_BITMASK but the
      % architecture is chosen to be analogue with
      % bicas.proc.L1L2.qual.get_quality_ZVs() so that it can be expanded in a
      % similar way if needed.

      assert(islogical(badDensityQrcbAr))

      nRecords = size(badDensityQrcbAr, 1);

      QrcbMap = bicas.proc.QrcbMap(nRecords);
      QrcbMap.add("BAD_DENSITY", badDensityQrcbAr);

      [QUALITY_FLAG, L3_QUALITY_BITMASK] = ...
        bicas.proc.qual.QRCB_arrays_to_quality_ZVs(...
        nRecords, QrcbMap, bicas.const.QRCS_L3_DENSITY_MAP);
    end



  end    % methods(Static)



end
