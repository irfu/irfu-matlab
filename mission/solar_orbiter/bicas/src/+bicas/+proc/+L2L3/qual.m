%
% Collection of code relating to quality ZVs for L2 to L3 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual   % < handle
  % PROPOSAL: Automatic test code.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function [QUALITY_FLAG, L3_QUALITY_BITMASK] = get_quality_ZVs_density(isBadDensityAr)
      % IMPLEMENTATION NOTE: Function is (as of 2023-12-18) in principle
      % more complicated than necessary w.r.t. L3_QUALITY_BITMASK but the
      % architecture is chosen to (1) be analogue with
      % bicas.proc.L1L2.qual.get_quality_ZVs().
      assert(islogical(isBadDensityAr))

      QrcFlagsMap = containers.Map();
      QrcFlagsMap(bicas.const.QRCID.BAD_DENSITY) = isBadDensityAr;

      [QUALITY_FLAG, L3_QUALITY_BITMASK] = ...
        bicas.proc.qual.QRC_flag_arrays_to_quality_ZVs(...
        size(isBadDensityAr, 1), QrcFlagsMap, bicas.const.QRCS_L3_MAP_DENSITY);
    end



  end    % methods(Static)



end
