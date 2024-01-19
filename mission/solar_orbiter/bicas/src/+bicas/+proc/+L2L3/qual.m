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



    function [QUALITY_FLAG, L3_QUALITY_BITMASK] = get_quality_ZVs_density(isBadDensityFpa)
      % IMPLEMENTATION NOTE: Function is (as of 2023-12-18) in principle
      % more complicated than necessary w.r.t. L3_QUALITY_BITMASK but the
      % architecture is chosen to (1) be analogue with
      % bicas.proc.L1L2.get_quality_ZVs().

      QrcFlagsMap = containers.Map();
      QrcFlagsMap(bicas.const.QRCID.BAD_DENSITY) = isBadDensityFpa;

      [QUALITY_FLAG, L3_QUALITY_BITMASK] = ...
        bicas.proc.qual.QRC_flag_arrays_to_quality_ZVs(...
        size(isBadDensityFpa, 1), QrcFlagsMap, bicas.const.QRC_SETTINGS_L3_DENSITY);
    end



  end    % methods(Static)



end
