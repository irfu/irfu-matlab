%
% Abstract superclass to bicas.proc.L1L2.cal.CalImpl to facilitate mocking in
% automated tests.
%
% NOTE: Might not implement all abstract methods for all subclass methods which
% actually used externally. This needs to be done if automated tests (which
% should use another subclass) require access to more methods.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) CalAbstract



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract, Access=public)



    bltsSamplesAVoltCa = calibrate_voltage_all(obj, ...
        dtSec, bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, ...
        zvcti, ufv)



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Abstract, Static)



    biasCurrentTm = calibrate_current_sampere_to_TM(currentSAmpere)



  end    % methods(Static)



end
