%
% Class whose instances represent one loaded dataset (CDF file).
%
% NOTE: Not to be confused with bicas.swm.InputDataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-09
%
classdef InputDataset
  % PROPOSAL: Make all properties/instance variables immutable. Modify by
  %           creating modified copy of class instance.
  %
  % PROPOSAL: Use abbreviation for bicas.InputDataset and bicas.OutputDataset.
  %   NOTE: Should be consistent
  %     PROPOSAL: DS=DataSet
  %       PRO: Consistent with DSI=DataSet ID
  %     IDS, ODS
  %     IPDS, OPDS



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################

  %====================
  % Mutable properties
  %====================
  properties
    % IMPLEMENTATION NOTE: Variables here are not immutable since some
    % values are normalized during processing (2x6x L1/L1R-->L2).

    % Struct with zVariables stored as FPAs.
    ZvFpa
    % Struct with zVariables stored as plain arrays.
    Zv
  end

  %======================
  % Immutable properties
  %======================
  properties(SetAccess=immutable)

    % Struct with zVariable fill values.
    ZvFv

    % Struct with global attributes.
    Ga

    filePath

  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = InputDataset(ZvFpa, Zv, ZvFv, Ga, filePath)
      obj.ZvFpa    = ZvFpa;
      obj.Zv       = Zv;
      obj.ZvFv     = ZvFv;
      obj.Ga       = Ga;
      obj.filePath = filePath;
    end



  end    % methods(Access=public)



end
