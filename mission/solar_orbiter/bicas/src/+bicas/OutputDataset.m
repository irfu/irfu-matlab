%
% Represents a concrete output dataset produced by processing, primarily the
% content.
%
% NOTE: Not to be confused with bicas.swm.OutputDataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef OutputDataset
  % PROPOSAL: Automatic test code.
  %
  % PROPOSAL: Use abbrevatiation.
  %   See bicas.InputDataset comments.
  %
  % PROBLEM: Accepts RctdCa, despite that this format is a function of
  %          processing, and should not be useful outside of processing.
  %   PROPOSAL: Constructor argument for array of RCTDs. -- IMPLEMENTED
  %     CON: Multiple callers have to convert format.
  %   PROPOSAL: Constructor converts RctdCaMap to array of RCTDs.
  %
  % PROPOSAL: Argument and field for DSI. Replace OutputDatasetsMap with array.
  %   PRO: Should not need OutputDatasetsMap keys like "DSR_cdf",
  %        "EFIELD_DSR_cdf", "SCI_cdf" which make up its own namespace.
  %     NOTE: LfrSwmProcessing, TdsSwmProcessing are used for multiple SWMs with
  %           different output DSIs for the same OutputDatasetsMap key
  %           "SCI_cdf".
  %     NOTE: This is identical to bicas.swm.OutputDataset.pfoid.
  %   NOTE: Must match DSI in corresponding SWM.
  %
  % TODO-DEC: How handle modifications before writing dataset? How handle
  %           modifications to GAs? Should this class represent output dataset
  %           before/after deriving GAs?
  %   NOTE: Using master CDFs is also a kind of modification of the output
  %         dataset (conceptually).
  %   PROPOSAL: Modify (create new) OutputDataset from old one.
  % PROPOSAL: Use class as argument to bicas.write_dataset_CDF().
  %   CON: May contain unimportant information for actual file-writing.
  %     Ex: RctdCa
  %     Ex: BICAS-internal metadata
  %       Ex: DSI (not currently)
  %       Ex: File path (not currently)




  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    Zv
    Ga
    RctdCa
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = OutputDataset(Zv, Ga, RctdCa)
      assert(isstruct(Zv))
      assert(isstruct(Ga))
      assert(iscell(RctdCa) & iscolumn(RctdCa))

      obj.Zv = Zv;
      obj.Ga = Ga;
      obj.RctdCa = RctdCa;
    end



  end    % methods(Access=public)



end
