%
% Immutable class which instances represent the destination of a signal, i.e.
% where the data should ultimately be stored in the output datasets, i.e.
% either:
% (1) an ASR (to determine which zVariable to store it in in the output
%     dataset), or
% (2) "nowhere" (when BDM is unknown).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SignalDestinationId
  % PROPOSAL: Better name.
  %   NOTE: Abbreviation SDID = Intended future abbreviation for "dataset
  %         ID"/DSI.
  %   --
  %   (output) dataset, zVariables
  %   destination
  %   ZDID = zVariable Destination ID
  %     CON: Does not contain "signal".
  %   DDID = Dataset Destination ID
  %     CON: Does not contain "signal".
  %   SZID = Signal zVariable ID
  %   PROPOSAL: Redefine to represent where in dataset to store signal.
  %
  % PROPOSAL: Re-define to include SSID.
  %   CON-PROPOSAL: Separate class for representing both SDID and SSID.
  %     NOTE: Resembles bicas.proc.L1L2.Routing.
  %   NOTE: Still needs SSID to represent where to store data.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=public)
    % ASID object or empty.
    asid

    % Whether destination is "NOWHERE", i.e. the signal does not have a
    % destination and should be ignored.
    isNowhere
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Constructor
    function obj = SignalDestinationId(value)
      assert(isscalar(value))

      if isequal(value, "NOWHERE")
        obj.asid      = [];
        obj.isNowhere = true;
      elseif bicas.sconst.is_ASID(value)
        obj.asid      = value;
        obj.isNowhere = false;
      else
        error('BICAS:Assertion:IllegalArgument', 'Illegal argument.')
      end
    end



  end    % methods(Access=public)



end
