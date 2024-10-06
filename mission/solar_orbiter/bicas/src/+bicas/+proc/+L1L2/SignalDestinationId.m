%
% Immutable class which instances represent an SDID with metadata.
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
    % ASID or empty.
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
    function obj = SignalDestinationId(asidOrNowhere)
      assert(isscalar(asidOrNowhere))

      if isstring(asidOrNowhere) && isequal(asidOrNowhere, "NOWHERE")
        obj.asid      = [];
        obj.isNowhere = true;
      elseif isa(asidOrNowhere, 'uint8')
        obj.asid      = asidOrNowhere;
        obj.isNowhere = false;
      else
        error('BICAS:Assertion:IllegalArgument', 'Illegal argument.')
      end
    end



    function isAsr = is_ASR(obj)
      isAsr = ~isempty(obj.asid);
    end



  end    % methods(Access=public)



end
