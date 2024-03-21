%
% Immutable class that represents a particular routing of signals via a
% particular BLTS given a particular demux mode and DLR setting:
% (1) where a physical signal comes from (which relates to how
%     it should be calibrated), and
% (2) how it should be stored in the dataset object (if at all).
%
% Immutable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Routing   % < handle
  % PROPOSAL: Better name
  %   source, destination
  %   demultiplexer
  %   DemultiplexerRouting
  %   Source(To)DestinationRouting
  %       PRO: Can use abbreviation SDR.



  %###################
  %###################
  % STATIC PROPERTIES
  %###################
  %###################
  properties(GetAccess=public, Constant)
    C = bicas.proc.L1L2.Routing.init_const();
  end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Where the physical signal in the BLTS ultimately comes from. This is
    % used to determine how the signal should be calibrated.
    Ssid

    % How the BLTS should be stored in the datasets.
    Sdid
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)


    % ARGUMENTS
    % =========
    % Syntax 1: Ssid
    %       Reuse Ssid.Asid for creating a corresponding SDID.
    % Syntax 2: Ssid, Sdid
    function obj = Routing(Ssid, varargin)
      assert(isa(Ssid, 'bicas.proc.L1L2.SignalSourceId'))
      obj.Ssid = Ssid;

      switch numel(varargin)
        case 0
          assert(Ssid.is_ASR())
          Sdid = bicas.proc.L1L2.SignalDestinationId(Ssid.Asid);
        case 1
          Sdid = varargin{1};
        otherwise
          error('BICAS:Assertion:IllegalArgument', ...
            'Illegal number of extra arguments.')
      end
      assert(isa(Sdid, 'bicas.proc.L1L2.SignalDestinationId'))
      obj.Sdid = Sdid;

    end



  end    % methods(Access=public)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Access=private, Static)



    function R = init_const()
      % PROPOSAL: Distinguish between different "channels" for 2.5V Ref
      %           and GND in the source (SSID).

      SSID = bicas.proc.L1L2.SignalSourceId.C;
      SDID = bicas.proc.L1L2.SignalDestinationId.C;
      R = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants(...
        @(Asid) (bicas.proc.L1L2.Routing(...
        bicas.proc.L1L2.SignalSourceId(Asid))));

      R.REF25V_TO_DC_V1    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V1);
      R.REF25V_TO_DC_V2    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V2);
      R.REF25V_TO_DC_V3    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V3);
      R.GND_TO_DC_V1       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V1);
      R.GND_TO_DC_V2       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V2);
      R.GND_TO_DC_V3       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V3);
      R.UNKNOWN_TO_NOWHERE = bicas.proc.L1L2.Routing(SSID.UNKNOWN, SDID.NOWHERE);
    end



  end



end
