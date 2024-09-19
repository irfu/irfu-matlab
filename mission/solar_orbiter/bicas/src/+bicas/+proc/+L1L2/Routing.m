%
% Immutable class that represents one particular routing of signals for one BLTS
% (but without specifying which BLTS). Given a particular demux mode and DLR
% setting, specify:
% (1) where the physical signal comes from (which among other things determines
%     how it should be calibrated), and
% (2) how it should be stored in the dataset object (if at all).
%
% Immutable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Routing
  % PROPOSAL: Better name
  %   PROPOSAL: Official abbreviation.
  %     "Routing" (!)
  %   PROBLEM: Class is not really related to the demultiplexer as the name
  %     "Routing" may imply. It really describes how to store (or not store)
  %     signals in datasets. Demuxer routing really refers to pairs (iBlts,
  %     SSID).
  %   --
  %   source, destination
  %   input, output
  %   demultiplexer, demuxer, demux
  %   mux(?)
  %   dataset, zVariables
  %   ~dataset representation
  %   --
  %   DemuxerRouting
  %   DemultiplexerRouting
  %   Source(To)DestinationRouting
  %       PRO: Can use abbreviation SDR, STDR.
  %   DSRP : DatasetRepresentation
  %     "DSR" already used abbreviation.
  %   ZVR = zVariable Representation.



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
    % Where the BLTS ultimately comes from.
    Ssid

    % How the BLTS should be stored in datasets.
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

      % Set SSID
      obj.Ssid = Ssid;

      % Set SDID
      switch numel(varargin)
        case 0
          assert(Ssid.is_ASR(), 'Can not use first argument to derive SDID.')
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



    function C = init_const()
      % PROPOSAL: Distinguish between different "channels" for 2.5V Ref
      %           and GND in the source (SSID).

      SSID = bicas.proc.L1L2.SignalSourceId.C;
      SDID = bicas.proc.L1L2.SignalDestinationId.C;

      C = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants(...
        @(Asid) (bicas.proc.L1L2.Routing(...
        bicas.proc.L1L2.SignalSourceId(Asid))));

      C.REF25V_TO_DC_V1    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V1);
      C.REF25V_TO_DC_V2    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V2);
      C.REF25V_TO_DC_V3    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V3);
      C.GND_TO_DC_V1       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V1);
      C.GND_TO_DC_V2       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V2);
      C.GND_TO_DC_V3       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V3);
      C.UNKNOWN_TO_NOWHERE = bicas.proc.L1L2.Routing(SSID.UNKNOWN, SDID.NOWHERE);
    end



  end



end
