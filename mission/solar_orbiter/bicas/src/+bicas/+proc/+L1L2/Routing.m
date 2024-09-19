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



end
