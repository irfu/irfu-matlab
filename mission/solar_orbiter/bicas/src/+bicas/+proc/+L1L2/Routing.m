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
    ssid

    % How the BLTS should be stored in datasets.
    sdid
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = Routing(ssid, sdid)
      assert(isa(ssid, 'uint8'))
      assert(isa(sdid, 'uint8'))

      % IMPLEMENTATION NOTE: Can not use these functions since
      % bicas.proc.sconst is initialized by calling this very constructor.
      %assert(bicas.proc.sconst.is_SSID(ssid) & isscalar(ssid))
      %assert(bicas.proc.sconst.is_SDID(sdid) & isscalar(sdid))

      obj.ssid = ssid;
      obj.sdid = sdid;
    end



  end    % methods(Access=public)



end
