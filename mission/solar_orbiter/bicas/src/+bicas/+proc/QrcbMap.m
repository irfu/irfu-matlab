%
% Class for map QRCID-->QRCB array.
%
% NOTE: Handle class due to internally using a containers.Map.
%
% IMPLEMENTATION NOTE: Seems to support equality, but unsure why. There are
% automated tests for equality.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcbMap < handle
  % PROPOSAL: Better name? Abbreviation?
  %     ~QRCID, ~QRCB, ~map
  %     ~arrays
  %     QrcbMap
  %       QRCBM=
  %     QrcBitMap
  %       CON: Sounds too much like bitmap (image).
  %     QrcbArraysMap
  %       QAM=
  %
  % PROPOSAL: Method which does not require pre-existing key and which does
  %           OR:ing.
  %   PRO: Can simplify union().



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Dependent)
    qrcidAr
  end
  properties(SetAccess=immutable)
    nRecords
  end    % properties(SetAccess=immutable)
  properties(SetAccess=immutable, GetAccess=private)
    Map
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods



    function qrcidAr = get.qrcidAr(obj)
      % IMPLEMENTATION NOTE: Sort to gain deterministic result which is good for
      % testing.
      qrcidAr = sort(string(obj.Map.keys'));
    end



  end
  methods(Access=public)



    function obj = QrcbMap(nRecords)
      assert(isscalar(nRecords) & isnumeric(nRecords) & (nRecords >= 0))

      obj.nRecords = nRecords;
      obj.Map      = containers.Map("KeyType", "char", "ValueType", "any");
    end



    % Add key-value pair. Must not pre-exist.
    function add(obj, qrcid, qrcbAr)
      assert(isstring(qrcid)   & isscalar(qrcid), "qrcid is not a scalar string.")
      assert(islogical(qrcbAr) & iscolumn(qrcbAr), "qrcbAr is not a logical column.")

      % Require not overwriting key.
      assert(~obj.Map.isKey(qrcid), "Object already has a key qrcid=""%s"".", qrcid)
      assert(size(qrcbAr, 1) == obj.nRecords)

      obj.Map(qrcid) = qrcbAr;
    end



    % Add QRCIDs with QRCB=false. useful for initialization.
    function add_false(obj, qrcidAr)
      assert(iscolumn(qrcidAr))

      qrcbAr = false([obj.nRecords, 1]);
      for qrcid = qrcidAr'
        obj.add(qrcid, qrcbAr)
      end
    end



    % Overwrite existing key-value pair. Must pre-exist.
    function set(obj, qrcid, qrcbAr)
      assert(isstring(qrcid)   & isscalar(qrcid))
      assert(islogical(qrcbAr) & iscolumn(qrcbAr))

      assert(obj.Map.isKey(qrcid))   % Require overwriting key.
      assert(size(qrcbAr, 1) == obj.nRecords)

      obj.Map(qrcid) = qrcbAr;
    end



    % NOTE: Method might be unused.
    % function remove(obj, qrcid)
    %   assert(isstring(qrcid) & isscalar(qrcid))
    %   assert(obj.Map.isKey(qrcid))
    %
    %   obj.Map.remove(qrcid);
    % end



    function qrcbAr = get(obj, qrcid)
      assert(isstring(qrcid) & isscalar(qrcid))
      if ~obj.has_QRCID(qrcid)
        error("Object does not contain qrcid=""%s"".", qrcid)
      end
      qrcbAr = obj.Map(qrcid);
    end



    function b = has_QRCID(obj, qrcid)
      % PROPOSAL: Better name.
      %   NOTE: Cf isKey().
      %   ~contains, has, is key
      %   ~QRCID
      assert(isstring(qrcid) & isscalar(qrcid))

      b = obj.Map.isKey(qrcid);
    end



    % Add one QRCBM to the current one. Non-overlapping keys are just added.
    % Overlapping keys have their values OR'ed together.
    %
    function union(obj, Qrcbm)
      assert(isa(Qrcbm, "bicas.proc.QrcbMap"))
      assert(obj.nRecords == Qrcbm.nRecords)

      newKeyCa = Qrcbm.Map.keys;
      for i = 1:numel(newKeyCa)
        key = newKeyCa{i};

        if obj.Map.isKey(key)
          % Merge key-value pairs.
          obj.Map(key) = obj.Map(key) | Qrcbm.Map(key);
        else
          % Add key-value pair.
          obj.Map(key) = Qrcbm.Map(key);
        end
      end
    end



    % Method for creating a very simple plot of the content of the object in a
    % separate figure.
    % NOTE: ONLY INTENDED FOR DEBUGGING!
    %
    function create_debug_figure(obj, tt2000Ar, figName)
      % PROPOSAL: Move to separate file for debug plot code.

      bicas.utils.assert_ZV_Epoch(tt2000Ar)
      assert(isstring(figName))
      assert(numel(tt2000Ar) == obj.nRecords)

      figure('WindowState', 'maximized', "Name", figName);
      tiledlayout(numel(obj.qrcidAr'), 1, "TileSpacing", "compact", "Padding", "none");
      for i = 1:numel(obj.qrcidAr')
        qrcid = obj.qrcidAr(i);
        h = nexttile;
        plot(tt2000Ar, obj.get(qrcid), '.-');
        h.YLim = [-0.05, 1.05];
        grid on

        legend(irf.graph.escape_str(qrcid))
      end
    end



  end    % methods(Access=public)



end
