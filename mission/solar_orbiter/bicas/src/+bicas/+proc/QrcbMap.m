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
  % PROPOSAL: Make add() permit pre-existing key.
  %   PRO: Can remove set().
  % PROPOSAL: Replace set() with method for OR:ing.
  % PROPOSAL: Replace add() and set() with method which does not require
  %           pre-existing key and which does OR:ing.
  %   PRO: Can simplify add_map().
  %
  % NOTE: add_map() is not entirely analogous to add(), but rather a combination
  % of add() and set().
  %
  % NOTE: Similar to bicas.utils.ZvMap.
  %   PROPOSAL: Implement using bicas.utils.ZvMap somehow?
  %     CON: ZVM does not enforce
  %          QRCID (scalar string)-->QRCB (logical; column array)
  %     CON: ZVM does not support add_map() with OR:ing for overlapping keys.



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
      assert(isstring(qrcid)   & isscalar(qrcid))
      assert(islogical(qrcbAr) & iscolumn(qrcbAr), "qrcbAr is not a logical column.")

      assert(~obj.Map.isKey(qrcid))    % Require not overwriting key.
      assert(size(qrcbAr, 1) == obj.nRecords)

      obj.Map(qrcid) = qrcbAr;
    end



    % Overwrite existing key-value pair. Must pre-exist.
    function set(obj, qrcid, qrcbAr)
      assert(isstring(qrcid)   & isscalar(qrcid))
      assert(islogical(qrcbAr) & iscolumn(qrcbAr))

      assert(obj.Map.isKey(qrcid))   % Require overwriting key.
      assert(size(qrcbAr, 1) == obj.nRecords)

      obj.Map(qrcid) = qrcbAr;
    end



    function remove(obj, qrcid)
      assert(isstring(qrcid) & isscalar(qrcid))
      assert(obj.Map.isKey(qrcid))

      obj.Map.remove(qrcid);
    end



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



    % Add one QRCB map to the current one. Non-overlapping keys are just added.
    % Overlapping keys have their values OR'ed together.
    %
    function add_map(obj, QrcbMap)
      assert(isa(QrcbMap, "bicas.proc.QrcbMap"))
      assert(obj.nRecords == QrcbMap.nRecords)

      newKeyCa = QrcbMap.Map.keys;
      for i = 1:numel(newKeyCa)
        key = newKeyCa{i};

        if obj.Map.isKey(key)
          % Merge key-value pairs.
          obj.Map(key) = obj.Map(key) | QrcbMap.Map(key);
        else
          % Add key-value pair.
          obj.Map(key) = QrcbMap.Map(key);
        end
      end
    end



  end    % methods(Access=public)



end
