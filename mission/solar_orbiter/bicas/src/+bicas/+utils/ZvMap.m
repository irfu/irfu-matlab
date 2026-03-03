%
% Map from key to ZV-like objects/arrays which all have the same number of
% rows/records. Keys do not have to be of the same type.
%
%
% DESIGN NOTES, RATIONALE
% =======================
% Is designed to potentially be *EXTENDED* to:
% * Collect all kinds of ZV-like variables (with the same number of records),
%   both FPAs and regular arrays.
% * Potentially be used recursively by itself emulating column size.
% * Potentially support indexing (subsasgn()) for individual values to be able
%   avoid separate operations for reading+modifying+overwriting value.
% * Support indexing (subsasgn(), subsref()) in the first dimension on
%   the entire object (i.e. for all keys, all constituent ZV-like values
%   simultaneously).
%
%
% NOTE: Currently (2025-07-09) only used for collecting data which is not obviously
% indexable over one dimension, in practice SDID.
%
% NOTE: Handle class.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ZvMap < handle
  % PROPOSAL: Automatic test code.
  % PROPOSAL: Rename "entries"--> ZVs.
  %
  % PROPOSAL: Use handle class which wraps a non-handle value.
  %   PRO: May be needed for implementing
  %    (1) subsasgn() for a specified key value (modify value without
  %        reading+modifying+overwriting with high performance.
  %    (2) subsref()/subsasgn() for all keys combined.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Dependent)
    nEntries
    keyCa
  end
  properties(SetAccess=private, GetAccess=public)
    nRecords
  end
  properties(SetAccess=private, GetAccess=private)
    Dict
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = ZvMap(nRecords)
      assert(isscalar(nRecords) & isnumeric(nRecords))

      obj.nRecords = nRecords;
      obj.Dict     = configureDictionary("cell", "cell");
    end



    % Asserts nonexisting key.
    function add(obj, key, value)
      % NOTE: Effectively excludes all char strings (row vectors) which are not
      % length one.
      assert(isscalar(key))
      assert(~obj.Dict.isKey({key}), "Key ""%s"" has already been added", key)
      assert(size(value, 1) == obj.nRecords)

      obj.Dict({key}) = {value};
    end



    % Asserts pre-existing key.
    function set(obj, key, value)
      % PROPOSAL: Check that the MC and size of the new value are identical to
      %           that of the old value.

      % NOTE: Effectively excludes all char strings (row vectors) which are not
      % length one.
      assert(isscalar(key))
      assert(obj.Dict.isKey({key}), "There is no key ""%s"".", key)
      assert(size(value, 1) == obj.nRecords)

      obj.Dict({key}) = {value};
    end



    function value = get(obj, key)
      value = obj.Dict({key});
      value = value{1};
    end



  end    % methods(Access=public)
  methods



    function nEntries = get.nEntries(obj)
      nEntries = obj.Dict.numEntries;
    end



    function keyCa = get.keyCa(obj)
      keyCa = obj.Dict.keys;
      keyCa = keyCa(:);
    end



  end



end
