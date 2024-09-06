%
% Basic map key-->value, where all values must have the same number of rows
% (size in first dimension). Values must be immutable.
%
% Intended for zVariables. Implementation may be upgraded to require FPAs.
% Some method names are chosen to be identical with dictionary.
%
%
% IMPLEMENTATION NOTE
% ===================
% bicas.utils.SameRowsMap.set_rows() can be slow *IF* the implementation stores
% data directly as (non-handle) values in containers.Map (old implementation)
% and presumably in dictionary, presumably since preallocation does not work.
% To avoid this, the implementation instead stores all values indirectly via
% handle class objects (bicas.utils.HandleWrapper), which (apparently) makes it
% possible to modify arrays without implicit copying by MATLAB, thus increasing
% performance. Does seem to work. Since the class's internal data structure
% uses handle objects, the class itself also has to be a handle class, to avoid
% that internal handle objects are shared between different instances of
% bicas.utils.SameRowsMap.
% --
% Enforces that the MATLAB class for keys is consistent, treating char
% array<>string as opposed to the dictionary implementation.
%
%
% RATIONALE
% =========
% Class was created for the purpose of managing sets of ZV-like arrays, where
% the first dimension in each array represent CDF records and is thus equal in
% size for all arrays. This could be useful for doing indexing operations and
% enforcing the same number of CDF records/rows. In practice, this seems to not
% be as useful as intended, and harder to implement than expected. The class is
% therefore not widely used. Still useful for maintaining a map ASR
% channel-->samples. /Erik P G Johansson 2023-10-02
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SameRowsMap < handle
  % PROPOSAL: Use syntactic sugar somehow.
  %   set, get, length, nRows, keys
  %
  %   PROPOSAL: Overload indexing: set/add
  %   PROPOSAL: Properties for nRows, keys
  %       CON: Unnecessary since they are read-only, and methods can be called
  %            without brackets.
  %       PRO: Values become part of the default human-readable text represention(?)
  %
  % PROPOSAL: Rename method length().
  %   PRO: Does not refer to length of object as in object emulating an array.
  %
  % TODO-DEC: How handle ~indexing (overload)? What should it be used for?
  %   PROBLEM: Would like to use indexing for (1) specifying variables, and
  %            (2) specifying rows.
  %       NOTE: Variable keys can both be strings and numbers, but never ranges or "end".
  %             Rows can be any 1D-numerical index.
  %   TODO-NI: How does this impact implementing SameSizeTypeMap using
  %            (composition) SameRowsMap?
  %   PROPOSAL: Use {} for of them and () for the other.
  %   PROPOSAL: Support both
  %       (1) (variableKey)
  %       (2) (:, <rowIndex>), including : and "end".
  %       PRO: Conceptually clean. Variables are another dimension.
  %       CON: Complicated to implement(?).
  %       NOTE: (:, <rowIndex>) must be supported for assignment.
  %
  % PROPOSAL: Make compatible with bicas.utils.FPArray.
  %   NOTE: Should already be, assuming that FPA implements subsref(),
  %         subsasgn() etc.
  %
  % NEED: Simultaneously store collections of same-rows ZVs and same-size&type
  %       ZVs and enforce assertions between them.
  %   PROPOSAL: Class(es) which support a hierarchy of collections.
  %       SameSizeType inside/under SameRows.
  %   PROPOSAL: Collections which are connected to each other (keep
  %       references).
  %
  % NOTE: Examples of usage for ZV collections:
  %   AsrMap (while building) (until bicas.utils.SameSizeTypeMap)
  %   AsrMap (final)          (until bicas.utils.SameSizeTypeMap)
  %   DemultiplexingCalibrationInput/DemultiplexingCalibrationOutput.Zv
  %   bicas.proc.L1L2.dc.calibrate_demux_voltages_subsequence(): Arguments which represent
  %       a subsequence (a shorter interval of CDF records).
  %
  % PROPOSAL: Fixed set of keys.
  % PROPOSAL: Immutable.
  %   PRO: Makes inheritance SameRowsMap-->SameSizeTypeMap natural.
  %   CON: Less convenient for building AsrMap which can be variable-number of
  %        variables.
  %       PROPOSAL: Build content using dictionary. Submit to constructor.
  %           CON: Can never extend to accept wider set of keys, e.g. objects.
  %           PRO: Easier to implement.
  %
  %
  %
  % TODO-DEC: How implement relationships between similar classes of variable
  %           collections?
  %   NEED: Shared test code.
  %   -
  %   PROPOSAL: Merge SameRowsMap and SameSizeTypeMap. Constructor
  %             argument determines which map value properties should be consistent.
  %       PROPOSAL: Submit function handles that convert variable to value
  %                 that should be consistent between variables.
  %           CON: Mixing consistent and inconsistent sizes does not fit with
  %                having one scheme for indexing.
  %       CON: Can not share all methods, since some methods only makes sense
  %            for one of the classes.
  %           Ex: SameSizeTypeMap:  size():  All dimensions.
  %           Ex: SameSizeTypeMap: ~class(): MATLAB class.
  %               CON: Should not be needed.
  %       NOTE: Might still want to set consistent values (e.g. nRows) in
  %             constructor to handle case of zero keys.
  %
  %   PROPOSAL: Specify how many (first) dimensions must be consistent.
  %       CON: Overkill for BICAS.
  %       PRO: Natural generalization.
  %       PROPOSAL: Specify max number of (non-size-one) dimensions.
  %       PROPOSAL: Constructor flag for whether to require consistent type.
  %           CON: Not compatible with using FPA.
  %               CON: No other scheme is either.
  %               CON-PROPOSAL: Argument FH for specifying non-size assertion(s)
  %                             on variables/objects.
  %                   CON: Can not reliably compare that quality (FHs) for different
  %                        instances.
  %           PRO: Can merge SameRowsMap and SameSizeTypeMap and covers more use cases.
  %           TODO-DEC: How compare maps? Include/exclude consistency flags?
  %           CON: Ugly to have different behaviour for different instances.
  %               CON: No more ugly than specifying number of consistent
  %                    dimensions.
  %               CON: Must specify extra assertions on arguments if wants to
  %                    test for this.
  %                   Ex: assert(Sstm.sameType & Sstm.consistentDims == Inf)
  %                   CON-PROPOSAL: Define subclasses for specific configurations.
  %   --
  %   PROPOSAL: SameSizeTypeMap is subclass of SameRowsMap.
  %       PRO: Can add methods size(), class() to subclass.
  %       CON: Subclass and superclass are identical w.r.t. reading, but not
  %            w.r.t. writing content. Subclass adds more assertions.
  %           PRO: assert(isa(x, 'SameRowsMap')) on an argument for reading
  %                makes sense, but not on an argument for writing to.
  %       PROPOSAL: Make immutable.
  %           CON: Can not implement set_rows()
  %                as method (without creating new object).
  %   PROPOSAL: SameSizeTypeMap, SameRowsMap are subclasses of the same
  %             abstract superclass.
  %       NOTE: Must be abstract superclass to avoid same superclass-subclass
  %             problems.
  %   PROPOSAL: SameSizeTypeMap composes/wraps SameRowsMap.
  %       CON: Must redefine all methods, just for wrapping.
  %   PROBLEM: How compare SameRowsMap with SameSizeTypeMap?
  %
  % PROPOSAL: Prefix class name with "H" to imply handle class.
  % PROPOSAL: Prefix all variables that point to handle objects "h".
  %   Ex: hSrm
  %   NOTE: Does not want to use capital "H" since could cause ambiguity w.r.t.
  %         abbreviations.
  %   CON: Can not make distinction between pointing to objects/structs and not.
  %       CON: Does not need to. Pointed-to objects are always objects.
  %
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TODO: Check that pre-allocation works when using subsasgn.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Access=private)
    Dict

    % IMPLEMENTATION NOTE: Can not use name "nRows" since it coincides with
    % method name.
    nRows2

    % MC for keys.
    % IMPLEMENTATION NOTE: Does not rely on the dictionary.types since
    % dictionary treats (1) char array == string, and (2) different number
    % types as equal.
    mcKeys
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % mcKeys
    %       MATLAB class for keys.
    % varargin
    %       initType == 'EMPTY':    Zero length.
    %       initType == 'CONSTANT':
    %           varargin{1} = Array of values
    %           varargin{2} = Cell array of keys which should have this
    %                         value.
    %
    % NOTE: To initialize with multiple keys with unique values, use both
    %       constructor and method "add".
    function obj = SameRowsMap(mcKeys, nRows, initType, varargin)
      assert(isnumeric(nRows) && nRows >= 0, 'nRows is not a non-negative number.')

      obj.nRows2 = nRows;
      obj.Dict   = configureDictionary(mcKeys, 'bicas.utils.HandleWrapper');
      obj.mcKeys = mcKeys;

      switch(initType)
        case 'EMPTY'
          % Zero initial keys and values.

          assert(numel(varargin) == 0)

        case 'CONSTANT'
          % Same (constant) initial value for all keys.

          assert(numel(varargin) == 2)
          value  = varargin{1};
          keysCa = varargin{2};
          % NOTE: Does not check for unique keys.
          for keyCa = keysCa(:)'
            bicas.utils.SameRowsMap.assert_legal_key(keyCa{1})
          end

          for keyCa = keysCa(:)'
            obj.add(keyCa{1}, value);
          end

        otherwise
          error('BICAS:Assertion:IllegalArgument', ...
            'Illegal argument initType="%s"', initType)
      end

    end



    % Number of variables inside the object. Unrelated to their size (e.g.
    % rows).
    function n = length(obj)
      n = obj.Dict.numEntries;
    end



    % NOTE: Method name chosen to be identical with dictionary.keys().
    function keysCa = keys(obj)
      keysCa = num2cell(obj.Dict.keys);
    end



    % Mostly for debugging.
    function valuesCa = values(obj)
      hwCa     = obj.Dict.values();
      valuesCa = arrayfun(@(x) (x.v), hwCa, 'UniformOutput', false);
      valuesCa = valuesCa(:);
    end



    % NOTE: Method name chosen to be identical with dictionary.isKey().
    % The name is therefore inconsistent with naming conventions.
    function isKey = isKey(obj, key)
      bicas.utils.SameRowsMap.assert_legal_key(key)

      isKey = obj.Dict.isKey(key);
    end



    % Add NEW key-value pair. Disallow overwriting.
    function add(obj, key, value)
      bicas.utils.SameRowsMap.assert_legal_key(key)
      assert(isa(key, obj.mcKeys))
      assert(~obj.Dict.isKey(key))
      assert(obj.nRows2 == size(value, 1), ...
        'The argument''s number of rows (%i) is not equal to the object''s number of rows (%i).', ...
        obj.nRows2, size(value, 1))

      obj.Dict(key) = bicas.utils.HandleWrapper(value);
    end



    % Use another SRM to overwrite selected rows in this SRM.
    %
    % IMPLEMENTATION NOTE: Method is important for speeding up LFR-SWF which
    % tends to be broken into subsequences of 1 record. This can be done by
    % pre-allocating and then overwriting parts.
    %
    %
    % ARGUMENTS
    % =========
    % Srm2
    %       bicas.utils.SameRowsMap. Must have the same set of keys and
    %       value types.
    % iRowsArray
    %       Column array. Same length as number of rows in Srm2 fields.
    %       Specifies the rows that shall be overwritten.
    %       NOTE: Can not use logical indexing.
    %
    %
    function set_rows(obj, Srm2, iRowsArray)
      % PROPOSAL: Support logical indexing.
      % PROPOSAL: Support using an 1-row SRM for overwriting N rows.
      % PROPOSAL: Implement method by overloading indexing notation:
      %           subsasgn.
      %   PRO: Using subsasgn() for assigning internal arrays might handle
      %        (1) logical indexing), (2) assigning multiple rows with
      %        1-row SRM.
      assert(isa(Srm2, 'bicas.utils.SameRowsMap'))
      assert(isnumeric(iRowsArray) && iscolumn(iRowsArray))
      assert(size(iRowsArray, 1) == Srm2.nRows)

      assert(bicas.utils.object_sets_isequaln(obj.keys, Srm2.keys))

      keysCa = obj.keys();
      for keyCa = keysCa(:)'
        key = keyCa{1};

        hw1 = obj.Dict(key);
        hw2 = Srm2.Dict(key);

        size1 = size(hw1.v);
        size2 = size(hw2.v);
        % IMPLEMENTATION NOTE: Using num2str(key) since it can handle
        % both strings and numbers.
        assert(isequal(size1(2:end), size2(2:end)    ), ...
          'Values for key="%s" have inconsistent sizes.', num2str(key))
        assert(isequal(class(hw1.v), class(hw2.v)), ...
          'Values for key="%s" have inconsistent MATLAB classes.', num2str(key))

        % IMPLEMENTATION NOTE: Unsure, but think that explicitly setting
        % second dimension to ":" makes the command handle any
        % dimensionalities (assuming that dimensionalities and sizes are
        % consistent).
        hw1.v(iRowsArray, :) = hw2.v(:, :);

        % IMPLEMENTATION NOTE: Does not need to set obj.Dict(key) since
        % using handle classes.
      end
    end



    % Syntactic sugar.
    % Support indexing with (): array = Srm(key)
    function varargout = subsref(obj, S)

      switch(S(1).type)
        case '()'
          assert(isscalar(S))
          assert(isscalar(S.subs), 'Illegal index. Must be exactly one argument.')

          if isnumeric(S.subs{1})
            assert(isscalar(S.subs{1}), 'Illegal index. Value must be scalar.')
          end

          % IMPLEMENTATION NOTE: Only intended for singular values,
          % whether strings or numbers. Should not support indices
          % like "1,2", colons, or "end".
          key = S.subs{1};
          assert(isa(key, obj.mcKeys))
          varargout = {obj.Dict(key).v};

        otherwise
          % CASE: {} or .

          % Fail for {}. % Call method/property.
          [varargout{1:nargout}] = builtin('subsref', obj, S);
      end
    end



    function nRows = nRows(obj)
      nRows = obj.nRows2;
    end



    function equals = eq(obj1, obj2)
      % IMPLEMENTATION NOTE: Must support the case of zero keys.

      assert(isa(obj2, 'bicas.utils.SameRowsMap'))

      if ~bicas.utils.object_sets_isequaln(obj1.keys, obj2.keys)
        equals = false;
        return
      elseif obj1.nRows ~= obj2.nRows
        equals = false;
        return
      elseif ~isequaln(obj1.mcKeys, obj2.mcKeys)
        equals = false;
        return
      end

      keysCa = obj1.Dict.keys;
      for i = 1:numel(keysCa)
        key = keysCa(i);

        value1 = obj1.Dict(key).v;
        value2 = obj2.Dict(key).v;

        % NOTE: NaN == NaN ==> Use isequaln().
        if ~isequaln(value1, value2) || ~isequaln(class(value1), class(value2))
          equals = false;
          return
        end
      end

      equals = true;
    end



    function unequal = ne(obj, other)
      unequal = ~obj.eq(other);
    end



    % Overload disp(Srm). Useful for debugging e.g. tests.
    %
    % Unclear if works for non-numeric, non-char key values. Probably not.
    %
    function s = disp(obj)
      keysCa = obj.Dict.keys();

      sCa = {sprintf('%i row(s)', obj.nRows())};
      for i = 1:numel(keysCa)
        key   = keysCa(i);
        value = subsref(obj, substruct('()', {key}));

        sCa{end+1} = sprintf('%s : %s (%s)', num2str(key), mat2str(value), class(value));
      end

      s = strjoin(sCa, '\n');
    end



  end    % methods(Access=public)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function assert_legal_key(key)
      assert(isnumeric(key) || isstring(key))
    end



  end    % methods



end
