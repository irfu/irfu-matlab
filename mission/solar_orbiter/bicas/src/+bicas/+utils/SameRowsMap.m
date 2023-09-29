%
% Basic map key-->value, where all values must have the same number of rows
% (size in first dimension). Values must be immutable.
%
% Intended for zVariables. May be upgraded to require FPAs.
% Some method names are chosen to be identical with containers.Map.
%
%
% IMPLEMENTATION NOTE
% ===================
% bicas.utils.SameRowsMap.setRows() can be slow if storing data directly as
% values in containers.Map, presumably since preallocation does not work.
% Therefore storing all values indirectly via handle class objects
% (bicas.utils.HandleWrapper), which (apparently) makes it possible to modify
% arrays without implicit copying by MATLAB, thus increasing performance. Does
% seem to work. Since the class's internal data structure uses handle objects,
% the class itself also has to be a handle class, to avoid that internal handle
% objects are shared between different instances of bicas.utils.SameRowsMap.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SameRowsMap < handle
    % PROPOSAL: Automatic test code.
    % PROPOSAL: Better name.
    %   Collection, ZV, Map, Set
    % PROPOSAL: Use syntactic sugar somehow.
    %   set, get, length, nRows, keys
    %   PROPOSAL: Overload indexing: set/add, get
    %   PROPOSAL: Properties for nRows, keys
    %
    % NOTE: Examples of usage for ZV collections:
    %   currentTmCa
    %   AsrMap (while building)
    %   AsrMap (final) / DemuxerOutput
    %   PreDc/PostDc.Zv
    %
    % PROPOSAL: Fixed set of keys.
    % PROPOSAL: Immutable.
    %   PRO: Makes inheritance SameRowsMap-->SameSizeTypeMap natural.
    %   CON: Less convenient for building AsrMap which can be variable-number of
    %        variables.
    %       PROPOSAL: Build content using containers.Map. Submit to constructor.
    %           CON: Can never extend to accept wider set of keys, e.g. objects.
    %           PRO: Easier to implement.
    %
    % PROPOSAL: Constrain MATLAB classes of values: numeric, logical(?)
    %   CON: Unnatural constraint. Assertion on same rows will reveal whether
    %        values are ~array/matrix-like or not.
    %
    % TODO-DEC: How implement relationships between similar classes of ZV
    %           collections?
    %   NEED: Shared test code.
    %   PROPOSAL: Merge SameRowsMap and SameSizeTypeMap. Constructor
    %             argument determines which map value properties should be consistent.
    %       PROPOSAL: Submit function handles that convert variable to value
    %                 that should be consistent between variables.
    %       CON: Can not share all methods, since some methods only makes sense
    %            for one of the classes.
    %           Ex: SameSizeTypeMap:  size():  All dimensions.
    %           Ex: SameSizeTypeMap: ~class(): MATLAB class.
    %               CON: Should not be needed.
    %       NOTE: Might still want to set consistent values (e.g. nRows) in
    %             constructor to handle case of zero keys.
    %   PROPOSAL: SameSizeTypeMap is subclass of SameRowsMap.
    %       PRO: Can add methods size(), class() to subclass.
    %       CON: Subclass and superclass are identical w.r.t. reading, but not
    %            w.r.t. writing content. Subclass adds more assertions.
    %           PRO: assert(isa(x, 'SameRowsMap')) on an argument for reading
    %                makes sense, but not on an argument for writing to.
    %       PROPOSAL: Make immutable.
    %           CON: Can not implement setRows()
    %                as method (without creating new object).
    %   PROPOSAL: SameSizeTypeMap, SameRowsMap are subclasses of
    %             same abstract superclass.
    %       NOTE: Must be abstract superclass to avoid same superclass-subclass
    %             problems.
    %   PROPOSAL: SameSizeTypeMap composes/wraps SameRowsMap.
    %       CON: Must redefine all methods, just for wrapping.
    %   PROBLEM: How compare SameRowsMap with SameSizeTypeMap?
    %
    % NEED: Simultaneously store collections of same-rows ZVs and same-size&type
    %       ZVs and enforce assertions between them.
    %   PROPOSAL: Class(es) which support a hierarchy of collections.
    %       SameSizeType inside/under SameRows.
    %   PROPOSAL: Collections which are connected to each other (keep
    %       references).
    %
    % TODO-DEC: How handle properties (method return values) which are undefined
    %           before setting any variables?
    %   Ex: nRows, size, class
    %   PROPOSAL: Set them in constructor.
    %       PRO: Makes values well-defined.
    %   PROPOSAL: Return special value if unknown: [], NaN etc.
    %   PROPOSAL: Assertion against calling function.
    %
    % PROBLEM: How compare sets of keys, when keys can be both strings and
    %          numbers?
    %   PROPOSAL: Convert numbers to strings. Then compare sets of strings.
    %       CON: Ugly.
    %       CON: Different MATLAB class representations of the same number
    %            convert to the same string.
    


    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(Access=private)
        % containers.Map
        Map

        % IMPLEMENTATION NOTE: Can not use name "nRows" since it coincides with
        % method name.
        nRows2
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = SameRowsMap(keyType, nRows, initType, varargin)
            assert(isnumeric(nRows) && nRows >= 0)
            
            obj.nRows2 = nRows;
            obj.Map = containers.Map('KeyType', keyType, 'ValueType', 'any');
            
            switch(initType)
                case 'empty'
                    assert(numel(varargin) == 0)
                    
                case 'constant'
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


        
        function n = length(obj)
            n = obj.Map.length;
        end


        
        % NOTE: Method name chosen to be identical with containers.Map.keys().
        function keysCa = keys(obj)
            keysCa = obj.Map.keys();
            keysCa = keysCa(:);
        end
        
        
        
        % Mostly for debugging.
        function valuesCa = values(obj)
            hwCa     = obj.Map.values();
            valuesCa = cellfun(@(x) (x.v), hwCa, 'UniformOutput', false);
            valuesCa = valuesCa(:);
        end


        
        % NOTE: Method name chosen to be identical with containers.Map.isKey().
        function isKey = isKey(obj, key)
            bicas.utils.SameRowsMap.assert_legal_key(key)

            isKey = obj.Map.isKey(key);
        end


        
        % Add NEW key-value pair. Disallow overwriting.
        function add(obj, key, value)
            bicas.utils.SameRowsMap.assert_legal_key(key)
            assert(~obj.Map.isKey(key))
            assert(obj.nRows2 == size(value, 1))
            
            obj.Map(key) = bicas.utils.HandleWrapper(value);
        end
        
        
        
        % Use another map to overwrite selected rows in this map.
        %
        % IMPLEMENTATION NOTE: Method is important for speeding up LFR-SWF which
        % tends to be broken into subsequences of 1 record. This can be done by
        % pre-allocating and then overwriting parts.
        %
        %
        % ARGUMENTS
        % =========
        % Map2
        %       bicas.utils.SameRowsMap. Must have the same set of keys and
        %       value types.
        % iRowsArray
        %       Column array. Same length as number of rows in Map2 fields.
        %       Specifies the rows that shall be overwritten.
        %       NOTE: Can not use logical indexing.
        %
        %
        function setRows(obj, Map2, iRowsArray)
            % NOTE: Could add support for other (future) custom-made types of
            %       Maps.
            % PROPOSAL: Support logical indexing.
            % PROPOSAL: Support using an 1-row SRM for overwriting N rows.
            % PROPOSAL: Implement method by overloading indexing notation:
            %           subsasgn.
            %   PRO: Using subsasgn() for assigning internal arrays might handle
            %        (1) logical indexing), (2) assigning multiple rows with
            %        1-row SRM.
            assert(isa(Map2, 'bicas.utils.SameRowsMap'))
            assert(isnumeric(iRowsArray) && iscolumn(iRowsArray))
            assert(size(iRowsArray, 1) == Map2.nRows)
            
            assert(bicas.utils.SameRowsMap.key_sets_equal(obj.keys, Map2.keys))
            
            keysCa = obj.keys();
            for keyCa = keysCa(:)'
                key = keyCa{1};
                
                hw1 = obj.Map(key);
                hw2 = Map2.Map(key);
                
                size1 = size(hw1.v);
                size2 = size(hw2.v);
                % IMPLEMENTATION NOTE: Using num2str(key) since it can handle
                % both strings and numbers.
                assert(isequal(size1(2:end),     size2(2:end)    ), ...
                    'Values for key="%s" have inconsistent sizes.', num2str(key))
                assert(isequal(class(hw1.v), class(hw2.v)), ...
                    'Values for key="%s" have inconsistent MATLAB classes.', num2str(key))

                % IMPLEMENTATION NOTE: Unsure, but think that explicitly setting
                % second dimension to ":" makes the command handle any
                % dimensionalities (assuming that dimensionalities and sizes are
                % consistent).
                hw1.v(iRowsArray, :) = hw2.v(:, :);
                
                % IMPLEMENTATION NOTE: Does not need to set obj.Map(key) since
                % using handle classes.
            end
        end



        function value = get(obj, key)
            value = obj.Map(key).v;
        end



        function nRows = nRows(obj)
            nRows = obj.nRows2;
        end
        
        
        
        function equals = eq(obj1, obj2)
            assert(isa(obj2, 'bicas.utils.SameRowsMap'))
            
            % IMPLEMENTATION NOTE: Must support the case of zero keys.
            
            if ~bicas.utils.SameRowsMap.key_sets_equal(obj1.keys, obj2.keys)
                equals = false;
                return
            elseif obj1.nRows ~= obj2.nRows
                equals = false;
                return
            elseif ~isequal(obj1.Map.KeyType, obj2.Map.KeyType)
                equals = false;
                return
            end
            
            keysCa = obj1.Map.keys;
            for i = 1:numel(keysCa)
                key = keysCa{i};
                
                % IMPLEMENTATION NOTE: Using methods for accessing data bypasses
                % the indirection introduced by using bicas.utils.HandleWrapper.
                value1 = obj1.get(key);
                value2 = obj2.get(key);
                
                % NOTE: NaN == NaN ==> Use isequaln().
                if ~isequaln(value1, value2) || ~isequal(class(value1), class(value2))
                    equals = false;
                    return
                end
            end
            
            equals = true;
        end
        
        
        
        function unequal = ne(obj, other)
            unequal = ~obj.eq(other);
        end



    end    % methods(Access=public)



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static, Access=public)



        function assert_legal_key(key)
            assert(isnumeric(key) || ischar(key))
        end
        
        
        
        function assert_legal_key_set(keysCa)
            keysCa = cellfun(@bicas.utils.SameRowsMap.normalize_key, keysCa, 'UniformOutput', false);
            irf.assert.castring_set(keysCa)
        end



        function equal = key_sets_equal(keySetCa1, keySetCa2)
            % PROPOSAL: Abandon normalization. Use for loops and isequaln().
            %   PROPOSAL: Implement via general-purpose function.
            %   PROPOSAL: Return pair of index arrays to matching elements in each set.
            %       TODO-DEC: How handle input arrays which are not sets
            %       (contain doubles)?
            keySetCa1 = cellfun(@bicas.utils.SameRowsMap.normalize_key, keySetCa1, 'UniformOutput', false);
            keySetCa2 = cellfun(@bicas.utils.SameRowsMap.normalize_key, keySetCa2, 'UniformOutput', false);
            
            equal = isempty(setxor(keySetCa1, keySetCa2));
        end



    end



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % IMPLEMENTATION NOTE: Implemented on the assumption that keys which
        % are strings/char and numeric can be mixed, which is not true in
        % the current implementation since containers.Map does not support
        % it.
        %
        % IMPLEMENTATION NOTE: MATLAB's set operations (setxor) seem to
        % support numbers and strings, but not mixed. Therefore using hack
        % to convert all to strings before set operations on keys.
        %
        % NOTE: Will not distinguish between different numeric types (e.g.
        % double and single).
        function key = normalize_key(key)
            % IMPLEMENTATION NOTE: Add prefix to string in both cases, in
            % case the original string coincides with string version of a
            % number (e.g. string key "3").
            if isnumeric(key)
                key = ['n', num2str(key)];
            else
                key = ['s', key];
            end
        end



    end



end
