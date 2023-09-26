%
% Basic map key-->value, where all values must have the same number of rows
% (size in first dimension). Values must be immutable.
%
% Intended for zVariables. May be upgraded to require FPAs.
% Some method names are chosen to be identical with containers.Map.
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
                        'Illegal argument initType=%s', initType)
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
            valuesCa = obj.Map.values();
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
            
            obj.Map(key) = value;
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
        %
        %
        function setRows(obj, Map2, iRowsArray)
            % NOTE: Could add support for other (future) custom-made types of
            % Maps.
            assert(isa(Map2, 'bicas.utils.SameRowsMap'))
            assert(isnumeric(iRowsArray) && iscolumn(iRowsArray))
            assert(bicas.utils.SameRowsMap.key_sets_equal(obj.keys, Map2.keys))
            
            keysCa = obj.keys();
            for keyCa = keysCa(:)'
                key = keyCa{1};
                
                value1 = obj.get(key);
                value2 = Map2.get(key);
                
                size1 = size(value1);
                size2 = size(value2);
                assert(isequal(size1(2:end ), size2(2:end )))
                assert(isequal(class(value1), class(value2)))

                % IMPLEMENTATION NOTE: Unsure, but think that explicitly setting
                % second dimension to ":" makes the command handle any
                % dimensionalities (assuming that dimensionalities and sizes are
                % consistent).
                value1(iRowsArray, :) = value2(:, :);
                
                obj.Map(key) = value1;
            end
        end



        function value = get(obj, key)
            value = obj.Map(key);
        end



        function nRows = nRows(obj)
            nRows = obj.nRows2;
        end
        
        
        
        function equals = eq(obj, other)
            assert(isa(other, 'bicas.utils.SameRowsMap'))
            
            % IMPLEMENTATION NOTE: Must support the case of zero keys.
            
            if ~bicas.utils.SameRowsMap.key_sets_equal(obj.keys, other.keys)
                equals = false;
                return
            elseif obj.nRows ~= other.nRows
                equals = false;
                return
            elseif ~isequal(obj.Map.KeyType, other.Map.KeyType)
                equals = false;
                return
            end
            
            keysCa = obj.Map.keys;
            for i = 1:numel(keysCa)
                key = keysCa{i};
                
                value1 = obj.get(key);
                value2 = other.get(key);
                
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
