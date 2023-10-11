%
% Class for storing one arbitrary array together with information on which
% elements are "fill positions" (elements which values should be ignored)
% without resorting to using fill values to represent them. Analogous to, and
% inspired by, JUICE/RPWI GS pipeline's class FillPositionsArray.
%
%
% RATIONALE
% =========
% Class should be useful for representing variables where individual elements
% can also be unknown (in practice, originate from CDF fill values). BICAS was
% originally not written with this in mind, but some has been retrofitted to use
% this class, but far from all. Using this class more should however be the
% natural approach to solve problems associated with fill values/unknown values
% in the future.
%
%
% ~BUG
% ====
% The implementation of subsref() seems to block instances of the (value) class
% from being modified. Unknown why.
%
%
% NAMING CONVENTION
% =================
% FP: Fill Position.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FPArray   % < handle
    % PROPOSAL: Shorter name.
    %   PRO: Frequently referenced.
    %       PRO: 2023-10-09: 98 occurrences of class name.
    %   PROPOSAL: FPArray
    %
    % PROPOSAL: Change name of "get_data()"
    %   PROPOSAL: array, data, data_array
    % PROPOSAL: Better name for method "get_non_FP_data".
    %   NOTE: Cf .get_data().
    %   ~get_non_FP
    %   ~NFP_array    ## cf. .array(fv)
    %   ~elements
    %   ~data
    %   ~clean
    %   ~valid
    %   ~actual
    %   ~real
    %   PROPOSAL: General terms: FP, non-FP data, FP + non-FP data.
    %
    % TODO-DEC: Tolerate cell arrays?
    %   NOTE: Should not be needed for BICAS.
    %
    % TODO-DEC: Naming convention for FPA variables?
    %   *Fpa, Fpa* Fpa_L3_QUALITY_BITMASK
    %   Fpa_L3_QUALITY_BITMASK
    %   ZvFpa_L3_QUALITY_BITMASK
    %   L3_QUALITY_BITMASK_Fpa
    %
    % PROPOSAL: Shorten "doubleNan" --> "dblNan"
    %
    % TODO-DEC: Naming convention for convenience methods for FPA<-->array.
    %   NEED: Short.
    %       PROBLEM: How represent that output is FPA or array, if kept short?
    %
    % PROPOSAL: Performance w.r.t. pre-allocation?
    %   TODO-DEC:  Is it a relevant question? Would such code use FPA?
    %   PROPOSAL: Use HandleWrapper internally.
    %       NOTE: Class needs to become handle class.
    %
    % PROBLEM: How handle specifying a FV in .get_data() that is identical to
    %          some non-FP elements?
    %   PROPOSAL: Assertion in method.
    %       PROPOSAL: Extra argument for whether to have assertion.
    %   PROPOSAL: Make it easy to check before calling method.
    %       PROPOSAL: Method .contains(value) .
    %
    % PROPOSAL: disp() method.
    %   PROPOSAL: Number of FPs and non-FPs of total (percent?).
    %
    % PROBLEM: Can not call e.g. isnumeric(Fpa), isfloat(Fpa), islogical().
    %   PROPOSAL: isnumeric(cast(0, Fpa.mc)) etc.
    %
    % PROBLEM: Can not create empty FPA for arbitrary MATLAB class dynamically (from MATLAB class string).
    %   PROBLEM: How specify in constructor?
    %       PROPOSAL: Syntax(size, 'EMPTY', mc)
    %           CON: Inconsistent compared to other constructor modes.
    %   PROPOSAL: Treat as special case: Do not initialize dataAr until
    %             non-empty.
    %       CON: More work.
    %           CON: Not feasible for implementing subsref, subsasgn, cat.
    %       CON: Ugly.
    %   PROPOSAL: Use eval inside constructor.
    %   PROPOSAL: Restrict supported data types.
    %       PROPOSAL: logical & numeric & char.
    %           PRO: Already done in constructor...
    
    
    
    %##############################
    %##############################
    % "STATIC" CONSTANT PROPERTIES
    %##############################
    %##############################
    properties(Constant)
        % List of all numeric MATLAB classes.
        MC_NUMERIC_CA = {...
            'single', 'double', ...
            'int8', 'uint8', 'int16', 'uint16', ...
            'int32', 'uint32', 'int64', 'uint64' ...
        };
    
        % Constant that is useful for setting FPA elements to be Fs (also for
        % multiple elements using indexing).
        FP_UINT8 = bicas.utils.FPArray.get_scalar_FP('uint8');
    end



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(GetAccess=private, SetAccess=private)
        % NOTE: Should be completely private, but in practice it is possible to
        % read "dataAr" under MATLAB R2019b. Unknown why. Property is
        % write-protected though. Update test w.r.t. to this if fixed.
        dataAr
    end
    
    properties(GetAccess=public, SetAccess=private)
        % Logical array of same size as dataAr. True<=>The corresponding
        % position in dataAr is a fill position where the value is irrelevant
        % and must be hidden from the user.
        fpAr
    end
    
    properties(GetAccess=public, SetAccess=immutable)
        % MATLAB class for the internal data.
        % NOTE: Immutable.
        mc
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        % ARGUMENTS
        % =========
        % dataAr
        %       Array of actual data. Limited to some data types.
        % fpDescriptionType
        %       String constant.
        % fpDescription
        %       Meaning determined by fpDescriptionType:
        %       'FILL_VALUE'
        %           Fill value. Same type as dataAr. Scalar.
        %           dataAr elements with this value will be identified as fill
        %           positions.
        %       'FILL_POSITIONS'
        %           Fill positions array. Logical. Same size as dataAr.
        %           True == Fill position.
        %       
        function obj = FPArray(dataAr, fpDescriptionType, varargin)

            % IMPLEMENTATION NOTE: Limiting the data types to those for which
            % one can use MATLAB's element-wise operations for matrices, e.g.
            % equality.
            % NOTE: Does not permit cell arrays.
            % NOTE: Needs to permit char arrays for metadata ZVs.
            assert(isnumeric(dataAr) || ischar(dataAr) || islogical(dataAr))
            irf.assert.castring(fpDescriptionType)

            % ===========
            % Assign fpAr
            % ===========
            switch(fpDescriptionType)
                case 'FILL_VALUE'
                    assert(numel(varargin) == 1)
                    fv = varargin{1};
                    assert(isscalar(fv))
                    assert(strcmp(class(fv), class(dataAr)), ...
                        'Fill value and data have different MATLAB classes.')

                    % NOTE: Array operation. Can not use isequaln(). Needs
                    % special case for NaN.
                    fpAr = (dataAr == fv) | (isnan(dataAr) & isnan(fv));
                    clear fv

                case 'FILL_POSITIONS'
                    assert(numel(varargin) == 1)
                    fpAr = varargin{1};
                    % NOTE: Assertions on variable come later.

                case 'NO_FILL_POSITIONS'
                    assert(numel(varargin) == 0)
                    fpAr = false(size(dataAr));
                    % NOTE: Assertions on variable come later.

                case 'ONLY_FILL_POSITIONS'
                    assert(numel(varargin) == 0)
                    fpAr = true(size(dataAr));
                    % NOTE: Assertions on variable come later.

                otherwise
                    error('Illegal argument "%s"', fpDescriptionType)
            end
            clear fpDescription
            assert(islogical(fpAr))
            assert(isequal(size(dataAr), size(fpAr)))

            % ====================
            % Assign object fields
            % ====================
            obj.dataAr = dataAr;
            obj.fpAr   = fpAr;
            obj.mc     = class(dataAr);
        end



        % Return array with fill positions filled in with specified fill value.
        %
        % ARGUMENTS
        % =========
        % fv
        %       Fill value to use for fill positions. Is optional for MATLAB
        %       classes for which a value can be automatially derived.
        function dataAr = get_data(obj, fv)
            
            switch(nargin)
                case 1
                    fv = bicas.utils.FPArray.get_cast_FV(obj.mc, obj.mc);
                case 2
                    assert(isscalar(fv))
                    assert(...
                        strcmp(class(fv), obj.mc), ...
                        'Argument fv has MATLAB class ("%s") which is inconsistent with the object''s MATLAB class ("%s").', ...
                        class(fv), obj.mc)
                otherwise
                    error('BICAS:Assertion', 'Illegal number of arguments.')
            end

            dataAr           = obj.dataAr;
            dataAr(obj.fpAr) = fv;
        end
        
        
        
        % Set all positions to FP.
        %
        % NOTE/BUG: Does not modify the object without explicitly using the
        % return value! Call Fpa = Fpa.set_FP() to effectively modify the
        % object.
%         function obj = set_FP(obj)
%             obj.fpAr(:) = true;
%         end
%       Function disabled since not used.
        
        
        
        % Convert FPA to other FPA using a specified ARRAY operation for its
        % data (array-->array).
        %
        % NOTE: Can not be used for combining data from multiple FPAs.
        %
        %
        % ARGUMENTS
        % =========
        % fhArrayOperation
        %       Function handle: outputArray = f(inputArray).
        % outputClass
        %       MATLAB class to which the operation output will be cast. Is the
        %       type of the new FPA.
        % fvBefore
        %       Fill value used for the input to the function.
        %
        %
        % RETURN VALUE
        % ============
        % Fpa
        %       FPA. Operation output array elements for fill positions will be
        %       ignored in the new FPA.
        %
        function Fpa = convert(obj, fhArrayOperation, fvBefore)
            % IMPLEMENTATION NOTE: Can not use
            %     (a) arbitrary fill values,
            %     (b) just use the values stored in obj.dataAr, or
            %     (c) bicas.utils.FPArray.MC_NUMERIC_CA,
            % since the array operation might trigger error for some values,
            % e.g. ~NaN (negating NaN, but NaN can not be interptered as
            % logical/boolean).
            
            assert(isa(fhArrayOperation, 'function_handle'))
            
            inputAr  = obj.get_data(fvBefore);
            outputAr = fhArrayOperation(inputAr);
            
            Fpa = bicas.utils.FPArray(...
                outputAr, 'FILL_POSITIONS', obj.fpAr);
        end
        
        
        
        % Convert FPA to FPA with other MATLAB class.
        %
        % ARGUMENTS
        % =========
        % fvBefore
        %       Optional, if value can be automatically derived.
        %       FV (before cast) that can survive the cast.
        %
        function Fpa = cast(obj, outputMc, fvBefore)
            
            switch(nargin)
                case 2
                    fvBefore = bicas.utils.FPArray.get_cast_FV(...
                        obj.mc, outputMc);
                case 3
                    % Do nothing. Use caller-supplied value of "fbBefore".
                otherwise
                    error('BICAS:AssertionError', 'Illegal number of arguments')
            end
            
            Fpa = obj.convert(...
                @(x) (cast(x, outputMc)), fvBefore);
        end
        
        
        
        % Set fill positions using values from another FPA, unless those
        % elements are also fill positions.
        %
        % NOTE: Creates new instance. ==> Not suitable for pre-allocation.
        function Fpa2 = complement(obj, Fpa1)
            % PROPOSAL: Better name.
            %   ~set fill positions
            %   ~fill in (holes)
            %   ~fill in (fill positions)
            %   ~complement
            %   PROPOSAL: Reverse roles of FPAs. Name it ~override, ~overlay.

            assert(strcmp(obj.mc, class(Fpa1.dataAr)))
            
            dataAr           = obj.dataAr;
            dataAr(obj.fpAr) = Fpa1.dataAr(obj.fpAr);
            fpAr             = obj.fpAr & Fpa1.fpAr;
            
            Fpa2 = bicas.utils.FPArray(dataAr, 'FILL_POSITIONS', fpAr);
        end
        
        
        
        % Return all elements which are not fill positions.
        %
        % RETURN VALUE
        % ============
        % ar
        %       1D column array. Internal array is first converted to column
        %       vector. Then fill positions are removed.
        function ar = get_non_FP_data(obj)
            % IMPLEMENTATION NOTE: Must convert to column array. Otherwise a 1D
            % vector that is not a column remains a column in that dimension.
            dataAr = obj.dataAr(:);   % Convert to column array.
            ar = dataAr(~obj.fpAr(:));
        end
        
        
        
        % Utility function
        function data = int2doubleNan(obj)
            assert(isinteger(obj.dataAr), 'FPA is not integer. It is of MATLAB class "%s".', obj.mc)
            
            Fpa  = obj.cast('double');
            data = Fpa.get_data(NaN);
        end
        
        
        
        % Utility function
        function data = logical2doubleNan(obj)
            assert(islogical(obj.dataAr))
            
            Fpa  = obj.cast('double');
            data = Fpa.get_data(NaN);
        end



        % Operator overloading: ==
        %
        % NOTE: Ignores values at fill positions.
        % NOTE: NaN counts as equal to itself.
        %
        %
        % RETURN VALUE
        % ============
        % r
        %       Scalar. Logical. Whether the two FPAs are equal.
        %       NOTE: Always scalar. Therefore can not do element-wise
        %       comparison.
        function r = eq(obj1, obj2)
            % Tentatively permit subclassing of this class, though requiring
            % identical classes. Subclass should be able to use this method to
            % e.g. implement its own counterpart method. Unclear what is the
            % best behaviour.
            if ~strcmp(class(obj1), class(obj2))
                r = false;
                
            elseif ~strcmp(class(obj1.dataAr), class(obj2.dataAr))
                r = false;
                
            elseif ~isequaln(obj1.fpAr, obj2.fpAr)
                % NOTE: Indirectly checks equal .dataAr sizes.
                r = false;
                
            else
                % CASE: .fpAr, sizes and types are equal.
                
                if isempty(obj1.dataAr)
                    r = true;
                else
                    % CASE: Arrays are not empty.
                    fv = obj1.dataAr(1);   % Arbitrary fill value.
                    % NOTE: NaN == NaN
                    r = isequaln(obj1.get_data(fv), obj2.get_data(fv));
                end
            end
        end



        % Operator overloading: ~=
        function r = ne(obj1, obj2)
            r = ~obj1.eq(obj2);
        end
        
        
        
        function r = isequaln(obj1, obj2)
            r = obj1.eq(obj2);
        end

        function r = isequal(obj1, obj2)
            error('BICAS:OperationNotImplemented', 'isequal() has not been implemented.')
        end



        % Indexing overloading: Array indexing for reading: Fpa(i, j, ...)
        function varargout = subsref(obj, S)
            switch S(1).type
                case '()'
                    dataAr = subsref(obj.dataAr, S);
                    fpAr   = subsref(obj.fpAr, S);

                    varargout = {bicas.utils.FPArray(...
                        dataAr, 'FILL_POSITIONS', fpAr)};

                case '.'
                    % Call method (sic!)
                    
                    % NOTE/BUG: This seems to prevent methods from modifying the
                    % (value) object for unknown reason. Can for example not
                    % implement a method that mutates obj.
                    %   function obj = set_FP(obj)
                    %       obj.fpAr(:) = true;
                    %   end
                    % when calling
                    %   FpFpa.set_FP();
                    % For that call, varargout == {obj} though.
                    % This works though:
                    %   FpFpa = FpFpa.set_FP();
                    [varargout{1:nargout}] = builtin('subsref', obj, S);

                otherwise
                    error('BICAS:Assertion', 'Does not support operation.')
            end
        end
        
        
        
        % Indexing overloading: Array indexing for writing: Fpa(i, j, ...) = ...
        %
        % Supports assigning from (1) either a (1) FPA, or (2) array, either of
        % which must have (a) the same MATLAB class, and (b) compatible size
        % (either identical or scalar).
        %
        %
        % PERFORMANCE
        % ===========
        % Testing (bicas.utils.FPArray___subsasgn_SpeedTest) implies
        % that preallocating a large FPA and then overwriting subsets using
        % subsasgn does not work. Time consumption per element grows with size
        % of FPA. Class should thus not be suitable for storing samples in the
        % processing step which updates pre-allocated global array of samples.
        %
        function Fpa1 = subsasgn(Fpa1, S, obj2)
            % TODO-DEC: Is it appropriate that non-FPAs can be used to assign
            %           FPAs? Could lead to bugs in transition from non-FPAs to
            %           FPAs.
            
            switch S(1).type
                case '()'
                    assert(isscalar(S))
                    % IMPLEMENTATION NOTE: Check that index is not some
                    % array-like objet, e.g. FPA itself. Could maybe support FPA
                    % in the future(?!!).
                    for i = 1:numel(S.subs)
                        x = S.subs{i};
                        assert(isnumeric(x) || islogical(x) || strcmp(x, ':'))
                    end
                    
                    if isa(obj2, 'bicas.utils.FPArray')
                        assert(strcmp(Fpa1.mc, obj2.mc), ...
                            'FPA to be assigned has a MATLAB class "%s" that incompatible with the assigned FPA''s MATLAB class "%s".', ...
                            obj2.mc, Fpa1.mc)
                        dataAr2 = obj2.dataAr;
                        fpAr2   = obj2.fpAr;
                    % DISABLE assignment using non-FPAs (e.g. regular arrays).
                    % else
                    %     assert(strcmp(Fpa1.mc, class(obj2)), ...
                    %         'Value to be assigned has a MATLAB class "%s" that is incompatible with the assigned FPA''s MATLAB class "%s".', ...
                    %         class(obj2), Fpa1.mc)
                    %     dataAr2 = obj2;
                    %     fpAr2   = false(size(obj2));
                    else
                        error('BICAS:AssertionError', 'Assigning FPA elements with non-FPA.')
                    end
                    
                    Fpa1.dataAr = subsasgn(Fpa1.dataAr, S, dataAr2);
                    Fpa1.fpAr   = subsasgn(Fpa1.fpAr,   S, fpAr2);

                otherwise
                    error('BICAS:Assertion', 'Unsupported operation.')
            end
        end



        % "Overload" size(Fpa, ...)
        function s = size(obj, varargin)
            s = size(obj.dataAr, varargin{:});
        end



        % "Overload" ndims(Fpa, ...)
        function n = ndims(obj, varargin)
            n = ndims(obj.dataAr, varargin{:});
        end



        % "Overload" isempty(Fpa)
        function n = isempty(obj)
            n = isempty(obj.dataAr);
        end



        % Overload "end" in indexing.
        function iEnd = end(obj, iDim, nDim)
            sz = size(obj.dataAr);

            if iDim < nDim
                iEnd = sz(iDim);
            else
                iEnd = prod(sz(iDim:end));
            end
        end



        % Overload <
        function r = lt(obj1, obj2)
            r = obj1.binary_operation_to_FPA(obj2, @(o1, o2) (o1 < o2));
        end

        % Overload >
        function r = gt(obj1, obj2)
            r = obj1.binary_operation_to_FPA(obj2, @(o1, o2) (o1 > o2));
        end

        % Overload <=
        function r = le(obj1, obj2)
            r = obj1.binary_operation_to_FPA(obj2, @(o1, o2) (o1 <= o2));
        end

        % Overload >=
        function r = ge(obj1, obj2)
            r = obj1.binary_operation_to_FPA(obj2, @(o1, o2) (o1 >= o2));
        end


        % "Overload" cat().
        % NOTE: First argument should NOT be the instance of the class!!!
        function Fpa = cat(iDim, varargin)
            % ASSERTIONS
            assert(isnumeric(iDim))
            FpaCa = varargin;
            for i = 1:numel(FpaCa)
                assert(isa(FpaCa{i}, 'bicas.utils.FPArray'))
            end
            mcCa = cellfun(@(Fpa) Fpa.mc, FpaCa, 'UniformOutput', false);
            assert(numel(unique(mcCa)) == 1)
            
            % "ALGORITHM"
            dataArCa = cellfun(@(Fpa) Fpa.dataAr, FpaCa, 'UniformOutput', false);
            dataAr   = cat(iDim, dataArCa{:});
            fpArCa   = cellfun(@(Fpa) Fpa.fpAr, FpaCa, 'UniformOutput', false);
            fpAr     = cat(iDim, fpArCa{:});

            Fpa = bicas.utils.FPArray(dataAr, 'FILL_POSITIONS', fpAr);
        end

        % "Overload" vertcat() = [... ; ...].
        function Fpa = vertcat(varargin)
            Fpa = cat(1, varargin{:});
        end

        % "Overload" horzcat() = [..., ...]
        function Fpa = horzcat(varargin)
            Fpa = cat(2, varargin{:});
        end



    end    % methods(Access=public)



    %##########################
    %##########################
    % PRIVATE INSTANCE METHODS
    %##########################
    %##########################
    methods(Access=private)
        


        % Helper function to make it easier to implement operator overloading
        % for binary operators.
        %
        % NOTE: Always outputs an FPA.
        function Fpa3 = binary_operation_to_FPA(obj1, obj2, fhCompare)
            % Concerning potential names for method:
            % "More specifically, an internal binary operation on a set is a
            % binary operation whose two domains and the codomain are the same
            % set."
            % /https://en.wikipedia.org/wiki/Binary_operation
            % =================================================================
            % PROPOSAL: Only compare data from .get_data()?
            %   PRO: Safer.
            %   CON: Need to determine safe fill value to use.
            %       PRO: Can not do for any data type.
            %   PROPOSAL: Add argument for FVs.
            % PROPOSAL: Make public.
            %   PRO: Compare .convert().
            
            
            if isa(obj2, 'bicas.utils.FPArray')
                dataAr = fhCompare(obj1.dataAr, obj2.dataAr);
                fpAr   = (obj1.fpAr | obj2.fpAr);
            else
                dataAr = fhCompare(obj1.dataAr, obj2);
                fpAr   = obj1.fpAr;
            end
            Fpa3 = bicas.utils.FPArray(dataAr, 'FILL_POSITIONS', fpAr);
        end
        
        
        
    end
    
    
    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static, Access=public)



        % Return scalar (1x1) FPA containing one FP for a specified MATLAB
        % class.
        function Fpa = get_scalar_FP(mc)
            fv = bicas.utils.FPArray.get_cast_FV(mc, mc);
            Fpa = bicas.utils.FPArray(fv, 'ONLY_FILL_POSITIONS');
        end
        
        
        
        % Wrapper around constructor. Effectively custom constructor.
        %
        % NOTE: Requires input values to be [0, 1, NaN].
        function Fpa = floatNan2logical(ar)
            assert(isfloat(ar))            
            assert(all(ismember(ar, [0,1]) | isnan(ar)))   % All elements are [0,1,NaN].
            
            floatNaN  = cast(NaN, class(ar));
            
            Fpa = bicas.utils.FPArray(...
                ar, 'FILL_VALUE', floatNaN).cast('logical');
        end
        
        
        
        function Fpa = floatNan2int(ar, fpaMc)
            assert(isfloat(ar))
            assert(isinteger(cast(0, fpaMc)))

            floatNan  = cast(NaN, class(ar));
            
            Fpa = bicas.utils.FPArray(...
                ar, 'FILL_VALUE', floatNan).cast(fpaMc);
        end
        
        
        
    end
    
    
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Return FV which will survive a cast from MATLAB class mc1 to mc2, for
        % some pre-defined MATLAB classes. 
        %
        % NOTE: FPAs as such can handle more classes than specified here.
        %
        % RATIONALE: Having this function saves time not having to figure out
        %            which MATLAB class to specify, which is often error prone.
        %
        function fv = get_cast_FV(mc1, mc2)
            isNum1 = ismember(mc1, bicas.utils.FPArray.MC_NUMERIC_CA);
            isNum2 = ismember(mc2, bicas.utils.FPArray.MC_NUMERIC_CA);
            isLog1 = strcmp(mc1, 'logical');
            isLog2 = strcmp(mc2, 'logical');

            if isNum1 && (isNum2 || isLog2)
                % NOTE: Can not use NaN for numeric-->logical.
                fv = cast(0, mc1);
            elseif isLog1 && (isNum2 || isLog2)
                fv = false;
            else
                % CASE: Can not determine MATLAB class and hence FV.
                %fv = [];
                error('BICAS:AssertionError:IllegalArgument', ...
                    'Can not derive value that survives typecasting from "%s" to "%s".', mc1, mc2)
            end
        end
        
        
        
    end



end
