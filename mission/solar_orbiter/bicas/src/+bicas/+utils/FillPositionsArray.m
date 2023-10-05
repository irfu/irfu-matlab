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
% NAMING CONVENTION
% =================
% FP: Fill Position.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FillPositionsArray   % < handle
    % TODO-DEC: Tolerate cell arrays?
    %   NOTE: Should not be needed for BICAS.
    %
    % TODO-DEC: Naming convention for FPA variables?
    %   *Fpa, Fpa* Fpa_L3_QUALITY_BITMASK
    %   Fpa_L3_QUALITY_BITMASK
    %   ZvFpa_L3_QUALITY_BITMASK
    %   L3_QUALITY_BITMASK_Fpa
    %
    % PROPOSAL: Better fpDescriptionType naming convention.
    %   PROPOSAL: SCREAMING_SNAKE_CASE. -- IMPLEMENTED
    %       PRO: More conventional (looks)
    %           PRO: Used by JUICE/RPWI GS pipeline.
    %       PRO: Stands out visually.
    %           PRO: Resembles convention for constants.
    %       PRO: Easier to do string replacement.
    %       Ex: FILL_VALUE, FV
    %       Ex: FILL_POSITIONS, FILL_POS, FP
    %       Ex: FP_TRUE, FP_FALSE, NO_FILL_POSITIONS, ONLY_FILL_POSITIONS
    %
    % PROPOSAL: Shorten "doubleNan" --> "dblNan"
    %
    % PROPOSAL: Convenience methods for common conversion FPA<-->array
    %   Ex:
    %       float array --> logical FPA
    %       logical FPA --> float array
    %   CON: Tests?
    %   TODO-DEC: Naming convention
    %       NEED: Short.
    %       PROBLEM: How represent that output is FPA or array, if kept short?
    %
    % PROPOSAL: Method for data class (type).
    %   PRO: Can use for assertions.
    %   TODO-DEC: Name?
    %       .class, .matlabClass
    %   PROPOSAL: Read-only property.
    %
    % PROPOSAL: Performance w.r.t. pre-allocation?
    %   TODO-DEC:  Is it a relevant question? Would such code use FPA?
    %   PROPOSAL: Use HandleWrapper internally.
    %       NOTE: Class needs to become handle class.
    %
    % PROBLEM: How support operation Fpa1(b) = Fpa2(b) ?
    %   Ex: Using option BIAS_HK_LFR_SCI in
    %       bicas.proc.L1L2.lfr.process_CDF_to_PreDC().
    %       Using BDM from LFR when possible, but BIAS HK when not.
    %   PROPOSAL: Support assigning index. -- IMPLEMENTED
    %       CON: Makes class mutable.
    %       PRO: More general.
    %       PRO: Should not be difficult to implement assignment 
    %            Fpa1(<index1>) = Fpa2 since can just directly pass on indexing
    %            to internal arrays.
    %            NOTE: Can already use indices for Fpa2 (subsref()).
    %   PROPOSAL: Method: Fpa3 = Fpa1.set_FPs(Fpa2) -- IMPLEMENTED
    %       NOTE: Can not be trivially reduced to indexing operation (subsasgn).
    %
    % PROPOSAL: Support more constructor modes:
    %   PROPOSAL: All FPs false. Requires only data.
    %   PROPOSAL: All FPs true. Requires only array size.
    %   PROPOSAL: Reorder arguments to have fpDescriptionType first.
    %       CON: Unnecessary. Still makes sense to always have data first.
    %
    % PROPOSAL: More static constructor wrapper methods
    %   PRO: Useful for automatic tests.
    %       PROPOSAL: doubleNan --> Specified int
    
    
    
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
        % MATLAB class for internal data.
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
        function obj = FillPositionsArray(dataAr, fpDescriptionType, varargin)

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

%                 case 'ONLY_FILL_POSITIONS'
%                     assert(numel(varargin) == 0)
%                     fpAr = true(size(dataAr));
%                     % NOTE: Assertions on variable come later.

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
        function dataAr = get_data(obj, fv)
            assert(isscalar(fv))
            assert(...
                strcmp(class(fv), obj.mc), ...
                'Argument fv has a MATLAB class ("%s") which is inconsistent with the object''s MATLAB class ("%s").', ...
                class(fv), obj.mc)

            dataAr           = obj.dataAr;
            dataAr(obj.fpAr) = fv;
        end
        
        
        
        % Convert FPA to other FPA using a specified ARRAY operation
        % (array-->array).
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
        function Fpa = convert(obj, fhArrayOperation, outputClass, fvBefore)
            % IMPLEMENTATION NOTE: Can not use (a) arbitrary fill values, or (b)
            % just use the values stored in obj.dataAr since the array operation
            % might trigger error for some values, e.g. ~NaN (negating NaN, but
            % NaN can not be interptered as logical/boolean).
            
            % TODO-DEC: Should the function (1) require the output to be on the
            % specified MATLAB class/type, or (2) cast the output to the
            % specified MATLAB class/type?
            assert(isa(fhArrayOperation, 'function_handle'))
            
            inputAr      = obj.get_data(fvBefore);
            outputAr     = fhArrayOperation(inputAr);
            castOutputAr = cast(outputAr, outputClass);
            
            Fpa = bicas.utils.FillPositionsArray(...
                castOutputAr, 'FILL_POSITIONS', obj.fpAr);
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
                    fvBefore = bicas.utils.FillPositionsArray.get_cast_FV(...
                        obj.mc, outputMc);
                case 3
                    % Do nothing
                otherwise
                    error('BICAS:AssertionError', 'Illegal number of arguments')
            end
            
            Fpa = obj.convert(...
                @(x) (cast(x, outputMc)), outputMc, fvBefore);
        end
        
        
        
        % Utility function
        function data = int2doubleNan(obj)
            assert(isinteger(obj.dataAr), 'FPA is not integer. It is of class "%s".', obj.mc)
            
            Fpa  = obj.cast('double');
            data = Fpa.get_data(NaN);
        end
        
        
        
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



        % Indexing overloading: Array indexing for reading: Fpa(i, j, ...)
        function varargout = subsref(obj, S)
            switch S(1).type
                case '()'
                    dataAr = subsref(obj.dataAr, S);
                    fpAr   = subsref(obj.fpAr, S);

                    varargout = {bicas.utils.FillPositionsArray(...
                        dataAr, 'FILL_POSITIONS', fpAr)};

                case '.'
                    % Call method (sic!)
                    [varargout{1:nargout}] = builtin('subsref', obj, S);

                otherwise
                    error('BICAS:Assertion', 'Does not support operation.')
            end
        end
        
        
        
        % Indexing overloading: Array indexing for writing: Fpa(i, j, ...) = ...
        %
        % NOTE: Function should currently be unused (except for tests).
        %
        % PERFORMANCE
        % ===========
        % Testing (bicas.utils.FillPositionsArray___subsasgn_SpeedTest) implies
        % that preallocating a large FPA and then overwriting subsets using
        % subsasgn does not work. Time consumption per element grows with size
        % of FPA. Class should thus not be suitable for storing samples in the
        % processing step which updates pre-allocated global array of samples.
        function Fpa1 = subsasgn(Fpa1, S, Fpa2)
            switch S(1).type
                case '()'
                    assert(isscalar(S))
                    assert(isa(Fpa2, 'bicas.utils.FillPositionsArray'))
                    assert(isequaln(Fpa1.mc, Fpa2.mc))
                    
                    % IMPLEMENTATION NOTE: Check that index is not some
                    % array-like objet, e.g. FPA itself. Could maybe support FPA
                    % in the future(?!!).
                    for i = 1:numel(S.subs)
                        x = S.subs{i};
                        assert(isnumeric(x) || islogical(x) || strcmp(x, ':'))
                    end
                    
                    Fpa1.dataAr = subsasgn(Fpa1.dataAr, S, Fpa2.dataAr);
                    Fpa1.fpAr   = subsasgn(Fpa1.fpAr,   S, Fpa2.fpAr);

                otherwise
                    error('BICAS:Assertion', 'Unsupported operation.')
            end
        end



        % Set those elements which are fill positions using values from another
        % FPA.
        %
        % NOTE: Creates new instance. ==> Not suitable for pre-allocation.
        function fpa2 = set_FPs(obj, fpa1)
            % PROPOSAL: Better name.
            % PROPOSAL: Replace by subsasgn(). See BOGIQ.

            assert(strcmp(obj.mc, class(fpa1.dataAr)))
            
            dataAr           = obj.dataAr;
            dataAr(obj.fpAr) = fpa1.dataAr(obj.fpAr);
            fpAr             = obj.fpAr & fpa1.fpAr;
            
            fpa2 = bicas.utils.FillPositionsArray(dataAr, 'FILL_POSITIONS', fpAr);
        end
        
        
        
        % "Overload" size(Fpa, ...)
        function s = size(obj, varargin)
            s = size(obj.dataAr, varargin{:});
        end



        % "Overload" ndims(Fpa, ...)
        function n = ndims(obj, varargin)
            n = ndims(obj.dataAr, varargin{:});
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



    end    % methods(Access=public)
    
    
    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static, Access=public)
        
        
        
        % Wrapper around constructor. Effectively custom constructor.
        %
        % NOTE: Requires input values to be [0, 1, NaN].
        function Fpa = floatNan2logical(ar)
            assert(isfloat(ar))            
            assert(all(ismember(ar, [0,1]) | isnan(ar)))   % All elements are [0,1,NaN].
            
            floatNaN  = cast(NaN, class(ar));
            
            Fpa = bicas.utils.FillPositionsArray(...
                ar, 'FILL_VALUE', floatNaN).cast('logical');
        end
        
        
        
        function Fpa = floatNan2int(ar, fpaMatlabClass)
            assert(isfloat(ar))
            assert(isinteger(cast(0, fpaMatlabClass)))

            floatNan  = cast(NaN, class(ar));
            
            Fpa = bicas.utils.FillPositionsArray(...
                ar, 'FILL_VALUE', floatNan).cast(fpaMatlabClass);
        end
        
        
        
    end
    
    
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Return FV which will survive a cast from MATLAB class mc1 to mc2.
        %
        % NOTE: This function defines the MATLAB classes between which one can
        %       typecast FPAs.
        %
        % RATIONALE: Having this function saved time not having to figure out
        %            which MATLAB class to specify, which is often error prone.
        function fv = get_cast_FV(mc1, mc2)
            isNum1 = ismember(mc1, bicas.utils.FillPositionsArray.MC_NUMERIC_CA);
            isNum2 = ismember(mc2, bicas.utils.FillPositionsArray.MC_NUMERIC_CA);
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
                error('BICAS:AssertionError:IllegalArgument', 'Can not derive value that survives typecasting from "%s" to "%s".', mc1, mc2)
            end
        end


            
    end



end
