%
% Class for storing one arbitrary array together with information on which
% elements are "fill positions" (elements which values should be ignored)
% without resorting to using fill values to represent them. Analogous to, and
% inspired by, JUICE/RPWI GS pipeline's class FillPositionsArray.
%
% Immutable.
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
    %   PROPOSAL: Support assigning index.
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
    %   PROPOSAL: Reorder arguments to have initType first.
    %
    % PROPOSAL: More static constructor wrapper methods
    %   PRO: Useful for automatic tests.
    %       PROPOSAL: doubleNan --> Specified int



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(GetAccess=private, SetAccess=immutable)
        % NOTE: Should be completely private, but in practice it is possible to
        % read "dataAr" under MATLAB R2019b. Unknown why. Property is
        % write-protected though. Update test w.r.t. to this if fixed.
        dataAr
    end


    
    properties(GetAccess=public, SetAccess=immutable)
        % Logical array of same size as dataAr. True<=>The corresponding
        % position in dataAr is a fill position where the value is irrelevant
        % and must be hidden from the user.
        fpAr
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
        %       'fill value'
        %           Fill value. Same type as dataAr. Scalar.
        %           dataAr elements with this value will be identified as fill
        %           positions.
        %       'fill positions'
        %           Fill positions array. Logical. Same size as dataAr.
        %           True == Fill position.
        %       
        function obj = FillPositionsArray(dataAr, fpDescriptionType, fpDescription)

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
                case 'fill value'
                    fillValue = fpDescription;
                    assert(isscalar(fillValue))
                    assert(strcmp(class(fillValue), class(dataAr)), ...
                        'Fill value and data have different MATLAB classes.')

                    % NOTE: Array operation. Can not use isequaln(). Needs
                    % special case for NaN.
                    fpAr = (dataAr == fillValue) | (isnan(dataAr) & isnan(fillValue));
                    clear fillValue

                case 'fill positions'
                    fpAr = fpDescription;
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
        end



        % Return array with fill positions filled in with specified fill value.
        function dataAr = get_data(obj, fillValue)
            assert(isscalar(fillValue))
            assert(...
                strcmp(class(fillValue), class(obj.dataAr)), ...
                'Argument fillValue has a MATLAB class ("%s") which is inconsistent with the object''s MATLAB class ("%s").', ...
                class(fillValue), class(obj.dataAr))

            dataAr           = obj.dataAr;
            dataAr(obj.fpAr) = fillValue;
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
        % fillValueBefore
        %       Fill value used for the input to the function.
        %
        %
        % RETURN VALUE
        % ============
        % Fpa
        %       FPA. Operation output array elements for fill positions will be
        %       ignored in the new FPA.
        %
        function Fpa = convert(obj, fhArrayOperation, outputClass, fillValueBefore)
            % IMPLEMENTATION NOTE: Can not use (a) arbitrary fill values, or (b)
            % just use the values stored in obj.dataAr since the array operation
            % might trigger error for some values, e.g. ~NaN (negating NaN, but
            % NaN can not be interptered as logical/boolean).
            
            % TODO-DEC: Should the function (1) require the output to be on the
            % specified MATLAB class/type, or (2) cast the output to the
            % specified MATLAB class/type?
            assert(isa(fhArrayOperation, 'function_handle'))
            
            inputAr      = obj.get_data(fillValueBefore);
            outputAr     = fhArrayOperation(inputAr);
            castOutputAr = cast(outputAr, outputClass);
            
            Fpa = bicas.utils.FillPositionsArray(...
                castOutputAr, 'fill positions', obj.fpAr);
        end
        
        
        
        % Convert to other MATLAB class.
        %
        % ARGUMENTS
        % =========
        % fillValueBefore
        %       Fill value (before cast) that can survive the cast.
        %
        function Fpa = cast(obj, outputClass, fillValueBefore)
            % PROPOSAL: Make fillValue argument optional for cases where
            %           fillValue could be automatically derived depending on
            %           conversion.
            
            Fpa = obj.convert(...
                @(x) (cast(x, outputClass)), ...
                outputClass, fillValueBefore);
        end
        
        
        
        % Utility function
        function data = int2doubleNan(obj)
            assert(isinteger(obj.dataAr), 'FPA is not integer. It is of class "%s".', class(obj.dataAr))
            
            fillValue = cast(0, class(obj.dataAr));
            Fpa  = obj.cast('double', fillValue);
            data = Fpa.get_data(NaN);
        end
        
        
        
        function data = logical2doubleNan(obj)
            assert(islogical(obj.dataAr))
            
            Fpa  = obj.cast('double', false);
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
                        dataAr, 'fill positions', fpAr)};

                case '.'
                    % Call method (sic!)
                    [varargout{1:nargout}] = builtin('subsref', obj, S);

                otherwise
                    error('BICAS:Assertion', 'Does not support operation.')
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



    end    % methods(Access=public)
    
    
    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static, Access=public)
        
        
        
        % Wrapper around constructor. Effectively custom constructor.
        function Fpa = floatNan2logical(ar)
            assert(isfloat(ar))
            Fpa = bicas.utils.FillPositionsArray(...
                ar, 'fill value', NaN).cast('logical', 0);
        end


            
    end



end
