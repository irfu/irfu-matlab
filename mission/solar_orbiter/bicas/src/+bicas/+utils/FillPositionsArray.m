%
% Class for storing one arbitrary array together with information on which
% positions are "fill positions" (elements which values should be ignored)
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
    % PROPOSAL: Convenience methods for converting
    %       float array --> logical FPA
    %       logical FPA --> float array
    %
    % PROPOSAL: Method for data class (type).
    %   PRO: Can use for assertions.
    %   TODO-DEC: Name?
    %       .class, .matlabClass
    %   PROPOSAL: Read-only property.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(GetAccess=private)
        dataAr

        % Logical array of same size as dataAr. True<=>The corresponding
        % position in dataAr is a fill position where the value is irrelevant
        % (and must be hidden from the user).
        fillPositionsAr
    end


    
    properties
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
            assert(isnumeric(dataAr) || ischar(dataAr) || islogical(dataAr))
            irf.assert.castring(fpDescriptionType)

            % ======================
            % Assign fillPositionsAr
            % ======================
            switch(fpDescriptionType)
                case 'fill value'
                    fillValue = fpDescription;
                    assert(isscalar(fillValue))
                    assert(strcmp(class(fillValue), class(dataAr)))

                    % NOTE: Array operation. Can not use isequaln()
                    fillPositionsAr = (dataAr == fillValue) | (isnan(dataAr) & isnan(fillValue));
                    clear fillValue

                case 'fill positions'
                    fillPositionsAr = fpDescription;
                    % NOTE: Assertions on variable come later.

                otherwise
                    error('Illegal argument "%s"', fpDescriptionType)
            end
            clear fpDescription
            assert(islogical(fillPositionsAr))
            assert(isequal(size(dataAr), size(fillPositionsAr)))

            % ====================
            % Assign object fields
            % ====================
            obj.dataAr          = dataAr;
            obj.fillPositionsAr = fillPositionsAr;
        end



        % Return array with fill positions filled in with specified fill value.
        function dataAr = get_data(obj, fillValue)
            assert(isscalar(fillValue))
            assert(...
                strcmp(class(fillValue), class(obj.dataAr)), ...
                'Argument fillValue has a MATLAB class ("%s") which is inconsistent with the object''s MATLAB class ("%s").', ...
                class(fillValue), class(obj.dataAr))

            dataAr                      = obj.dataAr;
            dataAr(obj.fillPositionsAr) = fillValue;
        end
        
        
        
        % Convert FPA to other FPA using a specified array operation
        % (array-->array)
        %
        % NOTE: Can not be used for combining data from multiple FPAs.
        %
        %
        % ARGUMENTS
        % =========
        % fhArrayOperation
        %       Function handle: outputArray = f(inputArray).
        % fillValue
        %       Fill value used for the input to the function.
        % outputType
        %       MATLAB class to which the operation output will be cast. Is the
        %       type of the new FPA.
        %
        %
        % RETURN VALUE
        % ============
        % FPA
        %       Operation output array elements for fill positions will be
        %       ignored in the new FPA.
        %
        function Fpa = convert(obj, fhArrayOperation, outputClass, fillValue)
            % IMPLEMENTATION NOTE: Can not use (a) arbitrary fill values, or (b)
            % just use the values stored in obj.dataAr since the array operation
            % might trigger error for some values, e.g. ~NaN (negating NaN, but
            % NaN can not be interptered as logical/boolean).
            
            % TODO-DEC: Should the function (1) require the output to be on the
            % specified MATLAB class/type, or (2) cast the output to the
            % specified MATLAB class/type?
            assert(isa(fhArrayOperation, 'function_handle'))
            
            inputAr      = obj.get_data(fillValue);
            outputAr     = fhArrayOperation(inputAr);
            castOutputAr = cast(outputAr, outputClass);
            
            Fpa = bicas.utils.FillPositionsArray(...
                castOutputAr, 'fill positions', obj.fillPositionsAr);
        end
        
        
        
        % Convert to other MATLAB class.
        %
        % ARGUMENTS
        % =========
        % fillValue
        %       Fill value that can survive the cast.
        %
        function Fpa = cast(obj, outputClass, fillValue)
            % PROPOSAL: Make fillValue argument optional for cases where
            %           fillValue could be automatically derived depending on
            %           conversion.
            %   CON: Failure would depend on hidden fill position values.
            
            Fpa = obj.convert(...
                @(x) (cast(x, outputClass)), ...
                outputClass, fillValue);
        end



        % Operator overloading: ==
        function r = eq(obj1, obj2)
            % Tentatively permit subclassing of this class, though requiring
            % identical classes. Subclass should be able to use this method to
            % e.g. implement its own counterpart method.
            if ~strcmp(class(obj1), class(obj2))
                r = false;
                
            elseif ~strcmp(class(obj1.dataAr), class(obj2.dataAr))
                r = false;
                
            elseif ~isequal(obj1.fillPositionsAr, obj2.fillPositionsAr)
                % NOTE: Indirectly checks equal .dataAr sizes.
                r = false;
                
            else
                % CASE: .fillPositionsAr, sizes and types are equal.
                
                if isempty(obj1.dataAr)
                    r = true;
                else
                    % CASE: Arrays are not empty.
                    fv = obj1.dataAr(1);   % Arbitrary fill value.
                    r = isequaln(obj1.get_data(fv), obj2.get_data(fv));
                end
            end
        end



        % Operator overloading: ~=
        function r = ne(obj1, obj2)
            r = ~obj1.eq(obj2);
        end



        % Indexing overloading: Array indexing, Fpa(i, j, ...)
        function varargout = subsref(obj, S)
            switch S(1).type
                case '()'
                    dataAr          = subsref(obj.dataAr, S);
                    fillPositionsAr = subsref(obj.fillPositionsAr, S);

                    varargout = {bicas.utils.FillPositionsArray(...
                        dataAr, 'fill positions', fillPositionsAr)};

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



    end    % methods(Access=public)



    methods



        % Property
        % Make property readonly.
        %
        % IMPLEMENTATION NOTE: If this method is absent, then the statement to
        % assign the property does NOT raise exception, while at the same time
        % doing nothing! Therefore better to have this function. Unsure if there
        % is a better way to implement the same functionality.
        function fpAr = get.fpAr(obj)
            fpAr = obj.fillPositionsAr;
        end



        % Property
        function fpAr = set.fpAr(obj, value)
            error('Illegal assignment of property.')
        end



    end



end
