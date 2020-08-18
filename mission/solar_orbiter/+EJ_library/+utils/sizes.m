%
% Simultaneously
% ** check the size of multiple variables
% ** whether specified dimensions match explicit sizes
% ** whether specified dimensions have identical but arbitrary sizes (e.g. across variables)
% ** return above arbitarily-sized dimensions, so that the caller can e.g. check them.
%
% Ex: Check that variables representing zVariables have explicitly specified non-record dimensions, and the same
% arbitrary number of records (size in first dimension). Return number of records.
%
%
% ARGUMENTS
% =========
% varargin: Arbitrary number of argument pairs below:
%   var            : variable value which size will be tested.
%   sizeConstraint : 1D vector with integers specifying the size of the corresponding argument "var".
%       Nonnegative integer : Explicit required size.
%       Negative integer    : Arbitrary dimension size which must match between all arguments "var".
%                             Must be numbered -1, -2, ... , -N
%       NaN                 : Arbitrary dimension size independent of other dimensions. No corresponding return
%                             value.
%
%
% RETURN VALUES
% =============
% condSatisfied : Logical. Whether variable sizes satisfy criteria.
% varargout     : Size of dimensions labelled with negative integers, in order -1, -2, ... .
%                 NOTE: Values are only guaranteed to be correct if condSatisfied==true.
%
%
% NOTES
% =====
% FINISHED (except assertions on input) but uncertain if
%   has sensible design
%   has sensible name
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-07-29, when broken out from of other code.
%
function [condSatisfied, varargout] = sizes(varargin)
    %   PROPOSAL: Somehow be able to state that a variable is a 1D vector, regardless of which index is not size one.
    %       PROPOSAL: sizeConstraints = {N}, one numeric value (N, negativeValue, NaN).
    %       PROPOSAL: sizeConstraints = {'1D vector', N}
    %       PROPOSAL: Prepend sizeConstraints argument with string constant "vector", "1D", "1D vector".
    %   ~CON/NOTE: Can not assert equal size for variables with arbitrary number of dimensions.
    %
    % PROPOSAL: Abolish NaN (force caller to use negative integer size).
    %   PRO: Simplifies implementation.
    %   PRO: Saves NaN for possible future special value for other functionality.
    %   CON: Clearer that NaN refers to arbitrary value.
    %
    % PROBLEM?: Variables which are so large as to potentially cause memory problems should not(?) be submitted
    % to varargin, since they are likely copied then, even if not modified.
    %   PROPOSAL: Enable way of just submitting size(X) instead.
    %       PROPOSAL: Add suffix argument string constant 'size'.
    %
    % PROPOSAL: Return error message for assertion function to use as error message.
    %   CON: Must add new return value which a direct caller of this function is unlikely to want.
    
    nArgs = numel(varargin);
    
    % Number of return values representing (potentially linked) dimension sizes.
    % NOTE: This should be LESS OR EQUAL to the number of arbitrary (potentially linked) dimension sizes.
    nOutputDims = nargout - 1;
    
    % Assign default value to varargout.
    % IMPLEMENTATION NOTE: Required to avoid secondary error, if exiting function before varargout has been
    % entirely assigned.
    varargout = num2cell(ones(1, nOutputDims) * NaN);
    %varargout = cell(1, nOutputDims);
    
    sizeArray       = [];
    sizeConstrArray = [];
    
    %=====================================================================================
    % Read arguments and store values in a data structure that is suitable for algorithm.
    %=====================================================================================
    iArgPair = 0;    % NOTE: iArgPair=1 represents the first argument pair.
    while true
        %====================
        % Read argument pair
        %====================
        if nArgs == 2*iArgPair
            break
        elseif nArgs >= 2*(iArgPair+1)
            iArgPair = iArgPair + 1;
        else
            error('sizes:Assertion', 'Number of arguments is not even.')
        end
        
        sizeArg           = size(varargin{2*iArgPair-1});
        sizeConstraintArg =      varargin{2*iArgPair  };
        EJ_library.assert.vector(sizeConstraintArg)
        
        %====================================================
        % Add argument pair values into same-sized 1D arrays
        %====================================================
        % Force column arrays.
        sizeArg           = sizeArg(:);
        sizeConstraintArg = sizeConstraintArg(:);
        
        % Pad the smallest arrays with ones (1) until both have same size
        % ---------------------------------------------------------------
        % NOTE: padarray pads in first dimension by default.
        %       'post' : Pads after the last array element along each dimension.
        nDiff = numel(sizeArg) - numel(sizeConstraintArg);
        if nDiff >= 1   sizeConstraintArg = padarray(sizeConstraintArg,  nDiff, 1, 'post');
        else            sizeArg           = padarray(sizeArg,           -nDiff, 1, 'post');
        end
        
        sizeArray       = [sizeArray;       sizeArg];
        sizeConstrArray = [sizeConstrArray; sizeConstraintArg];
    end
    
    
    
    %==========================================================================================================
    % ASSERTION: Arbitrary linked dimension sizes
    % -------------------------------------------
    % IMPLEMENTATION NOTE: Performs this check before explicit dimension sizes to increase chance of correctly
    % assigning varargout before returning.
    %==========================================================================================================
    sizeConstrSpecialValue = -1;
    iOutputDim = 0;
    while true
        b = (sizeConstrArray == sizeConstrSpecialValue);
        if ~any(b)
            % CASE: sizeConstrSpecialValue is not used as a special value.
            % ==> Stop trying other special values.
            break
        end
        
        uniqueSizes  = unique(sizeArray(b));
        nUniqueSizes = numel(uniqueSizes);
        
        if nUniqueSizes ~= 1
            condSatisfied = false;
            return
        end
        %                 assert(nUniqueSizes == 1, ...
        %                     EJ_library.assert.ASSERTION_EMID, ...
        %                     'One or several variables have different sizes for dimensions labelled %i.', sizeConstrSpecialValue)
        
        % Assign return values.
        % NOTE: Does not matter if assigning return values which are not stored by caller, but it is important
        % to that there are at least as many linked dimension sizes as there are return values.
        iOutputDim = -sizeConstrSpecialValue;
        varargout{iOutputDim} = uniqueSizes;
        
        sizeConstrArray(b) = NaN;    % Effectively remove already checked size constraints from later checks.
        
        sizeConstrSpecialValue = sizeConstrSpecialValue - 1;
    end
    assert(nOutputDims <= iOutputDim, ...
        'sizes:Assertion', ...
        'There are more requested output values, than there are arbitrary (potentially linked) dimension sizes.')
    
    
    
    %==============================================
    % ASSERTION: Explicitly stated dimension sizes
    %==============================================
    b = (sizeConstrArray >= 0);    % NOTE: Relies on that (NaN >= 0) returns false.
    if ~all(sizeConstrArray(b) == sizeArray(b))
        condSatisfied = false;
        return
    end
    %             assert(all(sizeConstrArray(b) == sizeArray(b), ...
    %                 EJ_library.assert.ASSERTION_EMID, ...
    %                 'At least one variable dimension size differs from an explicitly specified dimension size.')
    
    sizeConstrArray(b) = NaN;    % Effectively remove already checked size constraints from later checks.
    
    
    
    %===========================================
    % ASSERTION: Only size constraint NaN left.
    %===========================================
    if ~all(isnan(sizeConstrArray))
        error('sizes:Assertion', ...
            ['Size constraints contains negative numbers that can not be interpreted as constraints.'])
    end
    %             assert(all(isnan(sizeConstrArray)), ...
    %                 EJ_library.assert.ASSERTION_EMID, ...
    %                 ['Size constraints contains negative numbers that can not (are not supposed to)', ...
    %                 ' be interpreted as constraints.'])
    
    %==============================================================================
    % NOTE: Ignore size constraint NaN. Does not need to be explicitly checked.
    %==============================================================================
    
    condSatisfied = true;
end
