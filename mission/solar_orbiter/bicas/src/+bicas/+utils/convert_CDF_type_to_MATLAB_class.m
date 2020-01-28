function matlabClass = convert_CDF_type_to_MATLAB_class(cdfDataType, policy)
%
% Convert a string representing a data type in CDF (spdfcdfread, spdfcdfwrite) to a string representing a MATLAB class
% (type). Useful for correctly type casting or type checking data before writing to a CDF.
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% cdfDataType : String representing a data type. See "policy".
% policy      : String: 'Only CDF data types'   : cdfDataType is interpreted as a CDF data type
%                       'Permit MATLAB classes' : cdfDataType is interpreted as a CDF data type OR a MATLAB class.
% matlabClass : MATLAB class corresponding to cdfDataType.
%
%
% IMPLEMENTATION NOTE
% ===================
% The table for translating/converting data types is based upon table in "spdfcdfwrite.m" with the following
% modifications:
%    1) "tt2000" is NOT a MATLAB type but corresponds to MATLAB's int64.
%    2) "epoch", "epoch16" are not MATLAB types and are not included (yet).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-xx
%


%=================================================================================================================
% Translation table from which variables are derived
% --------------------------------------------------
% Left column     = Legal MATLAB types.
% Right column(s) = Legal types for scpdfcdfwrite which all correspond to the ONE MATLAB type in the left column.
%=================================================================================================================
DATA = {...
    'int8',    {'CDF_INT1', 'CDF_BYTE'};
    'int16',   {'CDF_INT2'}; ...
    'int32',   {'CDF_INT4'}; ...
    'int64',   {'CDF_INT8', 'CDF_TIME_TT2000', 'tt2000'}; ...
    'uint8',   {'CDF_UINT1'}; ...
    'uint16',  {'CDF_UINT2'}; ...
    'uint32',  {'CDF_UINT4'}; ...
    'single',  {'CDF_FLOAT',  'CDF_REAL4'}; ...
    'double',  {'CDF_DOUBLE', 'CDF_REAL8'}; ...
    'char',    {'CDF_CHAR',   'CDF_UCHAR'}; ...
    };
MATLAB_CLASSES = DATA(:, 1);
CDF_TYPES      = DATA(:, 2);



switch(policy)
    case 'Only CDF data types';   permitMatlabClasses = 0;
    case 'Permit MATLAB classes'; permitMatlabClasses = 1;
    otherwise error('BICAS:convert_CDF_type_to_MATLAB_class:Assertion:IllegalArgument', 'Illegal policy.');
end
        
for i = 1:size(DATA, 1)
    if any(strcmp(cdfDataType, CDF_TYPES{i}))
        matlabClass = MATLAB_CLASSES{i};
        return
    end
    if permitMatlabClasses && strcmp(cdfDataType, MATLAB_CLASSES{i})
        matlabClass = MATLAB_CLASSES{i};
        return
    end
end
error('BICAS:convert_CDF_type_to_MATLAB_class:Assertion:IllegalArgument', 'Does not recognize CDF variable type "%s".', cdfDataType)

end
