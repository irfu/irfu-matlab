% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-xx
%
% Convert a string representing a data type in CDF (spdfcdfread, spdfcdfwrite) to a string representing a MATLAB class
% (type). Useful for correctly type casting or type checking data before writing to a CDF.
%
% Data based upon table in spdfcdfwrite.m with modifications:
%    1) "tt2000" is NOT a MATLAB type but corresponds to int64.
%    2) "epoch", "epoch16" are not MATLAB types and are not included (yet).
%
function MATLAB_class = convert_CDF_type_to_MATLAB_class(cdf_data_type, policy)

%-------------------------------------------------------------------------------------------------------------
% Left column = Legal MATLAB types.
% Right column(s) = Legal types for scpdfcdfwrite which correspond to the ONE MATLAB type in the left column.
%-------------------------------------------------------------------------------------------------------------
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
    case 'Permit MATLAB classes'; permit_MATLAB_classes = 1;
    case 'Only CDF data types';   permit_MATLAB_classes = 0;
    otherwise error('convert_CDF_type_to_MATLAB_class:Assertion:IllegalArgument', 'Illegal policy.');
end
        
for i = 1:size(DATA, 1)
    if any(strcmp(cdf_data_type, CDF_TYPES{i}))
        MATLAB_class = MATLAB_CLASSES{i};
        return
    end
    if permit_MATLAB_classes && strcmp(cdf_data_type, MATLAB_CLASSES{i})
        MATLAB_class = MATLAB_CLASSES{i};
        return
    end
end
error('convert_CDF_type_to_MATLAB_class:Assertion:IllegalArgument', 'Does not recognize CDF variable type "%s".', cdf_data_type)

end
