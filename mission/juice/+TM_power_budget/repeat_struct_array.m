% (Somewhat) generic function for repeating a 1D struct array, and adding fixed offsets to every repetition. Offsets
% increase from zero, by a fixed increment for every copy of the original 1D struct array.
%
%
% Intended for building event sequences for "engine", by repeating an existing sequence.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% Seq0                  : Struct 1D array. Can be row or vector array.
% numFieldName          : String. Name of field in Seq0 which is a scalar (presumably representing time).
% nRepetitionLength     : How much the specified field should be incremented between each repetition.
% nRepetitions          : Positive scalar. Can not be zero.
% Seq                   : Struct 1D array (column).
%
% 
% Created 2018-02-28 by Erik Johansson, IRF Uppsala.
%
function StructArray = repeat_struct_array(StructArray0, numFieldName, nRepetitionOffset, nRepetitions)
% NEED: Consider possible change from struct array to array of structs.
%
% PROPOSAL: Replace by package of utility functions.
%   PRO: Easier to add/remove small functions.
%
% PROPOSAL: Split into functions for
%   (1) repeating struct arrays. and
%   (2) adding to scalar field in struct array.
%   CON: Can wait until actually needed.
%
% PROPOSAL: Assertion for nRepetitions being integer.



% ASSERTIONS
if mod(nRepetitions, 1)~=0    % Test for integer. Can handle +-Inf, NaN.
    error('repeat_struct_array:Assertion', 'nRepetitions=%g is not an integer.', nRepetitions)
end
if nRepetitions < 1
    error('repeat_struct_array:Assertion', 'Non-positive nRepetitions=%i.', nRepetitions)
end
if ~isfield(StructArray0, numFieldName)
    error('repeat_struct_array:Assertion', 'StructArray0 does not have field numFieldName="%s".', numFieldName)
end
if max([StructArray0.(numFieldName)]) >= nRepetitionOffset
    % NOTE: Assertion is not logically required. Only required from how function is meant to be used
    % (numeric field=timestamp; non-overlapping sequences).
    error('repeat_struct_array:Assertion', 'nRepetitionOffset=%g too short compared to contents of numeric field in struct array.', nRepetitionOffset)
end



StructArray0 = StructArray0(:);   % Convert into column vector.



% Repeat the struct array such as it is (without modifying any field).
StructArray = repmat(StructArray0, nRepetitions, 1);


%=====================================
% Add offset(s) to the selected field
%=====================================
% Create array of offsets, one per StructArray element (not per repetition). Ex: [0,0,0, 10,10,10, 20,20,20, ...].
offsetArray = repelem(nRepetitionOffset*[0:(nRepetitions-1)], numel(StructArray0));   % NOTE: Does not need to create COLUMN array.

% Create array of NEW numeric field values (one per StructArray element). Ex: [0,1,5, 10,11,15, 20,21,25, ...].
numFieldArray = [StructArray.(numFieldName)] + offsetArray;     % NOTE: Does not work for empty StructArray.

% Assign field in StructArray.
tCellArray = num2cell(numFieldArray);
[StructArray.(numFieldName)] = deal(tCellArray{:});

end
