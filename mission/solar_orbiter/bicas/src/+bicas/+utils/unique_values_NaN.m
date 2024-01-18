function uniqueValues = unique_values_NaN(A)
    % Return unique values in an arbitrary array, treating +Inf, -Inf, and NaN as equal to themselves.
    % (MATLAB's "unique" function does not do this.)
    %
    % ARGUMENTS
    % =========
    % A : Numeric array of any dimensionality.
    %
    %
    % RETURN VALUE
    % ============
    % uniqueValues : 1D array of unique values. Sorted.
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    % First created 2016-10-11
    %



    % IMPLEMENTATION NOTE: The MATLAB command "unique" has special behaviour which must be taken into account:
    % 1) Inf and -Inf are treated as equal to themselves.
    % 2) NaN is treated as if it is NOT equal itself. ==> Can thus return multiple instances of NaN.
    % 3) Always puts NaN at the end of the resulting vector (one or multiple NaN).



    uniqueValues = unique(A);

    % Remove all NaN unless it is found in the last component (save one legitimate occurrence of NaN, if there is any).
    % NOTE: Does work for empty matrices.
    uniqueValues(isnan(uniqueValues(1:end-1))) = [];

end
