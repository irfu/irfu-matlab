%
% Convert CDF TT2000 to MATLAB date vector.
%
%
% ARGUMENTS
% =========
% tt2000 : Scalar. Does not have to be int64.
%
%
% RETURN VALUE
% ============
% dateVec : MATLAB date vector 1x6. See MATLAB function "datevec".
%
function dateVec = tt2000_to_datevec(tt2000)
    assert(isscalar(tt2000))
    
    % OUT = spdfbreakdowntt2000(tt2000) returns the UTC date/time from CDF TT2000
    % time. OUT is an array with each row having nine (9) numerical values
    % for year, month, day, hour, minute, second, millisecond, microsecond
    % and nanosecond.
    v = spdfbreakdowntt2000(tt2000);
    
    dateVec = [v(1:5), v(6) + 1e-3*v(7) + 1e-6*v(8) + 1e-9*v(9)];
end