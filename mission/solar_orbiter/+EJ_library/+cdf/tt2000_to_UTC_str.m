%
% Convert CDF tt2000 value to UTC string with nanoseconds.
%
% NOTE: spdfparsett2000 seems to work well as an inverse without wrapper.
%
%
% ARGUMENT
% ========
% tt2000 : Scalar numeric value.
%
%
% RETURN VALUE
% ============
% utcStr : Example: 2020-04-01T01:23:45.678901234
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-04-03.
%
function utcStr = tt2000_to_UTC_str(tt2000)    
    % TODO-DECISION: How handle various needs for formats? Rounding, truncation?
    % PROPOSAL: Assertions on argument being int64 as they are in CDF?
    % PROPOSAL: Convert array.
    %   NOTE: Return value must be cell array.
    %   NOTE: Special case for empty array.
    %
    % NOTE: Should be analogous to any inverted conversion function.
    
    % NOTE: spdfbreakdowntt2000 (used here) is the inverse to spdfparsett2000.
    %
    %  spdfbreakdowntt2000 converts the CDF TT2000 time, nanoseconds since
    %               2000-01-01 12:00:00 to UTC date/time.
    %
    %     OUT = spdfbreakdowntt2000(tt2000) returns the UTC date/time from CDF TT2000
    %     time. OUT is an array with each row having nine (9) numerical values
    %     for year, month, day, hour, minute, second, millisecond, microsecond
    %     and nanosecond.    
    
    % NOTE: spdfbreakdowntt2000 can handle Nx1 arrays, where N>=1.    
    assert(isscalar(tt2000), 'Illegal argument tt2000. Must be scalar.')
    
    
    %v = spdfbreakdowntt2000(tt2000);
    %utcStr = sprintf('%04i-%02i-%02iT%02i:%02i:%02i.%03i%03i%03i', v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9));
    
    dateVec = EJ_library.cdf.tt2000_to_datevec(tt2000);
    utcStr = sprintf('%04i-%02i-%02iT%02i:%02i:%02.9g', dateVec(1:6));
end
