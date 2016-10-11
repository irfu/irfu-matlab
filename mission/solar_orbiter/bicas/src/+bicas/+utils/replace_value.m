% x = replace_value(x, before, after)   Replace every occurrence of a specific value (incl. NaN) in an array.
%
% General-purpose routine, but created specifically for the purpose of converting fill/pad values<--->NaN when
% reading/writing CDF files.
%
% NOTE: Handles NaN as any other value.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-05
%
function x = replace_value(x, before, after)
	if isnan(before)
		i = isnan(before);
	else
		i = (x==before);
	end
	x(i) = after;
end	
