function dnum = todatenum(tt2000Obj)
%TODATENUM Convert a cdftt2000 object to a UTC string.
%   N = TODATENUM(CDFEPOCHOBJ) retrieves the cdftt2000 object values
%
%   See also CDFTT2000, CDFWRITE.

s = size(tt2000Obj);
% 
dnum = [tt2000Obj.date];
if ~isempty(dnum)
    dnum = reshape(dnum, s);
end

