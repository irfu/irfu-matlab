function res = rmDoubleIndex(ii)
% rmDoubleIndex(ii)
% remove multiple occurence of index values
%
% $Id$
%
warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'irf_rm_double_idx')

ii = ii(:);
ii = sort(ii);

while 1
	dff = ii(2:end) - ii(1:end-1);
	if isempty(find(dff == 0)), break, end
	ii_nrep = [1; find(dff ~= 0) + 1];
	ii = ii(ii_nrep);
end

res = ii;
