function res = rmDoubleIndex(ii)
% rmDoubleIndex(ii)
% remove multiple occurence of index values
%
% $Id$
%

ii = ii(:);
ii = sort(ii);

while 1
	dff = ii(2:end) - ii(1:end-1);
	if isempty(find(dff == 0)), break, end
	ii_nrep = [1; find(dff ~= 0) + 1];
	ii = ii(ii_nrep);
end

res = ii;
