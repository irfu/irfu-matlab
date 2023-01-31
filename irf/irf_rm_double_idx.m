function res = irf_rm_double_idx(ii)
%IRF_RM_DOUBLE_IDX   Remove multiple occurence of index values
%
% res_idx = irf_rm_double_idx(ii)
%

ii = ii(:);
ii = sort(ii);

while 1
  dff = ii(2:end) - ii(1:end-1);
  if isempty(find(dff == 0, 1)), break, end
  ii_nrep = [1; find(dff ~= 0) + 1];
  ii = ii(ii_nrep);
end

res = ii;
