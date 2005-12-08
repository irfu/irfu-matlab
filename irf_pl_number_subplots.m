function ht=irf_pl_number_subplots(h)
% ht=irf_pl_number_subplots(h)
% 
% number subplots by putting in each subplot 'a', 'b', 'c', etc.
%
% h is handles to subplots

h=h(:);
abc='a':'z';
for j=1:length(h)
  axes(h(j));
  ht(j)=irf_pl_info(abc(j),gca,[0.01,.8]);
  set(ht(j),'fontsize',12);
end
