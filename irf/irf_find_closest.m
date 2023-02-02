function [t1new,t2new,ind1new,ind2new]=irf_find_closest(t1,t2)
%IRF_FIND_CLOSEST  Finds pairs that are closest to each other in two timeseries
%
% [t1new,t2new,ind1new,ind2new]=irf_find_closest(t1,t2);
%
% t1,t2 - vector with time instants
% t1new,t2new - the identified time instants that are closest each other
%               t1new(1)<->t2new(1)
%               t1new(2)<->t2new(2)
%               ...
%

t1_orig=t1;t2_orig=t2;
flag=1;
while flag
  flag_t1=zeros(size(t1));
  ind=interp1(t1,1:length(t1),t2,'nearest','extrap');
  flag_t1(ind)=1;
  flag_t2=zeros(size(t2));
  ind=interp1(t2,1:length(t2),t1,'nearest','extrap');
  flag_t2(ind)=1;
  ind_zeros_t1=find(flag_t1 == 0);
  ind_zeros_t2=find(flag_t2 == 0);
  if ~isempty(ind_zeros_t1)
    irf_log('proc',['Throwing away t1(' num2str(ind_zeros_t1(:)') ')']);
    t1(ind_zeros_t1)=[];
  elseif ~isempty(ind_zeros_t2)
    irf_log('proc',['Throwing away t2(' num2str(ind_zeros_t2(:)') ')']);
    t2(ind_zeros_t2)=[];
  else
    flag=0;
    break;
  end
end

t1new=t1;
t2new=t2;
ind1new=interp1(t1_orig,1:length(t1_orig),t1new,'nearest');
ind2new=interp1(t2_orig,1:length(t2_orig),t2new,'nearest');
