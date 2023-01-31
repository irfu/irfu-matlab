function [y]=irf_add(c1,x1,c2,x2)

%IRF_ADD   add two vectors
%
% [y]=irf_add(c1,x1,c2,x2)
% estimates y=c1*x1+c2*x2;
% c1,c2 - scalars
% x1,x2 - time series with column one being time
%

if isa(x1,'TSeries') % assumes x2 is TSeries and no. columns are the same
  x2 = x2.resample(x1);
  y = x1;
  y.data = c1*x1.data+c2*x2.data;
else
  if size(x1,2)>size(x2,2)
    irf_log('proc','WARNING: x1 has more columns than x2');
    if size(x2,2)==2, for j=3:size(x1,2), x2(:,j)=x2(:,2); end % assume x2 is time series of scalar
    elseif size(x2,2)==1, x2(:,2)=ones(size(x1))*x2(1,1);x2(:,1)=x1(:,1); % add x2(1,1) to all x1
    elseif size(x2,2)==size(x1,2)-1, x2=[x1(1,1) x2(1,:)]; % use only first row of x2, taking time from x1
    else
      irf_log('proc','ERROR: could not make inteligent guess what you are meaning.');return;
    end
    
  elseif size(x2,2)>size(x1,2)
    irf_log('proc','WARNING: x2 has more columns than x1');
    if size(x1,2)==2, for j=3:size(x2,2), x1(:,j)=x1(:,2); end % assume x2 is time series of scalar
    elseif size(x1,2)==1, x1(:,2)=ones(size(x2))*x1(1,1);x1(:,1)=x2(:,1); % add x2(1,1) to all x1
    elseif size(x1,2)==size(x2,2)-1, x1=[x2(1,1) x1(1,:)]; % use only first row of x2, taking time from x1
    else
      irf_log('proc','irf_add() ERROR: could not make inteligent guess what you are meaning.');return;
    end
  end
  
  if size(x1,1) ~= size(x2,1)
    irf_log('proc','interpolating x2 to x1.');
    x2=irf_resamp(x2,x1);
  end
  
  y=x1;
  
  y(:,2:end)=c1*x1(:,2:end)+c2*x2(:,2:end);
end