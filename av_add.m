function [y]=av_add(c1,x1,c2,x2)
% function [y]=av_add(c1,x1,c2,x2)
% estimates y=c1*x1+c2*x2;
% c1,c2 - scalars
% x1,x2 - time series with column one being time
global AV_DEBUG; if isempty(AV_DEBUG), debug=0; else debug=AV_DEBUG;end

if size(x1,2)>size(x2,2),
  disp('av_add(c1,x1,c2,x2) WARNING: x1 has more columns than x2');
  if size(x2,2)==2, for j=3:size(x1,2), x2(:,j)=x2(:,2); end % assume x2 is time series of scalar
  elseif size(x2,2)==1, ones(size(x1))*x2(1,1);x2(:,1)=x1(:,1); % add x2(1,1) to all x1 
  elseif size(x2,2)==size(x1,2)-1, x2=[x1(1,1) x2(1,:)]; % use only first row of x2, taking time from x1
  else, disp('av_add() ERROR: could not make inteligent guess what you are meaning.');return; 
  end
elseif size(x2,2)>size(x1,2),
  disp('av_add(c1,x1,c2,x2) WARNING: x2 has more columns than x1');
  if size(x1,2)==2, for j=3:size(x2,2), x1(:,j)=x1(:,2); end % assume x2 is time series of scalar
  elseif size(x1,2)==1, ones(size(x2))*x1(1,1);x1(:,1)=x2(:,1); % add x2(1,1) to all x1 
  elseif size(x1,2)==size(x2,2)-1, x1=[x2(1,1) x1(1,:)]; % use only first row of x2, taking time from x1
  else, disp('av_add() ERROR: could not make inteligent guess what you are meaning.');return; 
  end
end
if size(x1,1) ~= size(x2,1), 
 if debug == 1,disp('interpolating x2 to x1. av_add(c1,x1,c2,x2)');end
 x2=av_interp(x2,x1);
end

y=x1;
y(:,2:end)=c1*x1(:,2:end)+c2*x2(:,2:end);
