function [y]=av_add(c1,x1,c2,x2)
% function [y]=av_add(c1,x1,c2,x2)
% estimates y=c1*x1+c2*x2;
% c1,c2 - scalars
% x1,x2 - time series with column one being time
global AV_DEBUG 
if isempty(AV_DEBUG), debug=0; else debug=AV_DEBUG;end

if size(x1,2) ~= size(x2,2), disp('ERROR: Vectors not of the same dimension');return;end
if size(x1,1) ~= size(x2,1), 
 if debug == 1,disp('interpolating x2 to x1.av_add(c1,x1,c2,x2)');end
 x2=av_interp(x2,x1);
end

y=x1;
y(:,2:end)=c1*x1(:,2:end)+c2*x2(:,2:end);
