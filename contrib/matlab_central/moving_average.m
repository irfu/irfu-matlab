function [Y,Nsum] = moving_average(X,F,DIM)
%MOVING_AVERAGE   Smooths a vector through the moving average method.
%
%   Syntax:
%     [Y,Nsum] = moving_average(X,F,DIM);
%
%   Input:
%     X   - Vector or matrix of finite elements.
%     F   - Window semi-length. A positive scalar (default 0).
%     DIM - If DIM=1: smooths the columns (default); elseif DIM=2 the rows.
%
%   Output:
%     Y    - Smoothed X elements.
%     Nsum - Number of not NaN's elements that fixed on the moving window.
%            Provided to get a sum instead of a mean: Y.*Nsum.
%
%   Description:
%     Quickly smooths the vector X by averaging each element along with the
%     2*F elements at its sides. The elements at the ends are also averaged
%     but the extrems are left intact. With the windows size defined in
%     this way, the filter has zero phase.
%
%   Example:
%      x = 2*pi*linspace(-1,1)'; 
%      yn = cos(x) + 0.25 - 0.5*rand(size(x)); 
%      ys = moving_average(yn,4);
%      plot(x,[yn ys]), legend('noisy','smooth',4), axis tight
%
%   See also FILTER, RECTWIN and MOVING_AVERAGE2, NANMOVING_AVERAGE,
%   NANMOVING_AVERAGE2 by Carlos Vargas and RUNMEAN by Jos van der Geest.

% Copyright 2006-2008  Carlos Vargas, nubeobscura@hotmail.com
%	$Revision$  $Date$

%   Written by
%   M. in S. Carlos Adrián Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico,  march 2008
%
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

% Updated source: https://www.mathworks.com/matlabcentral/fileexchange/12276-moving-average-v3-1--mar-2008-/

%   2008 Mar. Use CUMSUM as RUNMEAN by Jos van der Geest, no more
%   subfunctions.

% LICENSE UNSPECIFIED

%% Error checking:
if ~nargin
 error('Moving_average:Inputs','There are no inputs.')
elseif nargin<2 || isempty(F)
 F = 0;
end
if F==0
 Y = X;
 return
end
F = round(F);
ndim = ndims(X);
if (ndim ~= 2)
 error('Moving_average:Inputs','Input is not a vector or matrix.')
end
[N,M] = size(X);
if nargin<3 || isempty(DIM)
 DIM = 1;
 if N == 1
  DIM = 2;
 end
end
if DIM == 2
 X = X.';
 [N,M] = size(X);
end
if 2*F+1>N
 warning('Moving_average:Inputs',... % bug fixed 06 Mar 2008
  'Window size must be less or equal as the number of elements.')
 Y = X;
 if DIM == 2
  Y = Y.';
 end
 return
end

%% Window width
Wwidth = 2*F + 1;  

%% Smooth the edges but with the first and last element intact
F2 = Wwidth - 2;
Nsumedge = repmat((1:2:F2)',1,M);
Y1 =        X(     1:F2,:);
Y2 = flipud(X(N-F2+1:N ,:));
Y1 = cumsum(Y1,1);     
Y2 = cumsum(Y2,1);
Y1 = Y1(1:2:F2,:)./Nsumedge;            
Y2 = Y2(1:2:F2,:)./Nsumedge;

%% Recursive moving average method 
% With CUMSUM trick copied from RUNMEAN by Jos van der Geest (12 mar 2008)
Y = [zeros(F+1,M); X; zeros(F,M)];
Y = cumsum(Y,1);
Y = Y(Wwidth+1:end,:)-Y(1:end-Wwidth,:);
Y = Y/Wwidth;

%% Sets the smoothed edges:
Y(    1:F,:) =        Y1;
Y(N-F+1:N,:) = flipud(Y2);

%% Get the number of elements that were averaged for each element: 
if nargout == 2
 Nsum = repmat(Wwidth,size(Y));
 Nsum(    1:F,:) = Nsumedge;
 Nsum(N-F+1:N,:) = flipud(Nsumedge);
 if DIM ==2
  Nsum = Nsum.';
 end
end

%% Return the correct size:
if DIM == 2
 Y = Y.';
end 

%% % Recursive moving average code before Jos trick:
% Y = X;          
% Y(F+1,:) = sum(X(1:Wwidth,:),1);    
% for n = F+2:N-F
%  Y(n,:) = sum([Y(n-1,:); X(n+F,:); -X(n-F-1,:)],1); 
% end
% Y = Y/Wwidth;
