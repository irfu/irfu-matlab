function [y]=irf_lmn(x,b,n,flag)
% function [y]=irf_lmn(x,b,n,flag)
% function [y]=irf_lmn(x,L,M,N)
%
% x - input vector consists of columns [time xx xy xz]
% b - magnetic field, [time bx by bz]
% n - normal direction [nx ny nz]
% flag - method of estimate, 'L' or 'N', see below
% L,M,N - [Lx Ly Lz],[Mx ..],..
%         if L,M,N are not orthogonal then L=(NxL)xN, M=NxL
%
% y - output vector in LMN coordinates [time y_L y_M y_N] or [y_L y_M y_N], depending on x
%
% flag='L'
%
%  L - along B
%  N - closest to n
%  M - NxL
% the coordinate system follows B and thus is not stationary if B changes
%
% flag='N'
%
%  N - along n
%  L - closest to the mean direction of B
%  M - NxL
%

irf_log('proc','WARNING: The sign of the M-component has been changed since 2012.03.01.')

if nargin ==3, flag_case='L'; end
if (nargin ==4)
  if ischar(flag), flag_case=flag;
  elseif isnumeric(flag),  L=b; M=n; N=flag; flag_case='LMN';
  end
end

switch flag_case
  case 'L'
    be=irf_resamp(b,x);

    nl=irf_norm(be); % along the B
    nn=irf_norm(irf_cross(be,irf_cross(n,be))); % closest to given vn vector
    nm=irf_cross(nn,nl); % in (n x b) direction

  case 'N'
    nn=irf_norm(n);
    nm=irf_norm(irf_cross(mean(b),nn));
    nl=irf_cross(nm,nn);

  case 'LMN'
    if abs(dot(cross(L,M),N)-norm(L)*norm(M)*norm(N))>1e-5
      disp('irf_lmn: L,M,N does not satisfy M=LxN,N=MxL! using L=(NxL)xN; M=NxL;');
      L=cross(cross(N,L),N);
      M=cross(N,L);
    end
    nl=irf_norm(L);nm=irf_norm(M);nn=irf_norm(N);
end

% estimate e in new coordinates
xn=irf_dot(x,nn,1);
xl=irf_dot(x,nl,1);
xm=irf_dot(x,nm,1);

switch size(x,2)
  case 3 % no time column
    y=x; y(:,1)=xl; y(:,2)=xm; y(:,3)=xn;
  otherwise % [time xx xy xz] or [time xx xy xz abs(x)]
    y=x; y(:,2)=xl; y(:,3)=xm; y(:,4)=xn;
end
