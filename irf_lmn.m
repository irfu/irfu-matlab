function [y]=av_lmn(x,b,n,flag)
% function [y]=av_lmn(x,b,n,flag)
% function [y]=av_lmn(x,L,M,N)
%
% x - input vector consists of columns [time y_N y_M y_L]
% b - magnetic field, [time bx by bz]
% n - normal direction [nx ny nz]
% flag - method of estimate, 'B' or 'N', see below
% L,M,N - [Lx Ly Lz],[Mx ..],..
%         if L,M,N are not orthogonal then L=(NxL)xN, M=LxN
%
% y - output vector in LMN coordinates [time y_L y_M y_N] or [y_L y_M y_N], depending on x
%
% flag='L'
%
%  L - along B
%  N - closest to n
%  M - LxN
% the coordinate system follows B and thus is not stationary if B changes
%
% flag='N'
%
%  N - along n
%  L - closest to the mean direction of B
%  M - LxN
%

if nargin ==3, flag_case='L';end
if (nargin ==4) & isstr(flag), flag_case=flag;end
if (nargin ==4) & isnumeric(flag), L=b;M=n;N=flag;flag_case='LMN';end

switch flag_case
case 'L'
  be=irf_resamp(b,x);

  nl=irf_norm(be); % along the B
  nn=irf_norm(irf_cross(be,irf_cross(n,be))); % closest to given vn vector
  nm=irf_cross(nl,nn); % in (b x n) direction

case 'N'
  nn=irf_norm(n);
  nm=irf_norm(irf_cross(mean(b),nn));
  nl=irf_cross(nn,nm);

case 'LMN'
  if abs(dot(cross(M,L),N)-norm(L)*norm(M)*norm(N))>1e-5, 
    disp('av_lmn: L,M,N does not satisfy M=LxN,N=MxL! using L=(NxL)xN; M=LxN;');
    L=cross(cross(N,L),N);
    M=cross(L,N);
  end
  nl=irf_norm(L);nm=irf_norm(M);nn=irf_norm(N);
end

% estimate e in new coordinates
  xn=irf_dot(x,nn,1);
  xl=irf_dot(x,nl,1);
  xm=irf_dot(x,nm,1);

  y=x; y(:,end-2)=xl;y(:,end-1)=xm;y(:,end)=xn;
