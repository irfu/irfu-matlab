function [j,divB]=c_4_j(r1,r2,r3,r4,b1,b2,b3,b4)
%C_4_J calculate current from using 4 spacecraft technique
%
%  [j,divB]=c_4_j(r1,r2,r3,r4,b1,b2,b3,b4)
%  estimates also divergence B as the error estimate
%
%  [j]=c_4_j(r1,r2,r3,r4,b1,b2,b3,b4)  Calculates only current
%
%  r1..r4 are row vectors
%         column 1     is time
%         column 2-4   are satellite position in km
%  b1..b4 are row vectors
%         column 1     is time b1 time is used for all interpolation of r1..r4 and b2..b4
%         column 2-4   is magnetic field components in nT
%  j      is row vector,
%         column 1     time
%         column 2-4   current, units A
%  divB   column 1     time
%         column 2     div(B)/mu0, units A
%
%   See also C_4_K
%
%  Reference: ISSI book  Eq. 14.16, 14.17
%
% $Id$

if nargin<8;    disp('Too few parameters. See usage:');help c_4_j;     return;end

%%%%%%%%%%%%%%%% Estimate first reciprical coordinates %%%%%%%%%%%%%%
%
% because usually r1..r4 is of less time resolution, it is more
% computer friendly first calculate k1..k4 and only after interpolate
% and not the other way around
for ic=1:4,eval(irf_ssub('R?=av_interp(r?,r1,''spline'');',ic)),end
[k1,k2,k3,k4]=c_4_k(R1,R2,R3,R4);

%%%%%%%%%%%%%%%% Do interpolation to b1 time series %%%%%%%%%%%%%%%%%%%%%%
for ic=1:4,eval(irf_ssub('B?=av_interp(b?,b1);',ic)),end
for ic=1:4,eval(irf_ssub('K?=av_interp(k?,b1);',ic)),end

% initialize matrix j and divB with right time column
j=B1(:,1:4);j(:,2:4)=0;
divB=B1(:,1:2);divB(:,2)=0;

% Calculate j and divB
for ic=1:4, eval(irf_ssub(  'divB(:,2)=divB(:,2)+dot(K?(:,2:4),B?(:,2:4),2);'  ,ic));   end
divB(:,2)=divB(:,2)/1.0e3*1e-9/(4*pi*1e-7); % to get right units

for ic=1:4, eval(irf_ssub('j(:,2:4)=j(:,2:4)+cross(K?(:,2:4),B?(:,2:4),2);',ic));   end
j(:,2:4)=j(:,2:4)/1.0e3*1e-9/(4*pi*1e-7);   % to get right units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0&size(B1,1)==1,
       strj=['j= ' num2str(norm(j(1,2:4)),3) ' [ ' num2str(j(1,2:4)/norm(j(1,2:4)),' %5.2f') '] A '];
       strdivB=['divB= ' num2str(divB(1,2),3) '] A '];
       disp(strj);disp(strdivB);
end

