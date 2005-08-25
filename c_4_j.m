function [j,divB,B,jxB]=c_4_j(r1,r2,r3,r4,b1,b2,b3,b4)
%C_4_J  Calculate current from using 4 spacecraft technique
%  in addition one can obtain average magnetic field and jxB values
%
%  [j,divB,B,jxB] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4)
%  [j,divB,B,jxB] = c_4_j('R?','B?')
%  Estimate also divergence B as the error estimate
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
%  B      - average magnetic field, sampled at b1 time steps [nT]
%  jxB    - j x B force [T A]
%
%   See also C_4_K
%
%  Reference: ISSI book  Eq. 14.16, 14.17
%
% $Id$

if nargin~=8 & nargin~=2
	disp('Too few parameters. See usage:');
	help c_4_j;     
	return
end
if nargin==2
	rs = r1;
	bs = r2;
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(rs,cl_id)); 
		eval(irf_ssub('r? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(bs,cl_id)); 
		eval(irf_ssub('b? =ttt;',cl_id)); clear ttt
	end
	clear bs rs
end

%%%%%%%%%%%%%%%% Estimate first reciprical coordinates %%%%%%%%%%%%%%
%
% because usually r1..r4 is of less time resolution, it is more
% computer friendly first calculate k1..k4 and only after interpolate
% and not the other way around
for ic=1:4,eval(irf_ssub('R?=irf_resamp(r?,r1,''spline'');',ic)),end
[k1,k2,k3,k4]=c_4_k(R1,R2,R3,R4);

%%%%%%%%%%%%%%%% Do interpolation to b1 time series %%%%%%%%%%%%%%%%%%%%%%
for ic=1:4,eval(irf_ssub('B?=irf_resamp(b?,b1);',ic)),end
B=0.25*B1+0.25*B2+0.25*B3+0.25*B4; % estimate mean value of B
for ic=1:4,eval(irf_ssub('K?=irf_resamp(k?,b1);',ic)),end

% initialize matrix j and divB with right time column
j=B1(:,1:4);j(:,2:4)=0;
divB=B1(:,1:2);divB(:,2)=0;

% Calculate j and divB
for ic=1:4, eval(irf_ssub(  'divB(:,2)=divB(:,2)+dot(K?(:,2:4),B?(:,2:4),2);'  ,ic));   end
divB(:,2)=divB(:,2)/1.0e3*1e-9/(4*pi*1e-7); % to get right units

for ic=1:4, eval(irf_ssub('j(:,2:4)=j(:,2:4)+cross(K?(:,2:4),B?(:,2:4),2);',ic));   end
j(:,2:4)=j(:,2:4)/1.0e3*1e-9/(4*pi*1e-7);   % to get right units [A]
jxB=irf_tappl(irf_cross(j,B),'*1e-9'); % to get units [T A]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0&size(B1,1)==1,
       strj=['j= ' num2str(norm(j(1,2:4)),3) ' [ ' num2str(j(1,2:4)/norm(j(1,2:4)),' %5.2f') '] A '];
       strdivB=['divB= ' num2str(divB(1,2),3) '] A '];
       disp(strj);disp(strdivB);
end

