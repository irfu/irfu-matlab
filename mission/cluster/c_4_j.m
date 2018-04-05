function [j,divB,B,jxB,divTshear,divPb]=c_4_j(r1,r2,r3,r4,b1,b2,b3,b4) %#ok<INUSL,INUSD>
%C_4_J  Calculate current from using 4 spacecraft technique
%  in addition one can obtain average magnetic field and jxB values
%
%  [j,divB,B,jxB,curvature,divTshear,divPb] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4)
%  [j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?')
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
%         column 2-4   current, units A/m^2
%  divB   column 1     time
%         column 2     div(B)/mu0, units A/m^2
%  B      - average magnetic field, sampled at b1 time steps [nT]
%  jxB    - j x B force [T A]
%         - jxB=(1/muo) ( (B div)B + grad (B^2/2) )= divTshear+divPb
%  divTshear = (1/muo) (B div) B.  the part of the divergence of stress 
%                                   associated with curvature units [T A/m^2]
%  divPb = (1/muo) grad(B^2/2). gradient of magnetic pressure
% 
%   See also C_4_K
%
%  Reference: ISSI book  Eq. 14.16, 14.17

% TODO fix that it works for vector inputs without time column!

if nargin~=8 && nargin~=2
	disp('Too few parameters. See usage:');
	help c_4_j;     
	return
end
if nargin==2
	rs = r1;
	bs = r2;
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(rs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('r? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(bs,cl_id)); %#ok<NASGU>
		eval(irf_ssub('b? =ttt;',cl_id)); clear ttt
	end
	clear bs rs
end

useTSeries = 0;
if isa(r1,'TSeries')
    useTSeries = 1;
end

% Estimate divB/mu0. unit is A/m2
[divB,B]=c_4_grad('r?','b?','div');
divB=irf_tappl(divB,'/1.0e3*1e-9/(4*pi*1e-7)'); % to get right units why 

% estimate current j [A/m2]
curl_B=c_4_grad('r?','b?','curl');
j=irf_tappl(curl_B,'/1.0e3*1e-9/(4*pi*1e-7)');   % to get right units [A/m2]

% estimate jxB force [T A/m2]
if useTSeries
    jxB=irf_tappl(cross(j,B),'*1e-9'); % to get units [T A/m2]
else
    jxB=irf_tappl(irf_cross(j,B),'*1e-9'); % to get units [T A/m2]
end

% estimate divTshear = (1/muo) (B*div)B [T A/m2]
BdivB=c_4_grad('r?','b?','bdivb');
divTshear=irf_tappl(BdivB,'/1e3*1e-9*1e-9/(4*pi*1e-7)');

% estimate divPb = (1/muo) grad (B^2/2) = divTshear-jxB
divPb=irf_add(-1,jxB,1,divTshear);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0 && size(b1,1)==1
       strj=['j= ' num2str(norm(j(1,2:4)),3) ' [ ' num2str(j(1,2:4)/norm(j(1,2:4)),' %5.2f') '] A/m^2 '];
       strdivB=['divB/mu0 = ' num2str(divB(1,2),3) ' A/m^2 '];
       disp(strj);disp(strdivB);
end

