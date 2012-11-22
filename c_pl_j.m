function [hout,jout,divBout,jparout,jperpout]=c_pl_j(varargin)
%C_PL_J   Calculate and plot current using 4 spacecraft technique
%
% [h,j,divB,jpar,jperp]=c_pl_j(r1,r2,r3,r4,b1,b2,b3,b4,[ref_sc],[f_cut])
% [h,j,divB,jpar,jperp]=c_pl_j('R?','B?',[ref_sc],[f_cut])
%
% Calculates and plots current from using 4 spacecraft technique.
%
% Input:
%   r1..4 : positions
%   b1..4 : magnetic fields
%   ref_sc : reference spacecraft [OPTIONAL]
%   f_cut : cutoff frequency fot the irf_lowpass filter [OPTIONAL]
%
% Output [OPTIONAL]:
%   h : axes handles
%   j : current vector
%   divB : divergence B as the error estimate
%   jpar,jperp : parallel and perpendicular components of J
%
% Example:
%   h = c_pl_j('R?','B?',1,5);
%
% See also C_4_J, IRF_JZ, IRF_DEC_PARPERP, IRF_LOWPASS
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

B1 = [];

if nargin==0, 
  help c_pl_j;return;
end
 if nargin~=2 && nargin~=3 && nargin~=4 && nargin~=8 && nargin~=9 && nargin~=10
	error('wrong number of input arguments')
end
 
args = varargin;

if nargin < 6
	% we have 2 string arguments
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(args{2},cl_id));  %#ok<NASGU>
		eval(irf_ssub('B? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(args{1},cl_id));  %#ok<NASGU>
		eval(irf_ssub('R? =ttt;',cl_id)); clear ttt
	end
	if length(args) > 2, args = args(3:end); 
	else args = ''; end
else
	% We have 8 arguments
	c_eval('B? = args{4+?};');
	c_eval('R? = args{?};');
	if length(args) > 8, args = args(5:end); 
	else args = ''; end
end

if ~isempty(args)
	ref_sc = args{1};
	if length(args) > 1, args = args(2:end); 
	else args = ''; end
else
	ref_sc = 1;
end
if ~isempty(args)
	fcut = args{1};
else
	fcut = .7;
end



c_eval('fhz?=50/(B?(51,1)-B?(1,1));')

%irf_lowpass filter B
for j=2:4,c_eval(['B?(:,' num2str(j) ')=irf_lowpass(B?(:,' num2str(j) '),fcut,fhz?);']),end

[jj,divB,B,jxB] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4); %#ok<ASGLU>

jpar = []; jperp = [];
c_eval('[jpar,jperp]=irf_dec_parperp(B?,jj);',ref_sc)

NPLOTS = 6;
h = 1:NPLOTS;

h(1) = irf_subplot(NPLOTS, 1, -1);
c_eval('irf_plot(irf_abs(B?));',ref_sc)
ylabel([' B_{<' num2str(fcut) 'Hz} GSE [nT]'])

h(2) = irf_subplot(NPLOTS, 1, -2);
jpar(:,2:end) = jpar(:,2:end)*1e9;
irf_plot(jpar);
ylabel('j_{||} [nA/m^2]')

h(3) = irf_subplot(NPLOTS, 1, -3);
jperp(:,2:end) = jperp(:,2:end)*1e9;
irf_plot(jperp);
ylabel('j_{\perp} GSE [nA/m^2]')

h(4) = irf_subplot(NPLOTS, 1, -4);
irf_plot(jxB);
ylabel('jxB [A/m^2 T]')

h(5) = irf_subplot(NPLOTS, 1, -5);
irf_plot(irf_integrate(jxB));
ylabel('\int jxB [A/m^2 T s]')

h(6) = irf_subplot(NPLOTS, 1, -6);
divB(:,2:end) = divB(:,2:end)*1e9;
irf_plot(divB);
ylabel('div(B) [nA/m^2]')

for j=1:6,set(h(j),'YLim',get(h(j),'YLim')*.99),end

for j=1:5,xlabel(h(j),''),set(h(j),'XTickLabel',''),end

irf_zoom(h,'x',B1([1 end],1)')

if nargout>0, hout = h; end
if nargout>1, jout = jj; end
if nargout>2, divBout = divB; end
if nargout>3, jparout = jpar; end
if nargout>4, jperpout = jperp; end
