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
% see also C_4_J, IRF_DEC_PARPERP, IRF_LOWPASS
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if nargin~=2 & nargin~=3 & nargin~=4 & nargin~=8 & nargin~=9 & nargin~=10
	error('wrong number of input arguments')
end

args = varargin;

if nargin < 6
	% we have 2 string arguments
	for cl_id=1:4
		ttt = evalin('caller',irf_ssub(args{2},cl_id)); 
		eval(irf_ssub('B? =ttt;',cl_id)); clear ttt
		ttt = evalin('caller',irf_ssub(args{1},cl_id)); 
		eval(irf_ssub('R? =ttt;',cl_id)); clear ttt
	end
	if length(args) > 2, args = args(3:end); 
	else, args = ''; end
else
	% We have 8 arguments
	c_eval('B? = args{4+?};');
	c_eval('R? = args{?};');
	if length(args) > 8, args = args(5:end); 
	else, args = ''; end
end

if length(args)>0
	ref_sc = args{1};
	if length(args) > 1, args = args(2:end); 
	else, args = ''; end
else
	ref_sc = 1;
end
if length(args)>0
	fcut = args{1};
	if length(args) > 1, args = args(2:end); 
	else, args = ''; end
else
	fcut = .7;
end



c_eval('fhz?=50/(B?(51,1)-B?(1,1));')

%irf_lowpass filter B
for j=2:4,c_eval(['B?(:,' num2str(j) ')=irf_lowpass(B?(:,' num2str(j) '),fcut,fhz?);']),end

[jj,divB] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4);

c_eval('[jpar,jperp]=irf_dec_parperp(B?,jj);',ref_sc)

h = c_pl_tx('av_abs(B?)');

axes(h(1))
c_eval('irf_plot(av_abs(B?));',ref_sc)
ylabel(['C' num2str(ref_sc) ' B_{<' num2str(fcut) 'Hz} GSE [nT]'])

axes(h(2))
irf_plot(jpar);
ylabel('j_{||} [A/m^2]')

axes(h(3))
irf_plot(jperp);
ylabel('j_{\perp} GSE [A/m^2]')

axes(h(4))
irf_plot(divB);
ylabel('div(B) [A/m^2]')

for j=1:4,set(h(j),'YLim',get(h(j),'YLim')*.99),end

for j=1:3,xlabel(h(j),''),set(h(j),'XTickLabel',''),end

legend(h(1),'x','y','z','tot','Location','NorthEastOutside')
legend(h(3),'x','y','z','Location','NorthEastOutside')
if nargout>0, hout = h; end
if nargout>1, jout = jj; end
if nargout>2, divBout = divB; end
if nargout>3, jparout = jpar; end
if nargout>4, jperpout = jperp; end
