function c_eval(ev_str,sc_list_1,sc_list_2)
%C_EVAL evaluate expression for list of spacecraft
%
% c_eval(ev_str,[sc_list_1],[sc_list_2])
%
% Input:
% ev_str - string to evaluate.
% '?' sign in ev_str is replaced by SC number from sc_list_1
% '!' sign in ev_str is replaced by SC number from sc_list_2
% sc_list - list of SC [optional], 1:4 is assumed when not given
% 
% Example:
% c_eval('R?=r?;C?=R?.^2;',2:4)
% is the same as R2=r2;C2=R2.^2;R3=r3;C3=R3.^2;...
% 
% c_eval('r!r?=irf_abs(r!r?);');
% is the same as
% r1r1=irf_abs(r1r1);r1r2=irf_abs(r1r2);...r4r4=irf_abs(r4r4);
%
% See also IRF_SSUB, EVALIN
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(1,3,nargin))

if nargin==0,
    help c_eval;
elseif nargin==1,
    sc_list_1=1:4;
    sc_list_2=1:4;
elseif nargin==2,
    sc_list_2=1:4;
end

if findstr(ev_str,'?'),
	if findstr(ev_str,'!'),
		for num1=sc_list_1,
			for num2=sc_list_2,
				evalin('caller', irf_ssub(ev_str, num1,num2)),
			end
		end
	else
		for cl_id=sc_list_1, evalin('caller', irf_ssub(ev_str, cl_id)), end
	end
else
	irf_log('fcal','nothing to substitute');
	evalin('caller', ev_str)
end

