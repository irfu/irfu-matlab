function c_eval(ev_str,sc_list)
%c_eval evaluate expression for list of spacecraft
% c_eval(ev_str,[sc_list])
%
% Input:
% ev_str - string to evaluate.
% '?' sign in ev_str is replaced by SC number from sc_list
% sc_list - list of SC [optional], 1:4 is assumed when not given
% 
% Example:
% c_eval('R?=r?;C?=R?.^2;',2:4)
%
% is the same as R2=r2;C2=R2.^2;R3=r3;C3=R3.^2;...
%
% $Id$
%
% See also av_ssub, evalin

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(1,2,nargin))

if nargin<2, sc_list=1:4; end

for cl_id=sc_list, evalin('caller', av_ssub(ev_str, cl_id)), end
