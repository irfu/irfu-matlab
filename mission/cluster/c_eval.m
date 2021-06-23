function c_eval(ev_str,sc_list_1,sc_list_2)
%C_EVAL evaluate expression for lists of values
%	most often used when different variables need to be constructed
%	for different spacecraft.
%
% c_eval(ev_str,[sc_list_1],[sc_list_2])
%
% Input:
% ev_str - string to evaluate.
% '?' sign in ev_str is replaced by values from list sc_list_1
% '!' sign in ev_str is replaced by values from list sc_list_2
% sc_list_1,sc_list_2 - list of values, when omitted then put to [1 2 3 4].
%           sc_list can be also cell vector, e.g. {'a','b','c'}
%
% Example:
%   c_eval('R?=r?;C?=R?.^2;',2:4)
%          is the same as R2=r2;C2=R2.^2;R3=r3;C3=R3.^2;...
%
%   c_eval('r!r?=irf_abs(r!r?);');
%          is the same as
%          r1r1=irf_abs(r1r1);r1r2=irf_abs(r1r2);...r4r4=irf_abs(r4r4);
%
%	c_eval('a?=2;',{'a','b','c'});
%			is the same as aa=2;ab=2;ac=2;
%
% See also IRF_SSUB, EVALIN

if nargin==0
  help c_eval;
elseif nargin==1
  sc_list_1=1:4;
  sc_list_2=1:4;
elseif nargin==2
  sc_list_2=1:4;
elseif nargin > 3
  irf_log('fcal','cannot be more than 3 input arguments')
  return
end

if strfind(ev_str,'?') %#ok<STRIFCND>
  if strfind(ev_str,'!') %#ok<STRIFCND>
    for num1=sc_list_1
      for num2=sc_list_2
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

