function y=av_q(question,variable_name_of_old_value, default_value);
% function y=av_q(question,variable_name_of_old_value, default_value);
% example:
%   y=av_q('How many 1/0? [%]>','y',0)
%   y=av_q('Large? yes/no [%]>','y','yes')
%   y=av_q('How much? [%]>','y',10)
% see av_ask

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'irf_ask')

if nargin < 3, disp('ERROR using av_q');help av_q;return;end
if nargin == 3, def=default_value;else def=[];end

if evalin('caller',['exist(''' variable_name_of_old_value ''')']),
  defvalue=evalin('caller', [variable_name_of_old_value]);
else
  defvalue=def;
end
if isempty(defvalue), defvalue=default_value;end

if isstr(default_value);
 question_to_ask=strrep(question,'%',defvalue);
 y=input(question_to_ask,'s');
 if isempty(y);y=defvalue;end
else
 if length(defvalue)>0, s=num2str(defvalue(1));else s='';end
 for i=2:length(defvalue);s=[s ' ' num2str(defvalue(i))];end
 question_to_ask=strrep(question,'%',s);
 q=input(question_to_ask,'s');
 if isempty(q);
  y=defvalue;
 else
  eval(['y=[' q '];']);
 end
end
