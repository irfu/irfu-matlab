function y=irf_ask(question,variable_name_of_old_value, default_value);
%IRF_ASK   Ask a question
%
% y=irf_ask(question,variable_name_of_old_value, default_value);
% example:
%   y=irf_ask('How many 1/0? [%]>','y',0)
%   y=irf_ask('Large? yes/no [%]>','y','yes')
%   y=irf_ask('How much? [%]>','y',10)
%
% $Id$

if nargin < 3, disp('ERROR using irf_ask');help irf_ask;return;end
if nargin == 3, def=default_value;else def=[];end

if evalin('caller',['exist(''' variable_name_of_old_value ''')'])==1,
  defvalue=evalin('caller', [variable_name_of_old_value]);
else
  defvalue=def;
end
if isempty(defvalue), defvalue=default_value;end

ask_for_input=1; % flags whether input was ok (does not generate errors)
while ask_for_input
  ask_for_input=0; % if everything ok, do not ask for input anymore
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
end