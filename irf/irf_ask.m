function y=irf_ask(question,variable_name_of_old_value, default_value, seconds_to_wait)
%IRF_ASK   Ask a question
%
% y=irf_ask(question,variable_name_of_old_value, default_value);
% example:
%   y=irf_ask('How many 1/0? [%]>','y',0)
%   y=irf_ask('Large? yes/no [%]>','y','yes')
%   y=irf_ask('How much? [%]>','y',10)

useWaitInput=0;
if nargin < 3, disp('ERROR using irf_ask');help irf_ask;return;end
if nargin >= 3, def=default_value;else, def=[];end
if nargin==4, useWaitInput = 1; end
if evalin('caller',['exist(''' variable_name_of_old_value ''')'])==1
	defvalue=evalin('caller', variable_name_of_old_value);
	if isnan(defvalue), defvalue=def;end
else
	defvalue=def;
end
if isempty(defvalue), defvalue=default_value;end

ask_for_input=1; % flags whether input was ok (does not generate errors)
while ask_for_input
	ask_for_input=0; % if everything ok, do not ask for input anymore
	if ischar(default_value)
		question_to_ask=strrep(question,'%',defvalue);
		if useWaitInput
			y=waitinput(question_to_ask,seconds_to_wait,'s');
			if isnan(y) % use default value;
				y=defvalue;
			end
		else
			y=input(question_to_ask,'s');
		end
		if isempty(y);y=defvalue;end
	elseif isnumeric(default_value)
		if ~isempty(defvalue), s=num2str(defvalue(1));else, s='';end
		for i=2:length(defvalue);s=[s ' ' num2str(defvalue(i))];end %#ok<AGROW>
		question_to_ask=strrep(question,'%',s);
		if useWaitInput
			q=waitinput(question_to_ask,seconds_to_wait,'s');
			if isnan(q) % use default value;
				q=defvalue;
			end
		else
			q=input(question_to_ask,'s');
		end
		if isempty(q)
			y=defvalue;
		else
			eval(['y=[' q '];']);
		end
	end
end
