% script av_ask
%
% if variable exist changes in question % to variable value
%                             otherwise % to default
% variable='a';default='';question='';av_ask

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_ask')

if exist(variable) & (isstr(default) == isstr(eval(variable)));
  qqwwqq=eval(variable);
else,
 qqwwqq=default;
end;
eval([variable '=irf_ask(question,''' variable ''',qqwwqq);']);
