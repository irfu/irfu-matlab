% script av_ask
%
% if variable exist changes in question % to variable value
%                             otherwise % to default
% variable='a';default='';question='';av_ask


if exist(variable) & (isstr(default) == isstr(eval(variable)));
  qqwwqq=eval(variable);
else,
 qqwwqq=default;
end;
eval([variable '=av_q(question,''' variable ''',qqwwqq);']);
