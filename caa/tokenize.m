function tokens = tokenize(str, delimiter)
%TOKENIZE tokenize the string
%   TOKENIZE(S) returns a cellarray of tokens in the string
%   S delimited by "white space" characters.
%
%   TOKENIZE(S,D) returns a cellarray of tokens in the string
%   S delimited by one of the characters in D.
%
%   $Revision$  $Date$
%
%   See also STRTOK
%
% $Id$
%

% Copyright 2003 Yuri Khotyaintsev (yuri@irfu.se)

if nargin < 2
	delimiter = ' ';
end

i = 0;
tstr = str;

while length(tstr) > 1
	[token,reminder]=strtok(tstr, delimiter);
	i = i + 1;
	tokens{i} = token;
	if isempty(reminder)
		% disp('break')
		break;
	else
		% disp(sprintf('%d: %s',i,token))
		tstr = reminder(2:end);
	end
end
