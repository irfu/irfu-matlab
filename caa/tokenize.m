function tokens = tokenize(str, delimiter)
%TOKENIZE tokenize the string
%   TOKENIZE(S) returns a cellarray of tokens in the string
%   S delimited by "white space" characters.
%
%   TOKENIZE(S,D) returns a cellarray of tokens in the string
%   S delimited by string D.
%
%   $Revision$  $Date$
%
%   See also STRTOK,STRFIND
%
% $Id$
%

% Copyright 2003 Yuri Khotyaintsev (yuri@irfu.se)

if nargin < 2
	delimiter = ' ';
end

ind=strfind(str,delimiter); % find all start indexes of delimeter
s_ind=[1 ind+length(delimiter)]; s_ind(end)=[]; % start indexes
e_ind=ind-1;                                    % end indexes

i_string=1;
for j=1:length(s_ind),
  str_val=str(s_ind(j):e_ind(j));
  if str_val,
    tokens{i_string}=str_val;
    i_string=i_string+1;
  end
end

