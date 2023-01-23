function tokens = tokenize(str, delimiter)
%TOKENIZE  extract tokens from string
%
%   TOKENIZE(S) returns a cellarray of tokens in the string
%   S delimited by "white space" characters.
%
%   TOKENIZE(S,D) returns a cellarray of tokens in the string
%   S delimited by string D.
%
%   See also STRTOK,STRFIND
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2)

if nargin < 2
  delimiter = ' ';
end

ind = strfind(str,delimiter); % find all start indexes of delimeter
s_ind = [1 ind+length(delimiter)]; % start indexes
e_ind = [ind-1 length(str)];       % end indexes

tokens = cell(1,length(s_ind));
i_string = 1;
for j=1:length(s_ind)
  str_val = str(s_ind(j):e_ind(j));
  if str_val
    tokens{i_string} = str_val;
    i_string = i_string + 1;
  else
    tokens(i_string) = [];
  end
end

