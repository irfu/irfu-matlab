function col = torow(vector)
% TOROW - outputs a row vector whatever the format of
%   the input vector. Returns zero if input is rank higher 
%   than 1 (matrices etc), or input if input is empty or
%   scalar.
%
% See also TOCOLUMN

% Anders.Eriksson@irfu.se 021208

if isempty(vector)                     % input is empty
  col = vector;
else
  s = size(vector);
  if (max(size(s)) < 3) && find(s == 1) % input is a vector (or scalar)
    if s(1) > 1                        % input is column vector
      col = vector';
    else                               % input is row vector (or scalar)
      col = vector;
    end
  else                                 % input is matrix or higher rank
    col = 0;
  end
end  
