function M = add_A2M(M,A)
%
%Input:
% A -matrix
% M -matrix
%
%Output:
% M -Matrix M when A is added to M in the end of M.
%
%Descrition of the function:
% Adds matrix A to the end of matrix M
%
%Using:
%
%Work method:
%
%Error:
% If A is zero, nothing is added
% If M is zero, M becomes A.
% If A is zerom, M is unchanged
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if A ~= 0

  if M == 0
    M = A;

  else
    [new_lines,col] = size(A);
    [old_lines,col] = size(M);
    from = old_lines + 1;
    to = old_lines + new_lines;
    M(from:to,:) = A;
  end

end
