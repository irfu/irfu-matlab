%
% Union for logicals
%
% Calculates the logical union of A and B.
% See also union
function R=cefUnion(A,B)



if(islogical(A) && islogical(B))

    R=(A+B)>0;
    
else
  error('Nonlogical union not implemented yet please see matlab union')   
    
end

