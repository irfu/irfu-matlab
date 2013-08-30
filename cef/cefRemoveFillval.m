% Replace fillvalues  
%
% I = CEFREMOVEFILLVAL(INP, FILLVALUE, NEWVALUE)
% 
% See also cefRemoveNan

% By josef hook jhook@rssd.esa.int

function u=cefRemoveFillval(i, fillval,newval)
   
   if nargin == 2
    newval = NaN;
   end  
    
valid = find(not(i==fillval));
  
u = repmat(newval, size(i));
u(valid) = i(valid);