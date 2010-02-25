function index=caa_locate_bisect(vec,value,indx2)
% Given the 2D array vec, returns the index such that value is between
% vec(index,indx2) and vec(index+1,indx2).  vec(:,indx2) must be monotonic.
% Uses bisection method based on Numerical Recipes' implementation.
% If the value is not in the range of vec, returns either 1 (below range)
% or length(vec) (above range).
%
% Input:
%     vec     vector to be searched
%     value   value to be searched for
%     indx2   which index in the 2D array to search
%  Output:
%     index   index such that value is between vec(index,indx2) and vec(index+1,indx2).
% 
%  Author:     Chris Cully, Swedish Institute of Space Physics, <chris@irfu.se>
%
 
  if value <= vec(1,indx2) 
       index=1;
       return;
  end
  
  [upper,ncols]=size(vec);
  if value >= vec(upper,indx2) 
       index=upper;
       return;
   end
     
   dir=vec(2,indx2) > vec(1,indx2);
   lower=0;
   while (upper-lower)>1
       index=floor((upper+lower)/2);
       if (value > vec(index,indx2)) == dir
           lower=index;
       else
           upper=index;
       end
   end
   index=lower;
   if index==0, index=1; end
end