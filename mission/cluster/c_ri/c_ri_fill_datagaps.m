function B_out = c_ri_fill_datagaps(B)
%
%Input:
% B -data,[time, Bx,By,By]
%
%Output:
% B_out -B-data without datagaps.
%
%Descrition of the function:
% Removes datagaps in a B-data matrix. A datagap is defined as 0.
% Removes only datagap of 1 sample lenght
%
%Using:
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
[b_end, c] = size(B);
B_out = B;

for b = 2:b_end-1
 for c = 2:4
  if B(b,c) == 0
  t1 = B(b-1,c);
  t2 = B(b+1,c);
  nr_no_zero = 0;
	if t1 ~=0
	nr_no_zero = nr_no_zero +1;
	end
	if t2~=0
	nr_no_zero = nr_no_zero +1;
	end
	if nr_no_zero ~= 0
	B_out(b,c) = (t1+t2)/nr_no_zero;
	end
  end
 end
end

