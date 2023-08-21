function B_c = c_ri_comple(time_line,B,s_per)
%
%Input:
%
%Output:
%
%Descrition of the function:
%
%Using:
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if s_per == -1

  B_c = 0;

else

  [i_end,c] = size(B);
  [t_end,c] = size(time_line);

  B_c(1:t_end,1) = time_line;
  B_c(1:t_end, 2:4) = 0;
  t_s = time_line(1);

  for i = 1:i_end
    t_line_pos = round((B(i,1)-t_s)/s_per);

    if t_line_pos > 0 && t_line_pos <= t_end
      B_c(t_line_pos,2:4) = B(i,2:4);
    end

  end
end
