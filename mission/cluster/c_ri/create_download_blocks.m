function d_b = create_download_blocks(start_time, Dt, block_size)
%
%Input:
%
%Output:.
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
start_time = toepoch(start_time);
end_time = start_time + Dt;
if Dt <= block_size
  d_b = [start_time Dt];

else

  time_l = start_time:block_size:end_time;
  time_l = time_l';
  l =length(time_l);

  if time_l(l) == end_time
    time_l = time_l(1:l-1);
  end

  [r,c] = size(time_l);
  dt_l(1:r) = block_size;
  dt_l = dt_l';

  d_b =[time_l dt_l];

  if time_l(r) < end_time
    d_b(r,2) = end_time - d_b(r,1);
  end

end
