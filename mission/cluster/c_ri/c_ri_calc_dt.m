function dt = c_ri_calc_dt(B,k,mode)
%
% dt = c_ri_calc_dt(B,k,mode)
%
%Input:
% B -timestaps in left vecter [timestamp | ? | ?], in seconds
% k -nr of points in row
% mode -burst or normal (b/n) b=67 Hz, n=22 Hz
%
%Output:
% dt -a calculation of the periodtime between samples
%
%Descrition of the function:
% Finds the period between samples, the number of points that the period i calculated
% over i k.
%
%Using:
%
%Work method:
%
%Error:
% sends -1 of dt not is found
%
%Description of variables:
% the number of steps to calculate a period
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
nr = 0;
i = 1;
[i_max, col] = size(B);

if i_max > k

  if mode == 'b'
    est_dt = 1/67;
  end

  if mode == 'n'
    est_dt = 1/22;
  end

  diff = est_dt/2;
  found_dt = 0;

  while nr < k && i < i_max && found_dt == 0

    B_old = B(i,1);
    i = i+1;
    B_new = B(i,1);

    if (B_new < B_old + est_dt +diff) && (B_new > B_old + est_dt - diff)

      if nr == 0
        start_nr = i-1;
      end
      nr = nr +1;

      if nr == k
        end_nr = i;
        found_dt = 1;
      end

    else
      nr = 0;
      end_nr = 0;
    end

  end


  if found_dt == 0
    dt = -1;
  else
    dt = (B(end_nr,1) - B(start_nr,1))/(k - 1);
  end

else
  dt = -1;
end
