function [B_vector, position] = find_vector_to_time(B, n, t, dt)
%
%Input: One matrix containing the time and the B-components, The positon
% where the last sample was used. t is the time where you want to sample
% and ds is the time between samles.
%
%Output: The vector att the sampling time. Zero-vector if there is no
% valid vector.
%
%Descrition of the function:
% Finds the vector corresponing to the time in input.
%
%Using:
%
%Work method:
% Start the searching
% at position n. If the time at n is earlier than time. Then there is a search
% forward. I no vector i found a zero-vector is returned and n is unchanged.
%
%Error:
% If there is no match in the time frame, a zero-vector is returned and n is
% unchanged
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

%to avoid reading outside the matrix, if matrix size is bigger than n, it is ok for n+1
if max(size(B)) > n

  % in range of time window
  if (B(n+1,1) > t - dt / 2) && ( B(n+1,1) < t + dt /2)
    B_vector = B(n+1,2:4);
    position = n+1;

    %out of range in time window
  else
    % sets default values
    B_vector = [0 0 0];
    position= n;

    %searching forward in time
    k = 0;
    % to
    while B(n + k,1) < t - dt /2  &&  max(size(B)) >= n + k + 1
      k = k + 1;

      % if position found within time frame
      if (B(n + k,1) > t - dt /2) && ( B(n+1,1) < t + dt /2)
        B_vector = B(n+k,2:4);
        position = n+k;
      end

    end
  end

  %if the end of a matrix i reached then the zero-vector and the end-position is sent.
else
  B_vector = [0 0 0];
  position = max(size(B)) + 1;
end
