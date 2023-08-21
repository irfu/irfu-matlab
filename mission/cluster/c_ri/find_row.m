function [start_row, end_row] = find_row(start_time, end_time, Btdt, startposition)
%
% [start_row, end_row] = find_row(start_time, end_time, Btdt, startposition)
%
%Input:
% start_time -in epochs
% end_time -in epochs
% Btdt - matrix with [time, dt] which are the start of an event
%         and the duration of that event.
% startposition -where to start to search
%
%Output:
% start_row -which row that corresponds to the start_time
% end_row -which row that corresponds to the end_time
%
%Descrition of the function:
% Finds the rows in the matrix that corresponds to start_time and end_time
%
%Using:
%
%Work method:
%
%Error:
% -1 is returned if there are no valid rows.
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

[B_length, col] = size( Btdt);
start_row = -1;
end_row = -1;

if startposition > B_length
  startposition = B_length;
end

for i = startposition:B_length
  B_start = Btdt(i,1);
  B_dt = Btdt(i,2);
  B_end = B_start + B_dt;

  if i ~= B_length
    B_n = Btdt(i+1,1);
  else
    B_end = end_time +1;
  end


  %finds the start row
  if (start_time <= B_start && start_row == -1 && B_start < end_time) || (start_time >= B_start && start_time < B_end)
    start_row = i;
  end

  %finds the end row
  if (end_time <= B_end && end_row == -1 && start_row ~= -1) || (B_end < end_time && end_time < B_n && start_row ~= -1)
    end_row = i;
    return
  end

end % for-loop

%if the end of the passing_MP is outside the possible data
if end_row == -1 && start_row ~= -1
  end_row = B_length;
end
