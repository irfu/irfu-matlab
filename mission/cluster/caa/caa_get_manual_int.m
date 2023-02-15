function [intervals, qualities] = caa_get_manual_int(st,dt,man_int)
%CAA_GET_MANUAL_INT  Get those manual intervals which intersect the data interval.
%                    Also retrieve matching L2 and L3 qualities for the intervals.
%
% ints_out = caa_get_manual_int(st, dt, man_int)
% [intervals, qualities] = caa_get_manual_int(st, dt, man_int)
%
% See also: C_CTL
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks.

narginchk(3, 3)

if isempty(man_int), error('Empty man_int'), end

intervals = [];
qualities = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find intersecting intervals and their qualities.

data_start_time = st;
data_end_time = st + dt;

problem_start_time = man_int(:, 1);
problem_end_time = problem_start_time + man_int(:, 2);

row_index = (problem_start_time >= data_start_time & problem_start_time < data_end_time) | ...
  (problem_end_time > data_start_time & problem_end_time <= data_end_time) | ...
  (problem_start_time <= data_start_time & problem_end_time >= data_end_time);

intervals = [man_int(row_index, 1),  man_int(row_index, 1)+man_int(row_index, 2)];
qualities = man_int(row_index, 3:4);

if nargout < 2
  intervals = [intervals, qualities];
end