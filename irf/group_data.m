function grouped_data = group_data(data,threshold)
% GROUP_DATA Group sorted data when the difference between consecutive
% elements is less than a specified threshold. The original idea for this
% function was to use it with an array of TT2000 times to group them
% into an array of TT2000 tints.
%
% Input:
%     data - sorted 1D numeric array.
%     threshold - number to determine how array is split and grouped.
%
% Output:
%     grouped_data - 2D array of grouped data, see example below.
%
% Example:
%     A = [1,2,4,5,10,11,13,28,30,31,47];
%     thr = 3;
%     grouped_data = group_data(A,thr)
%     ans =
%          1     5
%         10    13
%         28    31
%         47    47

% If invalid arguments given, show help
if (nargin ~= 2) || ~isvector(data)
  help irf.group_data;
  grouped_data = [];
  return;
end

% convert to row vector if column vector given
if iscolumn(data)
  data = data';
end

% Find where difference between consecutive elements exceeds threshold
split_idx = [1, find(diff(data) >= threshold) + 1, length(data) + 1];

% Preallocate output array
num_groups = length(split_idx) - 1;
grouped_data = zeros(num_groups, 2);

% Populate output array with start and end of each group
for i = 1:num_groups
  grouped_data(i, 1) = data(split_idx(i));       % start of group
  grouped_data(i, 2) = data(split_idx(i+1) - 1); % end of group
end
end