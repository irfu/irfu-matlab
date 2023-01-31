function result = caa_fill_ns_ops(data, ns_ops, bitmask_col, fillval)
%CAA_FILL_GAPS(data,te)  fill gaps in the of a dataset
%
% res = caa_fill_gaps(data,te)
%       append NaNs to the end of DATA untill TE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks.

narginchk(2, 3)

if size(data,1)<2
  irf_log('proc','cannot fill gaps (not enough samples)')
  return
end

fs = c_efw_fsample(data);
if fs <= 0
  irf_log('proc','cannot fill gaps (invalid sample frequency)')
  return
end

if nargin < 4
  fillval = NaN;
  if nargin < 3
    bitmask_col = size(data, 2) - 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_before_ns_ops = data(data(:, 1) < ns_ops(1), :);
data_after_ns_ops = data(data(:, 1) > ns_ops(2), :);

if ~isempty(data_before_ns_ops)
  fill_interval = data_before_ns_ops(end, 1):(1/fs):ns_ops(2);
else     % NS OPS active from start of data set
  fill_interval = data(1,1):(1/fs):ns_ops(2);
end
fill_interval = fill_interval(fill_interval >= ns_ops(1) & fill_interval <= ns_ops(2));
fill_interval = fill_interval(:);

fill_data = zeros(length(fill_interval), size(data, 2));
fill_data(:, 1) = fill_interval;
fill_data(:, 2:end) = fillval;
fill_data(:, bitmask_col) = 0;

result = [data_before_ns_ops; fill_data; data_after_ns_ops];

% Do not report 1 point gaps
if length(fill_interval) > 1
  irf_log('proc',...
    sprintf('filling %d gaps at %s', length(fill_interval), epoch2iso(ns_ops(1), 1)))
end
