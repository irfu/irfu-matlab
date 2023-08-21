function result = caa_set_bitmask_and_quality(result, time_int, bitmask, quality, bitmask_col, quality_col)
%CAA_SET_BITMASK_AND_QUALITY  sets the bitmask and quality factor for given intervals.
%
% Input:
%     result         data in Cluster AV format.
%     time_int       time intervals (ISDAT epoch).
%     bitmask        the bitmask to set for the given interval.
%     quality        the quality factor to set for the given interval.
%     bitmask_col    column of data which contains the bitmask [optional, default=second last]
%     quality_col    column of data which contains the quality factor   [optional, default=last]
%
%  Output:
%     result         the input data with bitmask and quality factor added for specified
%                    time intervals.
%
% Example:
%   E_WAKE = caa_set_bitmask_and_quality(diEs1p34, LOWAKE1p34, 2048, 6);
%
%  Author:     Mikael Lundberg, Swedish Institute of Space Physics, <milu@irfu.se>
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks.

narginchk(4,6)

if isempty(result) || isempty(time_int), return, end

time_int = sortrows(time_int, 1);      % Sort time intervals ascending according to column 1

if nargin < 6
  quality_col = size(result, 2);        % Default is last column (if not given).
  if nargin < 5
    bitmask_col = quality_col - 1;   % Default is second last column (if not given).
  end
end

fs = c_efw_fsample(result);
if fs <= 0
  irf_log('proc', 'Cannot set bitmask and quality (wrong sample freq.).')
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set bitmask and quality factors for given intervals.

% Obtain time interval boundaries around each (averaged) data point;
% +/- half the time given from sampling frequency:
delta_time = 1 / (2 * fs);

% If time is not monotonic, default to Mikael's code.
if min(result(2:end,1)-result(1:end-1,1)) < 0
  irf_log('proc','***** WARNING: Non-monotonic times detected. *****')
  data_time_lower = result(:, 1) - delta_time;
  data_time_upper = result(:, 1) + delta_time;

  for num = 1:size(time_int, 1)    % Set for each time interval in turn.
    % Get the row indices for the interval in which to set bitmask and quality factor.
    % The first column of time_int holds 'start_time' and second column 'end_time', for
    % each interval.
    problem_start = time_int(num, 1);
    problem_stop = time_int(num, 2);
    % If any part of a problem interval coincides with a data interval, then that data
    % interval is affected by the problem.
    % We first check if the start time of the current problem interval is inside a data interval.
    % Next, we check if the end time of the current problem interval is inside a data interval.
    % Lastly, we check the case where the current problem interval completely spans over a data
    % interval; i.e. width(current_problem_interval) >= width(data_interval).
    % In all of the above three cases, a matching data interval is affected and should be marked.
    row_index = (problem_start >= data_time_lower & problem_start < data_time_upper) | ...
      (problem_stop  >  data_time_lower & problem_stop <= data_time_upper) | ...
      (problem_start <= data_time_lower & problem_stop >= data_time_upper);


    % Mask the bits together, keeping any previously set bits.
    result(row_index, bitmask_col) = bitor(result(row_index, bitmask_col), bitmask);

    % Set the quality factor; given quality is maximum possible for region, so
    % set to the minimum of current and any previous.
    result(row_index, quality_col) = min(result(row_index, quality_col), quality);

  end

  return
end

% Time is monotonic, so use the bisection technique
for num = 1:size(time_int, 1)    % Set for each time interval in turn.
  % Trap NaNs in the time intervals
  if ~isfinite(time_int(num, 1)) || ~isfinite(time_int(num, 2))
    irf_log('proc','***** WARNING: NaNs in the problem interval log. *****')
    continue
  end

  % Get the row indices for the interval in which to set bitmask and quality factor.
  % The first column of time_int holds 'start_time' and second column 'end_time', for
  % each interval.
  lower=caa_locate_bisect(result,time_int(num, 1)+delta_time-1e-6,1);
  upper=caa_locate_bisect(result,time_int(num, 2)+delta_time-1e-6,1);
  if time_int(num, 1) > (result(lower,1)+delta_time) , lower=lower+1; end
  if upper < length(result(:,1))
    if time_int(num, 2) >= (result(upper+1,1)-delta_time) , upper=upper+1; end
  end

  if (upper-lower)==0
    if time_int(num, 1) > (result(lower,1)+delta_time), continue, end
    if time_int(num, 2) < (result(lower,1)-delta_time), continue, end
  end

  % Mask the bits together, keeping any previously set bits.
  result(lower:upper, bitmask_col) = bitor(result(lower:upper, bitmask_col), bitmask);

  % Set the quality factor; given quality is maximum possible for region, so
  % set to the minimum of current and any previous.
  result(lower:upper, quality_col) = min(result(lower:upper, quality_col), quality);


end
