function result = caa_set_bitmask_and_quality(data, time_int, bitmask, quality, bitmask_col, quality_col)
%CAA_SET_BITMASK_AND_QUALITY  sets the bitmask and quality factor for given intervals.
%
% Input:
%     data           Cluster AV format.
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

% ----------------------------------------------------------------------------
% <milu@irfu.se> wrote this file.  It is a derivative work from Yuri's 
% file caa_rm_blankt.m, which carried the following text.   Mikael Lundberg
%
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks.

error(nargchk(4,6,nargin))

result = data;
if isempty(data) || isempty(time_int), return, end

time_int = sort(time_int, 1);    % Possible BUG! Columns are not kept together, but sorted individually!! (ML)

if nargin < 6
   quality_col = size(data, 2);        % Default is last column (if not given).
   if nargin < 5
      bitmask_col = quality_col - 1;   % Default is second last column (if not given).
   end
end

disp('set_bm_and_q'), keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set bitmask and quality factors for given intervals.

for num = 1:size(time_int, 1)    % Set for each time interval in turn.
   % Get the row indices for the interval in which to set bitmask and quality factor.
   % The first column of time_int holds 'start_time' and second column 'end_time', for
   % each interval
   row_index = data(:, 1) >= time_int(num, 1) & data(:, 1) <= time_int(num, 2);
   
   % Mask the bits together, keeping any previously set bits.
   result(row_index, bitmask_col) = bitor(result(row_index, bitmask_col), bitmask);
   
   % Set the quality factor; given quality is maximum possible for region, so
   % set to the minimum of current and any previous.
   result(row_index, quality_col) = min(result(row_index, quality_col), quality);
   
end

    
