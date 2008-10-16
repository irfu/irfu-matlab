function result = caa_identify_problems(data, data_level, probe, spacecraft_id, bitmask_column, quality_column)
%CAA_IDENTIFY_PROBLEMS  identifies problem areas in data, and sets appropriate bitmask and
%                       quality flag for these areas.
%
%  Input:
%     data              the set of data to go through.
%     data_level        the level of input data, i.e. 2 or 3.
%     probe             identifier for single probe or probe pair; i.e. {1,2,3,4,12,32,34}
%     spacecraft_id     the number of the Cluster spacecraft; i.e. {1,2,3,4}
%     bitmask_column    the number of the column in 'data' which holds the bitmask.
%     quality_column    the number of the column in 'data' which holds the quality flag.
%
% Output:
%     result            the set of data from input, with bitmask and quality flag added.
%
%  Author:     Mikael Lundberg, Swedish Institute of Space Physics, <milu@irfu.se>
%
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks

error(nargchk(4, 6, nargin))

result = data;
if isempty(data), warning('Data is empty: Cannot identify problems in empty data set.'), return, end
[rows columns] = size(data);

if nargin < 5
   quality_column = columns;           % Default quality flag to the last column in 'data'.
   if nargin < 4
      bitmask_column = columns - 1;    % Default bitmask to the second last column in 'data'.
   end
end

if ( (bitmask_column <= 0 || bitmask_column > columns) || (quality_column <= 0 || quality_column > columns) )
   error('Wrong column index(es) given.')
end

if ( spacecraft_id <= 0 || spacecraft_id > 4 )
   error('Wrong spacecraft ID given.')
end

if isstr(probe) & (regexp(probe, '^([1-4]|12|32|34|1234|3234)$') ~= 1)
   error('Wrong probe combination.')
%elseif isnumeric(sensor) && (regexp(num2str(sensor), '^([1-4]|12|32|34|1234|3234)$') ~= 1)
%   error('Wrong probe combination.')
elseif ~( isa(probe, 'char') || isa(probe, 'numeric') )
   error('Wrong probe format.')
end

if ( (isstr(probe) && length(probe) > 1) || (isnumeric(probe) && probe > 10) )
	switch probe
		case {12, '12'}
			probe_list = [1, 2];
			probe_pair_list = 12;
		case {32, '32'}
			probe_list = [3, 2];
			probe_pair_list = 32;
		case {34, '34'}
			probe_list = [3, 4];
			probe_pair_list = 34;
		case {1234, '1234'}              % Do nothing?
		   probe_list = [1, 2, 3, 4];
		   probe_pair_list = [12, 34];
		case {3234, '3234'}              % Do nothing?
		   probe_list = [2, 3, 4];
		   probe_pair_list = [32, 34];
		otherwise
			error('Unknown probe.')
	end
%elseif isnumeric(probe) && (probe>0 && probe <=4)
%   probe_list = probe;
elseif isstr(probe) & (regexp(probe, '^[1-4]$') == 1)
   probe_list = str2num(probe);
   probe_pair_list = [];
else
	error('Unknown probe.')
end

if ( data_level <= 0 || data_level > 3 )
   error('Incorrect level of data.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

% Bitmask values; 2^(bit_number - 1):
BITMASK_RESET                    =  1;    % Bit 1
BITMASK_BAD_BIAS                 =  2;    % Bit 2
BITMASK_PROBE_SATURATION         =  4;    % Bit 3
BITMASK_LOW_DENSITY_SATURATION   =  8;    % Bit 4
BITMASK_SWEEP_DATA               =  16;   % Bit 5
BITMASK_BURST_DUMP               =  32;   % Bit 6
BITMASK_NS_OPS                   =  64;   % Bit 7
BITMASK_MANUAL_INTERVAL          =  128;  % Bit 8
BITMASK_SINGLE_PROBE_PAIR        =  256;  % Bit 9
BITMASK_ASYMMETRIC_MODE          =  512;  % Bit 10
BITMASK_SOLAR_WIND_WAKE          =  1024; % Bit 11
BITMASK_LOBE_WAKE                =  2048; % Bit 12
BITMASK_PLASMASPHERE_WAKE        =  4096; % Bit 13
BITMASK_WHISPER_OPERATING        =  8192; % Bit 14

% Quality factors:
QUALITY_RESET                    =  0;
QUALITY_BAD_BIAS                 =  0;
QUALITY_PROBE_SATURATION         =  0;
QUALITY_LOW_DENSITY_SATURATION   =  0;
QUALITY_SWEEP_DATA               =  0;
QUALITY_BURST_DUMP               =  0;
QUALITY_NS_OPS                   =  0;
QUALITY_MANUAL_INTERVAL          =  1;    % NOTE: Not used! Instead read from file and set on per-interval basis!
QUALITY_SINGLE_PROBE_PAIR        =  1;    % NOTE: Applies to L2 only.
QUALITY_ASYMMETRIC_MODE          =  2;    % NOTE: Applies to L2 only.
QUALITY_SOLAR_WIND_WAKE          =  3;
QUALITY_LOBE_WAKE                =  0;
QUALITY_PLASMASPHERE_WAKE        =  0;
QUALITY_WHISPER_OPERATING        =  3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify problem areas in data, and set bitmask and quality flag for them.

% NOTE:  In case of multiple problems present in the same region, the quality
%        flag should be set to the lowest of the levels involved.
%        Since this is the case, check for problems in order of ascending
%        maximum quality!

% Mark bad data during burst dump:
[ok, problem_intervals, msg] = c_load('BDUMP?', spacecraft_id);
if ok
   if ~isempty(problem_intervals)
	   irf_log('proc', 'marking bad data due to burst dump')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_BURST_DUMP, QUALITY_BURST_DUMP, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Mark bad bias around EFW reset
[ok, problem_intervals, msg] = c_load('BADBIASRESET?', spacecraft_id);
if ok
   if ~isempty(problem_intervals)
	   irf_log('proc', 'marking bad bias due to EFW reset')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_RESET, QUALITY_RESET, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Mark bad bias from bias current indication
for probe_id = probe_list
	[ok, problem_intervals, msg] = c_load(irf_ssub('BADBIAS?p!', spacecraft_id, probe_id));
	if ok
		if ~isempty(problem_intervals)
			irf_log('proc', ['marking bad bias on P' num2str(probe_id)])
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_BAD_BIAS, QUALITY_BAD_BIAS, ...
               bitmask_column, quality_column);
		end
	else irf_log('load', msg)
	end
	clear ok problem_intervals msg
end


% Mark probe saturation
for probe_id = probe_list
	[ok, problem_intervals, msg] = c_load(irf_ssub('PROBESA?p!', spacecraft_id, probe_id));
	if ok
		if ~isempty(problem_intervals)
			irf_log('proc', ['marking saturated P' num2str(probe_id)])
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_PROBE_SATURATION, QUALITY_PROBE_SATURATION, ...
               bitmask_column, quality_column);
		end
	else irf_log('load', msg)
	end
	clear ok problem_intervals msg
end
			

% Mark probe saturation due to low density
for probe_id = probe_list
	[ok, problem_intervals, msg] = c_load(irf_ssub('PROBELD?p!', spacecraft_id, probe_id));
	if ok
		if ~isempty(problem_intervals)
			irf_log('proc', ...
				['marking low density saturation on P' num2str(probe_id)])
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_LOW_DENSITY_SATURATION, QUALITY_LOW_DENSITY_SATURATION, ...
               bitmask_column, quality_column);
		end
	else irf_log('load', msg)
	end
	clear ok problem_intervals msg
end


% Mark NS_OPS
ns_ops = c_ctl('get', spacecraft_id, 'ns_ops');
if isempty(ns_ops)
   c_ctl('load_ns_ops', [c_ctl('get', 5, 'data_path') '/caa-control'])
	ns_ops = c_ctl('get', spacecraft_id, 'ns_ops');
end
if ~isempty(ns_ops)
   data_start_time = result(1, 1);
   data_time_span = result(end, 1) - data_start_time;
   ns_ops_intervals = caa_get_ns_ops_int(data_start_time, data_time_span, ns_ops, 'bad_data');
   
   if ~isempty(ns_ops_intervals), irf_log('proc', 'marking NS_OPS'), end
   for k = 1:size(ns_ops_intervals, 1)
      result = caa_fill_ns_ops(result, ns_ops_intervals(k, :));   % Recreate time interval and fill this data with NaN.
      result = caa_set_bitmask_and_quality(result, ns_ops_intervals(k, :), ...
         BITMASK_NS_OPS, QUALITY_NS_OPS, bitmask_column, quality_column);
   end
   clear ns_ops data_start_time data_time_span ns_ops_intervals
end


% Mark data from manual intervals
man_int = c_ctl('get', spacecraft_id, 'man_int');
if isempty(man_int)
   c_ctl('load_man_int', [c_ctl('get', 5, 'data_path') '/caa-control'])
	man_int = c_ctl('get', spacecraft_id, 'man_int');
end
if ~isempty(man_int)
   data_start_time = result(1, 1);
   data_end_time = result(end, 1);
   data_time_span = data_end_time - data_start_time;
   
   [manual_intervals, qualities] = caa_get_manual_int(data_start_time, data_time_span, man_int);
   
   if ~isempty(manual_intervals), irf_log('proc', 'marking manual intervals'), end
   for k = 1:size(manual_intervals, 1)
      result = caa_set_bitmask_and_quality(result, manual_intervals(k, 1:2), ...
         BITMASK_MANUAL_INTERVAL, qualities(k, data_level-1), bitmask_column, quality_column);
   end
   clear data_start_time data_time_span data_end_time manual_intervals qualities
end
clear man_int
      

% Mark sweeps
[ok, problem_intervals, msg] = c_load('SWEEP?', spacecraft_id);
if ok
	if ~isempty(problem_intervals)
		irf_log('proc', 'marking sweeps')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_SWEEP_DATA, QUALITY_SWEEP_DATA, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Mark single probe pair data
if (data_level == 2 & regexp(probe, '^(12|32|34)$'))
   irf_log('proc', 'marking single probe pair data')
   result(:, bitmask_column) = bitor(result(:, bitmask_column), BITMASK_SINGLE_PROBE_PAIR);
   result(:, quality_column) = min(result(:, quality_column), QUALITY_SINGLE_PROBE_PAIR);
end


% Mark plasmasphere wakes
for probe_id = probe_pair_list
   [ok, problem_intervals, msg] = c_load(irf_ssub('PSWAKE?p!', spacecraft_id, probe_id));            
   if ok
       if ~isempty(problem_intervals)
           irf_log('proc', 'marking plasmaspheric wakes')
           result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_PLASMASPHERE_WAKE, QUALITY_PLASMASPHERE_WAKE, ...
               bitmask_column, quality_column);
       end
   else irf_log('load', msg)
   end
end
clear ok problem_intervals msg


% Mark lobe wakes
for probe_id = probe_pair_list
   [ok, problem_intervals, msg] = c_load(irf_ssub('LOWAKE?p!', spacecraft_id, probe_id));
   if ok
   	if ~isempty(problem_intervals)
   		irf_log('proc', 'marking lobe wakes')
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_LOBE_WAKE, QUALITY_LOBE_WAKE, ...
               bitmask_column, quality_column);
   	end
   else irf_log('load', msg)
   end
end
clear ok problem_intervals msg


% Mark whisper pulses
[ok, problem_intervals, msg] = c_load('WHIP?', spacecraft_id);
if ok
	if ~isempty(problem_intervals)
		irf_log('proc', 'marking Whisper pulses')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_WHISPER_OPERATING, QUALITY_WHISPER_OPERATING, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Mark data from asymmetric mode
if ( data_level == 2 && strcmp(probe, '3234') )
   % Asymmetric mode, p32 and p34 present
   irf_log('proc', 'marking data from asymmetric mode')
   result(:, bitmask_column) = bitor(result(:, bitmask_column), BITMASK_ASYMMETRIC_MODE);
   result(:, quality_column) = min(result(:, quality_column), QUALITY_ASYMMETRIC_MODE);
end


% Mark data from solar wind wake
for probe_id = probe_pair_list
   [ok, wake_info, msg] = c_load(irf_ssub('WAKE?p!', spacecraft_id, probe_id));
   if ok
   	if ~isempty(wake_info)
   		irf_log('proc', 'marking data from solar wind wake')
   		num_wakes = size(wake_info, 1);
   		problem_intervals = zeros(num_wakes, 2);
   		for k = 1:num_wakes
   		   center_time = wake_info(k, 1);      % First column gives position, in time, of wake center.
   		   problem_intervals(k, :) = center_time + [-1 1];  % Interval is 1 sec on either side of center.
   		end
   		disp('CHECK wake_info (center times) and problem_intervals !!!'), keyboard
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_SOLAR_WIND_WAKE, QUALITY_SOLAR_WIND_WAKE, ...
               bitmask_column, quality_column);
   	end
   else irf_log('load', msg)
   end
end
clear ok problem_intervals msg wake_info num_wakes center_time
			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intervals_out = caa_parse_intervals_subfunc(intervals_in)

   num_intervals = size(intervals_in, 1);
	intervals_out = cell(1, num_intervals);
	for ii=1:num_intervals
	   intervals_out{ii} = fromepoch([intervals_in(ii, 1); intervals_in(ii,2)]);
	end