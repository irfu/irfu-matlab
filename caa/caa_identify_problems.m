function result = caa_identify_problems(data, probe, spacecraft_id, bitmask_column, quality_column)
%CAA_IDENTIFY_PROBLEMS  identifies problem areas in data, and sets appropriate bitmask and
%                       quality flag for these areas.
%
%  Input:
%     data              the set of data to go through.
%     probe             identifier for single probe or probe pair; i.e. {1,2,3,4,12,32,34}
%     spacecraft_id     the number of the Cluster spacecraft; i.e. {1,2,3,4}
%     bitmask_column    the number of the column in 'data' which holds the bitmask.
%     quality_column    the number of the column in 'data' which holds the quality flag.
%
% Output:
%     result            the set of data from input, with bitmask and quality flag added.
%
% $Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks

error(nargchk(3, 5, nargin))

result = data;
if isempty(data), warning('Data is empty: Cannot identify problems in empty data set.'), return, end
[rows columns] = size(data);

if nargin < 5
   quality_column = columns;           % Default quality flag to the last column in 'data'.
   if nargin < 4
      bitmask_column = columns - 1;    % Default bitmask to the second last column in 'data'.
   end
end

if ( (bitmask_column <= 0 && bitmask_column > columns) || (quality_column <= 0 && quality_column > columns) )
   error('Wrong column index(es) given.')
end

if ( spacecraft_id <= 0 || spacecraft_id > 4 )
   error('Wrong spacecraft ID given.')
end

if probe > 10
	switch probe
		case 12
			probe_list = [1, 2];
		case 32
			probe_list = [3, 2];
		case 34
			probe_list = [3, 4];
		otherwise
			error('Unknown probe.')
	end
elseif probe>0 && probe <=4
   probe_list = probe;
else
	error('Unknown probe.')
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
QUALITY_RESET                    =  1;
QUALITY_BAD_BIAS                 =  1;
QUALITY_PROBE_SATURATION         =  1;
QUALITY_LOW_DENSITY_SATURATION   =  1;
QUALITY_SWEEP_DATA               =  2;
QUALITY_BURST_DUMP               =  0;
QUALITY_NS_OPS                   =  1;
QUALITY_MANUAL_INTERVAL          =  1;
QUALITY_SINGLE_PROBE_PAIR        =  5;    % NOTE: Applies to L2 only.
QUALITY_ASYMMETRIC_MODE          =  8;
QUALITY_SOLAR_WIND_WAKE          =  8;
QUALITY_LOBE_WAKE                =  6;
QUALITY_PLASMASPHERE_WAKE        =  6;
QUALITY_WHISPER_OPERATING        =  7;    % NOTE: Applies to L2 only.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify problem areas in data, and set bitmask and quality flag for them.

% NOTE:  In case of multiple problems present in the same region, the quality
%        flag should be set to the lowest of the levels involved.
%        Since this is the case, check for problems in order of ascending
%        maximum quality!

% Remove bad data during burst dump:
[ok, problem_intervals, msg] = c_load('BDUMP?', spacecraft_id);
if ok
   if ~isempty(problem_intervals)
	   irf_log('proc', 'blanking bad data due to burst dump')
		%result = caa_rm_blankt(result, burst_dump_intervals);
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_BURST_DUMP, QUALITY_BURST_DUMP, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Remove bad bias around EFW reset
[ok, problem_intervals, msg] = c_load('BADBIASRESET?', spacecraft_id);
if ok
   if ~isempty(problem_intervals)
	   irf_log('proc', 'blanking bad bias due to EFW reset')
		%res = caa_rm_blankt(result,bbias);
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_RESET, QUALITY_RESET, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Remove bad bias from bias current indication
for probe_id = probe_list
	[ok, problem_intervals, msg] = c_load(irf_ssub('BADBIAS?p!', spacecraft_id, probe_id));
	if ok
		if ~isempty(problem_intervals)
			irf_log('proc', ['blanking bad bias on P' num2str(probe_id)])
			%res = caa_rm_blankt(res,bbias);
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_BAD_BIAS, QUALITY_BAD_BIAS, ...
               bitmask_column, quality_column);
		end
	else irf_log('load', msg)
	end
	clear ok problem_intervals msg
end


% Remove probe saturation
for probe_id = probe_list
	[ok, problem_intervals, msg] = c_load(irf_ssub('PROBESA?p!', spacecraft_id, probe_id));
	if ok
		if ~isempty(problem_intervals)
			irf_log('proc', ['blanking saturated P' num2str(probe_id)])
			%res = caa_rm_blankt(res,sa);
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_PROBE_SATURATION, QUALITY_PROBE_SATURATION, ...
               bitmask_column, quality_column);
		end
	else irf_log('load', msg)
	end
	clear ok problem_intervals msg
end
			

% Remove probe saturation due to low density
for probe_id = probe_list
	[ok, problem_intervals, msg] = c_load(irf_ssub('PROBELD?p!', spacecraft_id, probe_id));
	if ok
		if ~isempty(problem_intervals)
			irf_log('proc', ...
				['blanking low density saturation on P' num2str(probe_id)])
			%res = caa_rm_blankt(res,sa);
         result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_LOW_DENSITY_SATURATION, QUALITY_LOW_DENSITY_SATURATION, ...
               bitmask_column, quality_column);
		end
	else irf_log('load', msg)
	end
	clear ok problem_intervals msg
end


% Remove NS_OPS
%[ok, problem_intervals, msg] = c_load('NS_OPS?', spacecraft_id);
%if ok
%	if ~isempty(problem_intervals)
%		irf_log('proc', 'blanking NS_OPS')
%		%res = caa_rm_blankt(res,sweep);
%      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
%         BITMASK_NS_OPS, QUALITY_NS_OPS, ...
%            bitmask_column, quality_column);
%	end
%else irf_log('load', msg)
%end
%clear ok problem_intervals msg


% Remove data from manual intervals
%[ok, problem_intervals, msg] = c_load('MANUALINT?', spacecraft_id);
%if ok
%	if ~isempty(problem_intervals)
%		irf_log('proc', 'blanking data from manual intervals')
%		%res = caa_rm_blankt(res,sweep);
%      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
%         BITMASK_MANUAL_INTERVAL, QUALITY_MANUAL_INTERVAL, ...
%            bitmask_column, quality_column);
%	end
%else irf_log('load', msg)
%end
%clear ok problem_intervals msg


% Remove sweeps
[ok, problem_intervals, msg] = c_load('SWEEP?', spacecraft_id);
if ok
	if ~isempty(problem_intervals)
		irf_log('proc', 'blanking sweeps')
		%res = caa_rm_blankt(res,sweep);
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_SWEEP_DATA, QUALITY_SWEEP_DATA, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Remove single probe pair data
%[ok, problem_intervals, msg] = c_load('PROBEPAIR?', spacecraft_id);
%if ok
%	if ~isempty(problem_intervals)
%		irf_log('proc', 'blanking single probe pair data')
%		%res = caa_rm_blankt(res,whip);
%      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
%         BITMASK_SINGLE_PROBE_PAIR, QUALITY_SINGLE_PROBE_PAIR, ...
%            bitmask_column, quality_column);
%	end
%else irf_log('load', msg)
%end
%clear ok problem_intervals msg


% Remove plasmasphere wakes
[ok, problem_intervals, msg] = c_load(irf_ssub('PSWAKE?p!', spacecraft_id, probe));            
   if ok
       if ~isempty(problem_intervals)
           irf_log('proc', 'blanking plasmaspheric wakes')
           %res = caa_rm_blankt(res,wake,0,5);
           result = caa_set_bitmask_and_quality(result, problem_intervals, ...
            BITMASK_PLASMASPHERE_WAKE, QUALITY_PLASMASPHERE_WAKE, ...
               bitmask_column, quality_column);
       end
   else irf_log('load', msg)
   end
clear ok problem_intervals msg


% Remove lobe wakes
[ok, problem_intervals, msg] = c_load(irf_ssub('LOWAKE?p!', spacecraft_id, probe));
if ok
	if ~isempty(problem_intervals)
		irf_log('proc', 'blanking lobe wakes')
		%num_intervals = size(problem_intervals, 1);		
		%intervals = cell(1, num_intervals);
		%for ii=1:num_intervals
		%   intervals{ii} = fromepoch([problem_intervals(ii, 1); problem_intervals(ii,2)]);
		%end
		intervals = caa_parse_intervals_subfunc(problem_intervals);
		keyboard
		%res = caa_rm_blankt(res,wake,0,5);
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_LOBE_WAKE, QUALITY_LOBE_WAKE, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Remove whisper pulses
[ok, problem_intervals, msg] = c_load('WHIP?', spacecraft_id);
if ok
	if ~isempty(problem_intervals)
		irf_log('proc', 'blanking Whisper pulses')
		%res = caa_rm_blankt(res,whip);
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
         BITMASK_WHISPER_OPERATING, QUALITY_WHISPER_OPERATING, ...
            bitmask_column, quality_column);
	end
else irf_log('load', msg)
end
clear ok problem_intervals msg


% Remove data from asymmetric mode
%[ok, problem_intervals, msg] = c_load('ASYMM?', spacecraft_id);
%if ok
%	if ~isempty(problem_intervals)
%		irf_log('proc', 'blanking data from asymmetric mode')
%		%res = caa_rm_blankt(res,whip);
%      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
%         BITMASK_ASYMMETRIC_MODE, QUALITY_ASYMMETRIC_MODE, ...
%            bitmask_column, quality_column);
%	end
%else irf_log('load', msg)
%end
%clear ok problem_intervals msg


% Remove data from solar wind wake
%[ok, problem_intervals, msg] = c_load('SOLARWIND?', spacecraft_id);
%if ok
%	if ~isempty(problem_intervals)
%		irf_log('proc', 'blanking data from solar wind wake')
%		%res = caa_rm_blankt(res,whip);
%      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
%         BITMASK_SOLAR_WIND_WAKE, QUALITY_SOLAR_WIND_WAKE, ...
%            bitmask_column, quality_column);
%	end
%else irf_log('load', msg)
%end
%clear ok problem_intervals msg
			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intervals_out = caa_parse_intervals_subfunc(intervals_in)

   num_intervals = size(intervals_in, 1);
	intervals_out = cell(1, num_intervals);
	for ii=1:num_intervals
	   intervals_out{ii} = fromepoch([intervals_in(ii, 1); intervals_in(ii,2)]);
	end