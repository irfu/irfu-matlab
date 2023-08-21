function result = caa_identify_problems(result, data_level, probe, spacecraft_id, bitmask_column, quality_column, mask_type)
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
%     mask_type         0=default
%                       1=spacecraft potential P caa_export_cef() L2/3 processing
%                       2=internal burst PBurst L2
%                       3=internal burst EBurst L2
%                       4=internal burst BBurst L2 (BITMASK_RESET & BITMASK_PROBE_SATURATION only)
%
% Output:
%     result            the set of data from input, with bitmask and quality flag added.


%  Author:     Mikael Lundberg, Swedish Institute of Space Physics, <milu@irfu.se>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checks
narginchk(4, 7)
if nargin < 7, mask_type = 0; end

if size(result,1) <= 1
  warning('Short Data: Cannot identify problems in data set.'),
  return,
end
columns = size(result,2);
[data_start_time, data_time_span] = irf_stdt(result(1,1), result(end,1));

% For nonzero MASK_TYPE only some of the bitmask values
% and modified quality factors are used.
iburst=0;
sc_potential = 0; % Default 0. No sc potential processing (P).
switch mask_type
  case 0
  case 1, sc_potential = 1;              % spacecraft potential
  case 2, sc_potential = 1; iburst = 1;  % PBurst
  case 3, iburst = 1;                    % EBurst
  case 4, iburst = 1;                    % BBurst
  otherwise
end

qindex=sc_potential+1;                 % Index# in quality factors.
if nargin < 6
  quality_column = columns;           % Default quality flag to the last column in 'data'.
  if nargin < 5
    bitmask_column = columns - 1;    % Default bitmask to the second last column in 'data'.
  end
end
% Set quality to zero on P level 3 NaN values
if sc_potential && data_level == 3
  ix = isnan(result(:,2));
  result(ix,quality_column) = 0;
end

if ( (bitmask_column <= 0 || bitmask_column > columns) ||...
    (quality_column <= 0 || quality_column > columns) )
  error('Wrong column index(es) given.')
end

if ( data_level <= 0 || data_level > 3 )
  error('Incorrect level of data.')
end

if isempty(intersect(spacecraft_id,1:4))
  error('Wrong spacecraft ID given.')
end

if ~isa(probe, 'char') || isempty(probe)
  error('PROBE must be a string')
end

if ~regexp(probe, '^([1-4]|12|32|34|1234|3234)$')
  error('Wrong probe combination.')
end
probe_list = [];
probe_pair_list = [];
EFW_PROBE_PAIRS = [12,32,34,42];
if length(probe) > 1
  pTmp = probe;
  while ~isempty(pTmp)
    probe_list = [probe_list str2double(pTmp(1))]; pTmp(1) = []; %#ok<AGROW>
  end
  probe_list = sort(unique(probe_list));
  pTmp = probe;
  while ~isempty(pTmp)
    if length(pTmp)<2, break, end
    probe_pair_list = [probe_pair_list str2double(pTmp(1:2))]; pTmp(1:2) = []; %#ok<AGROW>
  end
  probe_pair_list = intersect(probe_pair_list,EFW_PROBE_PAIRS);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

% Bitmask values; 2^(bit_number - 1):
BITMASK_RESET                    =  1;       % Bit 1   PQuality 1 BBurstQuality 0
BITMASK_BAD_BIAS                 =  2;       % Bit 2   PQuality 1
BITMASK_PROBE_SATURATION         =  4;       % Bit 3   PQuality 1 BBurstQuality 0
BITMASK_LOW_DENSITY_SATURATION   =  8;       % Bit 4   PQuality 1
BITMASK_SWEEP_DATA               =  16;      % Bit 5
BITMASK_BURST_DUMP               =  32;      % Bit 6
BITMASK_NS_OPS                   =  64;      % Bit 7
BITMASK_MANUAL_INTERVAL          =  128;     % Bit 8
BITMASK_SINGLE_PROBE_PAIR        =  256;     % Bit 9
BITMASK_ASYMMETRIC_MODE          =  512;     % Bit 10
BITMASK_SOLAR_WIND_WAKE          =  1024;    % Bit 11
BITMASK_LOBE_WAKE                =  2048;    % Bit 12
BITMASK_PLASMASPHERE_WAKE        =  4096;    % Bit 13
BITMASK_WHISPER_OPERATING        =  8192;    % Bit 14  PQuality 0
BITMASK_HIGH_BIAS_SATURATION     =  16384;   % Bit 15  PQuality 2
BITMASK_BAD_DAC                  =  32768;   % Bit 16
BITMASK_PROBE_SHADOW             =  65536;   % Bit 17  PQuality 3

% Quality factors:
QUALITY_RESET                    =  [ 0 1 ]; % PQuality 1
QUALITY_BAD_BIAS                 =  [ 0 1 ]; % PQuality 1
QUALITY_PROBE_SATURATION         =  [ 0 1 ]; % PQuality 1
QUALITY_LOW_DENSITY_SATURATION   =  [ 0 1 ]; % PQuality 1
QUALITY_SWEEP_DATA               =  0;
QUALITY_BURST_DUMP               =  0;
QUALITY_NS_OPS                   =  0;
%QUALITY_MANUAL_INTERVAL          =  1;   % NOTE: Not used! Instead read from file and set on per-interval basis!
QUALITY_SINGLE_PROBE_PAIR        =  1;    % NOTE: Applies to L2 only.
QUALITY_ASYMMETRIC_MODE          =  2;    % NOTE: Applies to L2 only.
QUALITY_SOLAR_WIND_WAKE          =  3;
QUALITY_LOBE_WAKE                =  1;
QUALITY_PLASMASPHERE_WAKE        =  1;
QUALITY_WHISPER_OPERATING        =  [ 2 0 ]; % PQuality 0
QUALITY_HIGH_BIAS_SATURATION     =  [ 1 1 1 2]; % PQuality=2 (L3), =1 (L2), EQuality=1
QUALITY_BAD_DAC                  =  2;    % NOTE: Applies to L2 only.
QUALITY_PROBE_SHADOW             =  [ 1 3 2 ]; % PQuality 3, EQuality = 2 for L3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify problem areas in data, and set bitmask and quality flag for them.

% NOTE:  In case of multiple problems present in the same region, the quality
%        flag should be set to the lowest of the levels involved.
%        Since this is the case, check for problems in order of ascending
%        maximum quality!

if ~sc_potential && mask_type~=4
  % Mark bad data during burst dump:
  [ok, problem_intervals] = c_load('BDUMP?', spacecraft_id);
  if ok
    if ~isempty(problem_intervals)
      irf_log('proc', 'marking bad data due to burst dump')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
        BITMASK_BURST_DUMP, QUALITY_BURST_DUMP, ...
        bitmask_column, quality_column);
    end
  end
  clear ok problem_intervals msg
end

% Mark probe saturation due to internal burst spike filter
if iburst
  %    'iburst spike'
  [ok, problem_intervals] = c_load(irf_ssub('SPIKE?', spacecraft_id));
  if ok
    if ~isempty(problem_intervals)
      irf_log('proc', irf_ssub('marking saturated probe spike IB?',spacecraft_id))
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
        BITMASK_PROBE_SATURATION, QUALITY_PROBE_SATURATION(qindex), bitmask_column, quality_column);
    end
  end
  clear ok problem_intervals msg
end

% Mark bad bias around EFW reset
[ok, problem_intervals] = c_load('BADBIASRESET?', spacecraft_id);
if ok
  if ~isempty(problem_intervals)
    irf_log('proc', 'marking bad bias due to EFW reset')
    result = caa_set_bitmask_and_quality(result, problem_intervals, ...
      BITMASK_RESET, QUALITY_RESET(qindex), bitmask_column, quality_column);
  end
end
clear ok problem_intervals msg

if mask_type~=4
  % Mark bad bias from bias current indication
  for probe_id = probe_list
    if isnan(probe_id) % commissioning NaN probes handling
      return
    end
    [ok, problem_intervals] = c_load(irf_ssub('BADBIAS?p!', spacecraft_id, probe_id));
    if ok
      if ~isempty(problem_intervals)
        irf_log('proc', ['marking bad bias on P' num2str(probe_id)])
        result = caa_set_bitmask_and_quality(result, problem_intervals, ...
          BITMASK_BAD_BIAS, QUALITY_BAD_BIAS(qindex), bitmask_column, quality_column);
      end
    end
    clear ok problem_intervals msg
  end

  % Mark bad bias and high bias saturation from NS_OPS
  ns_ops = c_ctl('get', spacecraft_id, 'ns_ops');
  if isempty(ns_ops)
    c_ctl('load_ns_ops', [c_ctl('get', 5, 'data_path') '/caa-control'])
    ns_ops = c_ctl('get', spacecraft_id, 'ns_ops');
  end
  if ~isempty(ns_ops)
    ns_ops_intervals = [caa_get_ns_ops_int(data_start_time, data_time_span, ns_ops, 'bad_bias')' ...
      caa_get_ns_ops_int(data_start_time, data_time_span, ns_ops, 'spec_bias')']';
    if ~isempty(ns_ops_intervals)
      irf_log('proc', 'marking bad bias from NS_OPS')
      result = caa_set_bitmask_and_quality(result, ns_ops_intervals, ...
        BITMASK_BAD_BIAS, QUALITY_BAD_BIAS(qindex), bitmask_column, quality_column);
    end
    ns_ops_intervals = caa_get_ns_ops_int(data_start_time, data_time_span, ns_ops, 'high_bias');
    if ~isempty(ns_ops_intervals)
      irf_log('proc', 'marking high bias saturation from NS_OPS')
      result = caa_set_bitmask_and_quality(result, ns_ops_intervals, ...
        BITMASK_HIGH_BIAS_SATURATION, ...
        QUALITY_HIGH_BIAS_SATURATION((data_level-2)*2 + qindex), ...
        bitmask_column, quality_column);
    end
    clear ns_ops ns_ops_intervals
  end
end

% Mark probe saturation
for probe_id = probe_list
  [ok, problem_intervals] = c_load(irf_ssub('PROBESA?p!', spacecraft_id, probe_id));
  if ok
    if ~isempty(problem_intervals)
      irf_log('proc', ['marking saturated P' num2str(probe_id)])
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
        BITMASK_PROBE_SATURATION, QUALITY_PROBE_SATURATION(qindex), bitmask_column, quality_column);
    end
  end
  clear ok problem_intervals msg
end

if mask_type==4, return, end

% Mark probe saturation due to low density
for probe_id = probe_list
  [ok, problem_intervals] = c_load(irf_ssub('PROBELD?p!', spacecraft_id, probe_id));
  if ok
    if ~isempty(problem_intervals)
      irf_log('proc', ...
        ['marking low density saturation on P' num2str(probe_id)])
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
        BITMASK_LOW_DENSITY_SATURATION, QUALITY_LOW_DENSITY_SATURATION(qindex), bitmask_column, quality_column);
    end
  end
  clear ok problem_intervals msg
end

if ~sc_potential
  % Mark NS_OPS
  ns_ops = c_ctl('get', spacecraft_id, 'ns_ops');
  if isempty(ns_ops)
    c_ctl('load_ns_ops', [c_ctl('get', 5, 'data_path') '/caa-control'])
    ns_ops = c_ctl('get', spacecraft_id, 'ns_ops');
  end
  if ~isempty(ns_ops)
    ns_ops_intervals = [caa_get_ns_ops_int(data_start_time, data_time_span, ns_ops, 'bad_data')' ...
      caa_get_ns_ops_int(data_start_time, data_time_span, ns_ops, 'bad_tm')']';

    if ~isempty(ns_ops_intervals), irf_log('proc', 'marking NS_OPS'), end
    for k = 1:size(ns_ops_intervals, 1)
      result = caa_fill_ns_ops(result, ns_ops_intervals(k, :));   % Recreate time interval and fill this data with NaN.
      result = caa_set_bitmask_and_quality(result, ns_ops_intervals(k, :), ...
        BITMASK_NS_OPS, QUALITY_NS_OPS, bitmask_column, quality_column);
    end
    clear ns_ops ns_ops_intervals
  end

  % Mark data from manual intervals
  man_int = c_ctl('get', spacecraft_id, 'man_int');
  if isempty(man_int)
    c_ctl('load_man_int', [c_ctl('get', 5, 'data_path') '/caa-control'])
    man_int = c_ctl('get', spacecraft_id, 'man_int');
  end
  if ~isempty(man_int)
    [manual_intervals, qualities] = caa_get_manual_int(data_start_time, data_time_span, man_int);

    if ~isempty(manual_intervals), irf_log('proc', 'marking manual intervals'), end
    for k = 1:size(manual_intervals, 1)
      result = caa_set_bitmask_and_quality(result, manual_intervals(k, 1:2), ...
        BITMASK_MANUAL_INTERVAL, qualities(k, data_level-1), bitmask_column, quality_column);
    end
    clear manual_intervals qualities
  end
  clear man_int

  % Mark sweeps
  [ok, problem_intervals] = c_load('SWEEP?', spacecraft_id);
  if ok
    if ~isempty(problem_intervals)
      irf_log('proc', 'marking sweeps')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
        BITMASK_SWEEP_DATA, QUALITY_SWEEP_DATA, ...
        bitmask_column, quality_column);
    end
  end
  clear ok problem_intervals msg

  % Mark single probe pair data
  if (data_level == 2 && ~isempty(regexp(probe, '^(12|32|34)$', 'once')))
    irf_log('proc', 'marking single probe pair data')
    result(:, bitmask_column) = bitor(result(:, bitmask_column), BITMASK_SINGLE_PROBE_PAIR);
    result(:, quality_column) = min(result(:, quality_column), QUALITY_SINGLE_PROBE_PAIR);
  end

  % Mark plasmasphere wakes
  for probe_id = probe_pair_list
    [ok, problem_intervals] = c_load(irf_ssub('PSWAKE?p!', spacecraft_id, probe_id));
    if ok
      if ~isempty(problem_intervals)
        irf_log('proc', 'marking plasmaspheric wakes')
        result = caa_set_bitmask_and_quality(result, problem_intervals, ...
          BITMASK_PLASMASPHERE_WAKE, QUALITY_PLASMASPHERE_WAKE, ...
          bitmask_column, quality_column);
      end
    end
  end
  clear ok problem_intervals msg

  % Mark lobe wakes
  for probe_id = probe_pair_list
    [ok, problem_intervals] = c_load(irf_ssub('LOWAKE?p!', spacecraft_id, probe_id));
    if ok
      if ~isempty(problem_intervals)
        irf_log('proc', 'marking lobe wakes')
        result = caa_set_bitmask_and_quality(result, problem_intervals, ...
          BITMASK_LOBE_WAKE, QUALITY_LOBE_WAKE, ...
          bitmask_column, quality_column);
      end
    end
  end
  clear ok problem_intervals msg

  % Mark nonsinusoidal wakes
  for probe_id = probe_pair_list
    [ok, problem_intervals] = c_load(irf_ssub('NONSINWAKE?p!', spacecraft_id, probe_id));
    if ok
      if ~isempty(problem_intervals)
        irf_log('proc', 'marking nonsinusoidal wakes')
        result = caa_set_bitmask_and_quality(result, problem_intervals, ...
          BITMASK_LOBE_WAKE, QUALITY_LOBE_WAKE, ...
          bitmask_column, quality_column);
      end
    end
  end

  clear ok problem_intervals msg
end

% Mark whisper pulses
[ok, problem_intervals] = c_load('WHIP?', spacecraft_id);
if ok
  if ~isempty(problem_intervals)
    irf_log('proc', 'marking Whisper pulses')
    result = caa_set_bitmask_and_quality(result, problem_intervals, ...
      BITMASK_WHISPER_OPERATING, QUALITY_WHISPER_OPERATING(qindex), bitmask_column, quality_column);
  end
end
clear ok problem_intervals msg

if ~sc_potential
  % Mark saturation due to probe shadow for SAA=90 deg
  [ok, problem_intervals] = c_load('SAASADI?', spacecraft_id);
  if ok && ~isempty(problem_intervals)
    if data_level==3
      qTmp = QUALITY_PROBE_SHADOW(3);
      % Extend intervals for 2 sec at each side to catch the sfit time stamp
      problem_intervals(:,[1 3 5 7]) = problem_intervals(:,[1 3 5 7]) - 2;
      problem_intervals(:,[2 4 6 8]) = problem_intervals(:,[2 4 6 8]) + 2;
    else, qTmp = QUALITY_PROBE_SHADOW(1);
    end
    irf_log('proc', 'marking saturation due to probe shadow')
    allPairs = [12 34 32 42];
    for probe_id = probe_pair_list
      iP = find(probe_id==allPairs);
      result = caa_set_bitmask_and_quality(result, problem_intervals(:,iP*2-[1 0]), ...
        BITMASK_PROBE_SHADOW, qTmp, bitmask_column, quality_column);
    end
  end
  clear ok problem_intervals msg qTmp

  % Mark data from asymmetric mode
  if ( data_level == 2 && strcmp(probe, '3234') )
    % Asymmetric mode, p32 and p34 present
    irf_log('proc', 'marking data from asymmetric mode')
    result(:, bitmask_column) = bitor(result(:, bitmask_column), BITMASK_ASYMMETRIC_MODE);
    result(:, quality_column) = min(result(:, quality_column), QUALITY_ASYMMETRIC_MODE);
  end

  % Mark data with suspected bad DAC settings
  if data_level == 2
    for probe_loop={'12' '32' '34'}
      [ok, problem_intervals] = c_load(irf_ssub('BADDAC?p!', spacecraft_id,probe_loop{1}));
      if ok && ~isempty(problem_intervals)
        irf_log('proc', 'marking bad DAC interval')
        result = caa_set_bitmask_and_quality(result, problem_intervals, ...
          BITMASK_BAD_DAC, QUALITY_BAD_DAC,bitmask_column, quality_column);
      end
    end
    clear ok problem_intervals
  end

  % Mark data from solar wind wake
  for probe_id = probe_pair_list
    [ok, wake_info] = c_load(irf_ssub('WAKE?p!', spacecraft_id, probe_id));
    if ok
      if ~isempty(wake_info)
        irf_log('proc', 'marking data from solar wind wake')
        num_wakes = size(wake_info, 1);
        problem_intervals = zeros(num_wakes, 2);
        for k = 1:num_wakes
          center_time = wake_info(k, 1);      % First column gives position, in time, of wake center.
          problem_intervals(k, :) = center_time + [-1 1];  % Interval is 1 sec on either side of center.
        end
        result = caa_set_bitmask_and_quality(result, problem_intervals, ...
          BITMASK_SOLAR_WIND_WAKE, QUALITY_SOLAR_WIND_WAKE, ...
          bitmask_column, quality_column);
      end
    end
  end
  clear ok problem_intervals msg wake_info num_wakes center_time
end % ~sc_potential

% Mark high bias saturation
for probe_id = probe_pair_list
  [ok, problem_intervals] = c_load(irf_ssub('HBIASSA?p!', spacecraft_id, probe_id));
  if ok
    if ~isempty(problem_intervals)
      irf_log('proc', 'marking high bias saturations')
      result = caa_set_bitmask_and_quality(result, problem_intervals, ...
        BITMASK_HIGH_BIAS_SATURATION, ...
        QUALITY_HIGH_BIAS_SATURATION((data_level-2)*2 + qindex),...
        bitmask_column, quality_column);
    end
  end
end
clear ok problem_intervals msg

end % FUNCTION