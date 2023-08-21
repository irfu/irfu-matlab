%
% Class to collect at least some short CDF time functions.
%
%
% CONVENTIONS
% ===========
% TT2000WOLS = ~TT2000 WithOut Leap Seconds.
%              Number of nanoseconds since roughly TT2000 epoch, excluding leap
%              seconds.
%              NOTE: Epoch is somewhat arbitrary and should not be relied
%              upon.
%              NOTE: Has no representation for time during positive leap
%              seconds. Has double representation for time during negative leap
%              seconds.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-10-21
%
classdef time
  % PROPOSAL: Should be package.
  %
  % PROPOSAL: Add more CDF time functions.
  %   Ex: TT2000_to_datevec, TT2000_to_UTC_str.
  %   CON: Longer call names (package+class+name).
  % PROPOSAL: Split into function files.
  %   CON: TT2000WOLS functions need constants.
  %
  % PROPOSAL: Wrapper for spdfparsett2000: UTC string(s) --> TT2000
  %   PRO: Can have easier to remember name: UTC_string_to_TT2000
  %   PRO: Can have assertion for failing to parse.
  %       Ex: spdfparsett2000({'2020-01-01'}) == (int64(-inf)+3)
  %   PRO: Can force input size to equal output size.
  % PROPOSAL: Wrapper for spdfcomputett2000: date vector --> TT2000
  %   PRO: Can have easier to remember name: datevec9_to_TT2000
  %   PRO: Can handle case of zero timestamps which spdfcomputett2000 can not.
  %   PRO: Can handle case of date vector on illegal class causing
  %        spdfcomputett2000 to give the WRONG return value!!

  % PROPOSAL: Wrapper for spdfbreakdowntt2000: TT2000 --> date vector
  %   PRO: Can have easier to remember name: TT2000_to_datevec9
  %   PRO: Can handle case of zero timestamps which spdfbreakdowntt2000 can not.
  %
  % PROPOSAL: Function for rounding to nearest non-leap second.


  properties(Constant, Access=private)

    % Approximately epoch of TT2000 as SDN. Must be integer.
    TT2000WOLS_EPOCH_SDN = datenum([2000,1,1]);

    % Range of days around epoch for which TT2000WOLS can be represented.
    % Must be less than corresponding range for TT2000.
    TT2000WOLS_DAY_RANGE = 290*365;

  end



  properties(Constant, Access=public)
    % Empirical: spdfparsett2000() returns this value for non-empty strings
    % when it can not interpret it as a UTC string.
    % PROPOSAL: Better name?
    % TODO-NI: Used by spdf code for other purposes too?
    %TT2000_CAN_NOT_INTERPRET_STR = int64(-9223372036854775805);
    TT2000_CAN_NOT_INTERPRET_STR = spdfparsett2000('A');
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % tt2000 : Column array. TT2000 values. Must not be during leap second.
    function tt2000wols = TT2000_to_TT2000WOLS(tt2000)
      v = int64(spdfbreakdowntt2000(tt2000));
      assert(v(6) < 60, 'Can not convert for positive leap second.')

      % IMPLEMENTATION NOTE: datenum only accepts double arrays.
      % NOTE: Not SDN.
      dayNbr = datenum(double(v(:,1:3))) - irf.cdf.time.TT2000WOLS_EPOCH_SDN;

      % int64 if v is int64.
      tt2000wols = irf.utils.mixed_radix.MRD_to_integer(...
        [fliplr(v(:, 4:9)), dayNbr], ...
        [1000, 1000, 1000, 60, 60, 24, irf.cdf.time.TT2000WOLS_DAY_RANGE]');

    end



    function tt2000 = TT2000WOLS_to_TT2000(tt2000wols)

      mrd = irf.utils.mixed_radix.integer_to_MRD(...
        tt2000wols, ...
        [1000, 1000, 1000, 60, 60, 24, irf.cdf.time.TT2000WOLS_DAY_RANGE]');

      % NOTE: Not SDN.
      dayNbr = mrd(:,7);
      if tt2000wols < 0
        % Correction needed because of how the day number is extracted
        % as if it were a digit in a mixed radix number.
        dayNbr = dayNbr - irf.cdf.time.TT2000WOLS_DAY_RANGE;
      end
      % SDN
      sdn = dayNbr + irf.cdf.time.TT2000WOLS_EPOCH_SDN;

      % IMPLEMENTATION NOTE: datevec only accepts double.
      dateVec = datevec(double(sdn));

      v = [dateVec(:, 1:3), fliplr(mrd(:, 1:6))];
      tt2000 = spdfcomputett2000(double(v));
    end



  end

end
