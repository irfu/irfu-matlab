%
% Class which measures wall time between two points in code (typically within
% the same function) and logs the result in a standardized way.
%
% NOTE: "timekeeping" is a single word in English (no whitespace).
%   /https://www.thefreedictionary.com/timekeeping
%   Therefore capitalizing it only once as an identifier.
%
%
% EXAMPLE LOG OUTPUT
% ==================
% 2024-01-05T15:38:20 -- DEBUG -- SPEED -- CODE_NAME: Start
% 2024-01-05T15:38:20 -- DEBUG -- SPEED -- CODE_NAME: Stop: 0.000648 [s] (wall time)
% 2024-01-05T15:38:20 -- DEBUG -- SPEED -- CODE_NAME: Start
% 2024-01-05T15:38:20 -- DEBUG -- SPEED -- CODE_NAME: Stop: 0.000746 [s] (wall time), 7.46e-05 [s/gadget], 10 [gadgets], 3.73e-05 [s/bin], 20 [bins], 7.46e-07 [s/byte], 1000 [bytes]
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Timekeeper < handle
  % PROPOSAL: Shorter default strings.
  %   PROPOSAL: Remove "SPEED"
  % PROBLEM: Not specifying "wall time" for all measures of time.
  %   CON: Implicit.
  % PROPOSAL: Make possible to add arbitrary log values to stop log message.
  %   Ex: Saturation detection: CDF records per flag change



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    LOG_LEVEL = 'debug'
  end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Access=private)
    L
    tTicToc
    timeKeepingIsRunning
    codeName
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % codeName
    %       String. Some kind of brief name for the code section being
    %       tested, e.g. function name. Is used for logging.
    % L
    %       bicas.Logger object.
    function obj = Timekeeper(codeName, L)
      assert(ischar(codeName))
      assert(isa(L, 'bicas.Logger'))

      obj.L                    = L;
      obj.timeKeepingIsRunning = true;
      obj.codeName             = codeName;

      obj.log(sprintf('%s: Start', codeName))

      % IMPLEMENTATION NOTE: Should call tic() as late as possible to make
      % the timekeeping more accurate.
      obj.tTicToc = tic();
    end



    % Stop the timekeeping and log the result, including time per some unit
    % (and total number of units), e.g. bins, records, snapshots, bytes.
    % Object is unusable after calling this method.
    %
    %
    % ARGUMENTS
    % =========
    % varargin
    %       {i*2-1} = nUnits
    %       {i*2  } = unitName
    %       nUnits
    %           Arbitrary positive quantity which has increased from zero
    %           during the timekeeping. Should effectively quantify the
    %           magnitude of (some subset of) processing that has taken
    %           place during the timekeeping.
    %       unitName
    %           String. Unit of argument "nUnits" in singular. "s" will be
    %           added automatically for plural.
    %
    function stop_log(obj, varargin)
      % IMPLEMENTATION NOTE: Should call toc() as early as possible to
      % make the timekeeping more accurate.
      %
      % PROPOSAL: Log over multiple rows, one per quantity.
      %   NOTE: Total time still needs special treatment. ==> Separate row.
      %   CON: Two rows even when just one quantity.
      wallTimeSec = toc(obj.tTicToc);

      assert(obj.timeKeepingIsRunning);
      obj.timeKeepingIsRunning = false;

      % ARGUMENT ASSERTIONS
      % Even number of vargargin arguments.
      assert(mod(numel(varargin), 2) == 0)
      nQuantities = numel(varargin) / 2;

      sCa = {};
      for i = 1:nQuantities
        nUnits   = varargin{i*2 - 1};
        unitName = varargin{i*2 - 0};

        assert(isnumeric(nUnits))
        assert(nUnits >= 0)
        assert(ischar(unitName))

        % IMPLEMENTATION NOTE: nUnits might be an integer. Must convert
        % to double for division to work.
        %   Ex: bicas.proc.dsr.get_downsampling_bins(): nBins
        nUnits = double(nUnits);

        % NOTE: Adds "s" after unit to get plural.
        sCa{end+1} = sprintf('%g [s/%s], %g [%ss]', ...
          wallTimeSec/nUnits, unitName, nUnits, unitName);
      end
      s = [', ', strjoin(sCa, ', ')];

      obj.log_stop_timekeeping(wallTimeSec, s)
    end



  end    % methods(Access=public)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Method for logging. All log messages should pass though this function
    % rather than calling obj.L directly.
    function log(obj, s)
      % NOTE: The output adds "SPEED" in a way that resembles the logging
      % prefix convention specified by the RCS ICD. This makes it easier
      % to filter the log file w.r.t. speed information.
      obj.L.logf(...
        bicas.utils.Timekeeper.LOG_LEVEL, ...
        sprintf('SPEED -- %s', s))
    end



    % NOTE: Does not stop the timekeeping.
    % ARGUMENTS
    % =========
    % sExtra
    %       String that is appended directly after end of default string.
    %       Must either be empty or beginning with appropriate prefix e.g.
    %       comma+whitespace to make sense.
    function log_stop_timekeeping(obj, wallTimeSec, sExtra)
      s = sprintf(...
        '%s: Stop: %g [s] (wall time)%s', ...
        obj.codeName, wallTimeSec, sExtra);

      obj.log(s)
    end



  end    % methods(Access=private)



end
