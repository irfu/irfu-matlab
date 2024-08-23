%
% Store of all input information for processing that is common for all
% L1/L1R-->L2 LFR+TDS processing (i.e. demultiplexing and calibration).
%
%
% IMPLEMENTATION NOTE: bltsSamplesTm and lrx
% ==========================================
% DCIP always represents (has variables/elements for) all five BLTS's,
% despite that only three are used at any given time. The channels not used
% are set to NaN. Which ones are actually used can be switched at any given
% time due to lrx changing.
%
% This handling does in principle differ from the handling of other "data
% parameters" (BW, LSF, isTdsCwf, freqHz, hasSwfFormat, etc.). LRX
% determines which channels contain actual data. The BLTS index into the
% array is used to determine whether one has DC diff or AC diff data.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DemultiplexingCalibrationInput
  % PROPOSAL: Shorten DemultiplexingCalibrationInput/Output
  %   PROPOSAL: Demultiplexing -->
  %             Demuxing
  %   PROPOSAL: Calibration -->
  %             Calib
  %   DemuxingCalibrationInput/Output
  %   DemuxingCalibInput/Output
  %
  % PROPOSAL: Refactor code to only represent three data channels (those
  %           BLTS's which (nominally) contain actual data for the given
  %           timestamp). Store iBlts for the channels that represent BLTS 2/4
  %           and BLTS 3/5 so that one can use the demultiplexer routings.
  %   CON: Has to split these channels/arrays to set output dataset ZVs which
  %        are greater in number.
  %       CON: Should be easy to implement.
  %   CON: This increases the assumptions which the ~demultiplexer code and
  %        the DCIP class make.
  %       CON: No, it does not. Which?
  %           CON: There are always 3 channels worth of data. BLTS 2-3/4-5 are
  %                different from BLTS 1
  %                and can change meaning due to LRX.
  %   CON: Must invent a new concept (type of index) to represent these three
  %        data channels.
  %       Ex: New abbreviation/name.
  %       Ex: Convention for indexing individual channels.
  %
  %   PROPOSAL: Split BLTSs into two separate groups (variables):
  %       (1) Always BLTS 1:                   Zv.blts1Sampmles   (Nx1)
  %       (2) Shifts between BLTS 2-3 and 4-5: Zv.blts2345Samples (Nx2)
  %       NOTE: BLTS 2-3 ARE NOT ALWAYS DIFFS!!! (They are always DC though.)
  %       CON: Can not iterate over all 5 BLTSs.
  %       TODO-DEC: Name of index BLTS 2-3/4-5 variable?
  %       PROPOSAL: Separate iBltsArray (Nx2) for labelling the data source.
  %           NOTE: Does not want to store Routing objects for each timestamp
  %                 but using iBltsArray, that can be derived when needed.
  %   PROPOSAL: Keep channes in single array (Nx3).
  %       PROPOSAL: Separate iBltsArray (Nx3) for labelling the data source.
  %           CON: iBltsArray(:, 1) == 1 always.
  %
  % PROPOSAL: Separate representations (variables) for
  %           (1) CWF and (2) SWF (snapshots) data.
  %   PRO: Could clarify code.
  %   PROBLEM: How handle the absent form of data? There will always be
  %            timestamps and other ZVs.
  %   PROPOSAL: Variable for unused data is empty, e.g. 0x0.
  %       CON: Violates convention/assumption that all Zv.* variables have
  %            same number of rows.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    Zv
    Ga
    hasSwfFormat
    isLfr
    isTdsCwf
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)

    function obj = DemultiplexingCalibrationInput(Zv, Ga, hasSwfFormat, isLfr, isTdsCwf)

      irf.assert.struct(Zv, ...
        {'Epoch', 'bltsSamplesTm', 'freqHz', 'nValidSamplesPerRecord', ...
        'bdmFpa', 'isAchgFpa', 'dlrFpa', ...
        'iLsf', ...
        'QUALITY_BITMASK', 'QUALITY_FLAG', 'SYNCHRO_FLAG', ...
        'DELTA_PLUS_MINUS', 'CALIBRATION_TABLE_INDEX', ...
        'ufv', 'lrx', 'BW'}, {});
      bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(Zv);
      assert(size(Zv.bltsSamplesTm, 3) == 5)
      assert(isa(Zv.freqHz,    'double' ))
      assert(isa(hasSwfFormat, 'logical'))
      assert(isa(isLfr,        'logical'))
      assert(isa(isTdsCwf,     'logical'))

      obj.Zv           = Zv;
      obj.Ga           = Ga;
      obj.hasSwfFormat = hasSwfFormat;
      obj.isLfr        = isLfr;
      obj.isTdsCwf     = isTdsCwf;
    end

  end    % methods(Access=public)

end
