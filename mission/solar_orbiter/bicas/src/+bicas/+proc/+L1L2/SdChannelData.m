%
% Class which stores
% * one channel of data (array of samples; nRecords x nSpr), and
% * VSQBs (nRecords x 1).
% for one SDID/SSID value.
%
% The class interface itself (~syntactic sugar) emulates a column array with
% support for addition and subtraction despite that the underlying storage of
% samples is a 2D array. Every row represents data for a CDF record (for a given
% channel). This makes reconstruction of missing channels more natural while the
% class itself automatically derives new VSQBs (one per row) for reconstructed
% channels under the hood. Since the class is meant to be used for
% reconstructing missing channels, it is meant to represent data both before and
% after reconstruction, and possibly demultiplexing.
%
%
% SD = Source/Destination?
%      Signal Destination? (as in SDID)
%
% NOTE: Effectively uses NaN as fill values (provides function for it).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SdChannelData
  % PROPOSAL: Better name.
  %   ~ASR
  %   ~SDID
  %   Signal Destination
  %   channel
  %   *one* channel, single channel
  %   samples (NOTE: Contains not only samples; contains VSQB too.)
  %   data
  %   metadata
  %   demultiplexing input & output
  %     CON: Might also use class for deriving saturation.
  %   DIOC = DemultiplexingIOChannel
  %   DCIO = DemultiplexingChannelIO
  %   CHND = ChannelData
  %   OCD = OneChannelData
  %   SCD, SCHD = SingleChannelData
  %     SCD is substring of SCDA, isCdag.
  %   SCH = SingleChannel
  %     SCH is substring of ischar.
  %   --
  %   NOTE: If changing name of this class, then should also change name of
  %         SdcdDict=bicas.proc.L1L2.SdChannelDataDict.
  %
  % PROPOSAL: Rename samplesAr to include unit.
  %   PRO: Class intended for being used for reconstruction.
  %   CON: Might be used when deriving saturation before calibration in the
  %        future.
  %
  % TODO-NI: Performance for large arrays? Need internal handle objects?
  %   Cf. bicas.utils.FPArray.
  %
  % PROPOSAL: Use FPAs for samples.
  %   CON(?): The size of the object is not the same as the size of samplesAr.
  %           Can therefore not *directly* reuse FPA fill positions as fill
  %           positions for this class.
  %
  % PROPOSAL: Include length of each snapshot.
  %   PRO: No longer requires variable-length snapshots to be padded with NaN.
  %   NOTE: Requires checking that lengths of snapshots are consistent when
  %         combining multiple objects (plus/minus).
  % PROPOSAL: Add SDID/SSID.
  %   PRO: SSID (source) used when deriving VSIBs (when calling
  %         bicas.proc.L1L2.dc.get_VSIB_5xBLTS_NEW()).
  %   PRO: SDID (destination) used when reconstructing.
  %   PRO: SDID/SSID och SDCD är ofta argument tillsammans till funktioner.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private)
    % NxM array. double. NaN represents missing data. Supports both CWF (Nx1)
    % and SWF (NxM). Could in principle be used for samples in either TM units
    % or calibrated units.
    samplesAr

    % Nx1 array.
    vsibAr
  end
  properties(Dependent)
    % Number of rows with at least one NaN in the underlying data.
    % Nx1 array. Logical.
    %
    % NOTE: Must have same size as object (column array), despite being a
    % function of samplesAr.
    bWholeRowIsNan
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods
    % TODO-NI: Why is this within in a separate method statement? Seems to be
    %          public anyway. Because it defines the behaviour of a dependant
    %          field variable?



    function bWholeRowIsNan = get.bWholeRowIsNan(obj)
      % IMPLEMENTATION NOTE: Must require ALL elements to be NaN in order to
      % make reconstruction algorithm work correctly due to using NaN for
      % padding snapshots. TDS-RSWF snapshots are variable-length meaning that
      % the unused elements are NaN, but those unused elements can not be used
      % for determining whether (1) a channel record should be reconstructed
      % (assigned a value), or (2) be used for reconstructing other channels
      % (read the value).

      % PROPOSAL: Abolish function/property. Derive when invoking reconstruction
      %           algorithm instead.
      %   CON: Not obvious what property represents given that there are
      %        multiple samples per property value.

      bWholeRowIsNan = all(isnan(obj.samplesAr), 2);
      assert(iscolumn(bWholeRowIsNan))
    end



  end
  methods(Access=public)



    function obj = SdChannelData(samplesAr, vsibAr)
      assert(isfloat(samplesAr))
      assert(islogical(vsibAr))

      irf.assert.sizes(...
        samplesAr, [-1, NaN], ...
        vsibAr,    [-1])

      obj.samplesAr = samplesAr;
      obj.vsibAr    = vsibAr;
    end



    % Indexing overloading: Array indexing for reading.
    function varargout = subsref(obj, S)
      switch S(1).type
        case '()'
          % Index object as if it were a column vector.
          assert(isscalar(S))
          assert(isscalar(S(1).subs))

          ib        = S(1).subs{1};
          samplesAr = obj.samplesAr(ib, :);
          vsibAr    = obj.vsibAr(   ib, :);
          % IMPLEMENTATION NOTE: Specifying ":" for second index for vsibAr is
          % necessary for ensuring always returning a column vector, despite
          % that it is a column vector already.

          varargout = {bicas.proc.L1L2.SdChannelData(samplesAr, vsibAr)};

        case '.'
          % Call method (sic!)
          [varargout{1:nargout}] = builtin('subsref', obj, S);

        otherwise
          error('BICAS:Assertion', 'Unsupported operation.')
      end
    end



    % Indexing overloading: Array indexing for writing: Sdcd(i) = ...
    %
    %
    % PERFORMANCE
    % ===========
    % TODO: Investigate. Cf. bicas.utils.FPArray.
    %
    function Sdcd1 = subsasgn(Sdcd1, S, Sdcd2)
      assert(isa(Sdcd2, 'bicas.proc.L1L2.SdChannelData'))

      switch S(1).type
        case '()'
          assert(isscalar(S))
          assert(isscalar(S(1).subs))

          ib = S(1).subs{1};

          Sdcd1.samplesAr(ib, :) = Sdcd2.samplesAr;
          Sdcd1.vsibAr(   ib)    = Sdcd2.vsibAr;

        otherwise
          error('BICAS:Assertion', 'Unsupported operation.')
      end
    end



    % "Overload" size(obj, ...)
    function s = size(obj, varargin)
      s = size(obj.vsibAr, varargin{:});
    end



    % Operator overloading.
    function Sdcd3 = plus(Sdcd1, Sdcd2)
      samplesAr3 = Sdcd1.samplesAr + Sdcd2.samplesAr;
      vsibAr3    = Sdcd1.vsibAr    | Sdcd2.vsibAr;

      Sdcd3 = bicas.proc.L1L2.SdChannelData(samplesAr3, vsibAr3);
    end



    % Operator overloading.
    function Sdcd3 = minus(Sdcd1, Sdcd2)
      samplesAr3 = Sdcd1.samplesAr - Sdcd2.samplesAr;
      vsibAr3    = Sdcd1.vsibAr    | Sdcd2.vsibAr;

      Sdcd3 = bicas.proc.L1L2.SdChannelData(samplesAr3, vsibAr3);
    end



  end    % methods(Access=public)



end
