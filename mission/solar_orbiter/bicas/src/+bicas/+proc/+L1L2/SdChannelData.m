%
% Class which stores
% * one channel of data (array of samples; nRecords x nSpr), and
% * VSQBs (nRecords x 1).
%
% The class interface itself (~syntactic sugar) emulates a column array with
% support for addition and subtraction despite that the underlying storage of
% samples is a 2D array. Every row represents data for a CDF record (for a given
% channel). This makes reconstruction of missing channels more natural while the
% class itself automatically derives new VSQBs (one per row) for reconstructed
% channels under the hood. Since the class is meant to be used for
% demultiplexing, it is meant to represent data both before and demultiplexing.
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
  %   samples (not only samples; contains vsqb too)
  %   data
  %   metadata
  %   PROPOSAL: Analogous to dictionary class.
  %
  % PROPOSAL: Rename vsqbAr/vstbAr if none is correct.
  %   PROPOSAL: VSB = Voltage Saturation Bit
  %             VSIB = Voltage Saturation Implementation/Intermediate Bit
  %   PROPOSAL: Globally replace VSTB-->VSIB=Voltage Saturation Intermediate Bit
  %
  % TODO-NI: Performance for large arrays? Need internal handle objects?
  %   Cf. bicas.utils.FPArray.
  %
  % TODO: Methods for deriving/reconstructing other channels.
  %   Needs (scalar?) SSID?
  %   PROPOSAL: Implement using operator overloading: +, - (only two needed).
  %
  % PROPOSAL: Use FPAs for samples.
  %   CON(?): The size of the object is not the same as the size of samplesAr.
  %           Can therefore not *directly* reuse FPA fill positions as fill
  %           positions for this class.
  %
  % PROBLEM: Current implementation leads to bug in reconstruction algorithm:
  %          Emulating column array while bIsNan represents there being at least
  %          one NaN per row.
  %   --
  %   PROPOSAL: Refactor to emulate same size as samplesAr (2D) but still only
  %             contain column array of saturation bits.
  %     PRO: Can abolish bIsNan.
  %     PRO: Needed for reconstruction algorithm: If row contains both non-NaN and
  %          NaN, then reconstruction algorithm thinks entire row is NaN, i.e.
  %          contains no data.
  %       Ex: TDS-RSWF snapshots do not fill entire row. The unused part is NaN.
  %     CON: Reconstruction algorithm requires linear indexing. ==> Intermediate
  %          SDCDs are linear!!!!
  %       Ex: A3(bDerive3) = fh12to3(A1(bDerive3), A2(bDerive3));
  %   --
  %   PROPOSAL: Refactor to emulate same size as samplesAr (2D) and keep 2D
  %             array of VSQBs (instead of column array).
  %     CON: Requires more memory.
  %     CON: Abandons idea that SDCD could contain other information on the form
  %          of one scalar per row (unless duplicates them to one per element).
  %   --
  %   PROPOSAL: Refactor to emulate column array but (1) include counters for
  %             valid number of samples per row, and (2) assumes/asserts that
  %             all (valid) samples on a row are either non-NaN or NaN.
  %     PRO: Records are "fundamental" and reconstruction should work on
  %          records.
  %       PRO: One saturation bit per record, makes records "fundamental".
  %       PRO: Arbitrary indexing operations, e.g. logical indexing, linear
  %            indexing make no sense unless index=record.
  %     CON: Must assert all valid samples on a row are either non-NaN or NaN.
  %     CON: When adding, subtracting: Must assert equal number of valid samples
  %          per row.
  %     CON: Makes class more conceptually complex.
  %       CON: Incorporates number of valid samples per row in class,
  %            and replaces the corresponding external variable?
  %         CON: One such variable per channel, not one globally.
  %     CON: Using variable for number of valid samples makes implementation
  %          resemble bicas.utils.FPArray.
  %     CON: Can not handle there being both fill values and non-fill values
  %          inside snapshot.
  %       CON: Should never happen.
  %       CON: Calibration should turn one fill value into fill values for
  %            entire snapshot.
  %       CON: If it happens, it could only be taken advantage of when there is
  %            redundant channel data anyway (rare).
  %   --
  %   PROPOSAL: bIsNan <=> All elements are NaN (not: at least one).
  %             -- IMPLEMENTED



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private)
    % NxM array. double. NaN represents missing data. Supports both CWF (Nx1)
    % and SWF (NxM).
    samplesAr

    % Nx1 array.
    % IMPLEMENTATION NOTE: Not completely in accordance with definition to name
    % this VSQB (quality bit in datasets), but it is probably better than VSTB
    % (threshold bit), or at least if windowing is done before this stage. In
    % short, there is no satisfying pre-defined abbreviation for this bit.
    vsqbAr
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
      % make reconstruction algorithm work correctly.
      % TDS-RSWF snapshots are variable-length meaning that the unused elements
      % are NaN, but those unused elements can not be used for determining
      % whether a channel record should be reconstructed (assigned a value) or
      % can be used for reconstructing other channels (read the value).

      % PROPOSAL: Abolish function/property. Derive when invoking reconstruction
      %           algorithm instead.
      %   CON: Not obvious what property represents given that there are
      %        multiple samples per property value.

      bWholeRowIsNan = all(isnan(obj.samplesAr), 2);
      assert(iscolumn(bWholeRowIsNan))
    end



  end
  methods(Access=public)



    function obj = SdChannelData(samplesAr, vsqbAr)
      assert(isfloat(samplesAr))
      assert(islogical(vsqbAr))

      irf.assert.sizes(...
        samplesAr, [-1, NaN], ...
        vsqbAr,    [-1])

      obj.samplesAr = samplesAr;
      obj.vsqbAr    = vsqbAr;
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
          vsqbAr    = obj.vsqbAr(   ib, :);
          % IMPLEMENTATION NOTE: Specifying ":" for second index for vsqbAr is
          % necessary for ensuring always returning a column vector, despite
          % that it is a column vector already.

          varargout = {bicas.proc.L1L2.SdChannelData(samplesAr, vsqbAr)};

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
          Sdcd1.vsqbAr(   ib)    = Sdcd2.vsqbAr;

        otherwise
          error('BICAS:Assertion', 'Unsupported operation.')
      end
    end



    % "Overload" size(obj, ...)
    function s = size(obj, varargin)
      s = size(obj.vsqbAr, varargin{:});
    end



    % Operator overloading.
    function Sdcd3 = plus(Sdcd1, Sdcd2)
      samplesAr3 = Sdcd1.samplesAr + Sdcd2.samplesAr;
      vsqbAr3    = Sdcd1.vsqbAr    | Sdcd2.vsqbAr;

      Sdcd3 = bicas.proc.L1L2.SdChannelData(samplesAr3, vsqbAr3);
    end



    % Operator overloading.
    function Sdcd3 = minus(Sdcd1, Sdcd2)
      samplesAr3 = Sdcd1.samplesAr - Sdcd2.samplesAr;
      vsqbAr3    = Sdcd1.vsqbAr    | Sdcd2.vsqbAr;

      Sdcd3 = bicas.proc.L1L2.SdChannelData(samplesAr3, vsqbAr3);
    end



  end    % methods(Access=public)



end
