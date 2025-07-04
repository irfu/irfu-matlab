%
% Class which stores
% * one channel of data (array of samples; nRecords x nSpr), and
% * VSIBs (nRecords x 1).
% for one SDID/SSID value.
%
% The class interface itself (~syntactic sugar) emulates a column array with
% support for addition and subtraction despite that the underlying storage of
% samples is a 2D array. Every row represents data for a CDF record (for a given
% channel). This makes reconstruction of missing channels more natural while the
% class itself automatically derives new VSIBs (one per row) for reconstructed
% channels under the hood. Since the class is meant to be used for
% reconstructing missing channels, it is meant to represent data both before and
% after reconstruction, and possibly demultiplexing.
%
%
% NOTE: Effectively uses NaN as fill values (provides function for it).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SingleChannelData
  % PROPOSAL: Better name.
  %   ~ASR
  %   ~SDID
  %   Signal Destination
  %   channel
  %   single channel
  %   samples (NOTE: Contains not only samples; contains VSIB too.)
  %   data
  %   metadata
  %   demultiplexing input & output
  %     CON: Might also use class for deriving saturation.
  %   --
  %   DIOC = DemultiplexingIOChannel
  %   DCIO = DemultiplexingChannelIO
  %   CHD, CHND = ChannelData
  %   OCD = OneChannelData
  %   SCD, SCHD = SingleChannelData -- IMPLEMENTED
  %     SCD is substring of SCDA, isCdag.
  %     CON: Too generic?
  %   SCH = SingleChannel
  %     SCH is substring of ischar.
  %   --
  %   NOTE: If changing name of this class, then should also change name of
  %         Cdac=bicas.proc.L1L2.ChannelDataAsrCollection.
  %
  % TODO-NI: Performance for large arrays? Need internal handle objects?
  %   Cf. bicas.utils.FPArray.
  %   PROPOSAL: Class itself as handle object?
  %
  % PROPOSAL: Use FPAs for samples.
  %   CON(?): The size of the object is not the same as the size of samplesAr.
  %           Can therefore not *directly* reuse FPA fill positions as fill
  %           positions for this class.
  %
  % PROPOSAL: Include length of each snapshot.
  %   PRO: No longer requires variable-length snapshots to be padded with NaN.
  %   CON: Uses more memory (duplicates information since snapshot lengths are
  %        the same for all channels).
  %   NOTE: Requires checking that lengths of snapshots are consistent when
  %         combining multiple objects (plus/minus).
  %   CON-PROPOSAL: Generic class for jagged array, for elements of "any"
  %                 MATLAB class.
  %     CON: Can not handle fill values like FPAs.
  %       PRO: Can not replace with FPAs in the long term unless
  %            (1) the class itself uses FPAs, or
  %            (2) implements fill positions itself.
  %
  % PROPOSAL: Add SDID/SSID (one value).
  %   PRO: SSID (source) is used when calibrating.
  %   PRO: SSID (source) is used when deriving VSIBs (when calling
  %         bicas.proc.L1L2.dc.get_VSIB_5xBLTS_NEW()).
  %   PRO: SDID (destination) used when reconstructing missing channel data.
  %   PRO: SDID/SSID and SCHD are often together as function arguments.
  %   CON: In CDAC: SDID is constant (within every SCHD separately) and
  %        specified externally in the corresponding CDAC key.
  %   PROPOSAL: Add both SSID and SDID.
  %     NOTE: SSID is unknown for reconstructed data.



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



    function obj = SingleChannelData(samplesAr, vsibAr)
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

          varargout = {bicas.proc.L1L2.SingleChannelData(samplesAr, vsibAr)};

        case '.'
          % Call method (sic!)
          [varargout{1:nargout}] = builtin('subsref', obj, S);

        otherwise
          error('BICAS:Assertion', 'Unsupported operation.')
      end
    end



    % Indexing overloading: Array indexing for writing: Schd(i) = ...
    %
    %
    % PERFORMANCE
    % ===========
    % TODO: Investigate. Cf. bicas.utils.FPArray.
    %
    function Schd1 = subsasgn(Schd1, S, Schd2)
      assert(isa(Schd2, 'bicas.proc.L1L2.SingleChannelData'))

      switch S(1).type
        case '()'
          assert(isscalar(S))
          assert(isscalar(S(1).subs))

          ib = S(1).subs{1};

          Schd1.samplesAr(ib, :) = Schd2.samplesAr;
          Schd1.vsibAr(   ib)    = Schd2.vsibAr;

        otherwise
          error('BICAS:Assertion', 'Unsupported operation.')
      end
    end



    % "Overload" size(obj, ...)
    function s = size(obj, varargin)
      s = size(obj.vsibAr, varargin{:});
    end



    % Operator overloading.
    function Schd3 = plus(Schd1, Schd2)
      samplesAr3 = Schd1.samplesAr + Schd2.samplesAr;
      vsibAr3    = Schd1.vsibAr    | Schd2.vsibAr;

      Schd3 = bicas.proc.L1L2.SingleChannelData(samplesAr3, vsibAr3);
    end



    % Operator overloading.
    function Schd3 = minus(Schd1, Schd2)
      samplesAr3 = Schd1.samplesAr - Schd2.samplesAr;
      vsibAr3    = Schd1.vsibAr    | Schd2.vsibAr;

      Schd3 = bicas.proc.L1L2.SingleChannelData(samplesAr3, vsibAr3);
    end



  end    % methods(Access=public)



end
