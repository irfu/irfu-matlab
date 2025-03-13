%
% Stores one channel of data (array of samples; nRecords x nSpr) plus VSTBs
% (nRecords x 1). The class itself emulates a *column* array, both for CWF and
% SWF(!) data to make reconstruction of missing channels more natural. Every row
% represents data for a CDF record (for a given channel).
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
  %   samples (not only samples)
  %   data
  %   PROPOSAL: Analogous to dictionary class.
  %
  % PROPOSAL: Rename vsqbAr/vstbAr if none is correct.
  %
  % PROPOSAL: Automatic test code: plus(), minus()
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
    % IMPLEMENTATION NOTE: Not completely in accordance with definition to call
    % this VSQB (quality bit in datasets), but it is probably better than VSTB
    % (threshold bit), or at least if windowing is done before this stage. In
    % short, there is no satisfying pre-defined abbreviation for this bit.
    vsqbAr
  end
  properties (Dependent)
    % Nx1 array. Logical.
    % NOTE: Must have same size as object (column array), despite being a
    % function of samplesAr.
    bFp
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



    function bFp = get.bFp(obj)
      bFp = any(isnan(obj.samplesAr), 2);
      assert(iscolumn(bFp))
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
