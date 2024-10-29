%
% Stores one channel of data (samples) plus VSTBs. The class models a *column*
% array, both for CWF and SWF(!) data to make reconstruction of missing
% channels more natural. Every row represents data for a CDF record (for a
% given channel).
%
% SD = Source/Destination?
%      Signal Destination? (as in SDID)
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
    % NxM array. double. NaN represents missing data.
    samplesAr

    % Nx1 array.
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
          assert(isscalar(S))
          assert(isscalar(S(1).subs))

          ib = S(1).subs{1};
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
