%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctDataTdsCwf < bicas.proc.L1L2.cal.rct.RctData



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    factorsIvpt
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = RctDataTdsCwf(filePath)
      obj@bicas.proc.L1L2.cal.rct.RctData(filePath)

      % NOTE: RCT contains no TFs and data is therefore trivial to use as it is
      % in the RCT.
      RctRawData = bicas.proc.L1L2.cal.rct.RctDataTdsCwf.read_RCT(filePath);
      obj.factorsIvpt = RctRawData.factorsIvpt;
    end



    function log_RCT(obj, L)

      L.logf(bicas.proc.L1L2.cal.rct.RctData.RCT_DATA_LL, ...
        'TDS CWF calibration factors: %s [ivolt/TM]', ...
        bicas.proc.L1L2.cal.utils.vector_string('%g', obj.factorsIvpt));
    end



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % NOTE: TDS CWF cwfFactorsIvpt are already inverted (can be seen from
    % units).
    function RctRawData = read_RCT(filePath)

      Do = dataobj(filePath);

      try
        % NOTE: Undocumented in CDF: zVar CALIBRATION_TABLE is
        % volt/count for just multiplying the TDS signal (for this kind
        % of data). Is not a frequency-dependent transfer function.

        % ASSUMPTION: Exactly 1 CDF record.
        % IMPLEMENTATION NOTE: Does not want to rely on dataobj's
        % special behaviour for 1 record case
        % ==> Remove leading singleton dimensions, much assertions.

        factorsIvpt = shiftdim(Do.data.CALIBRATION_TABLE.data);

        % ASSERTIONS: Check CDF array sizes, no change in format.
        irf.assert.sizes(factorsIvpt, [3,1])

        D = [];
        D.factorsIvpt = factorsIvpt;

      catch Exc1
        Exc2 = MException(...
          'BICAS:FailedToReadInterpretRCT', ...
          ['Error when interpreting calibration file (TDS team''s', ...
          ' LFM CWF RCT for BIAS/BICAS) "%s"'], ...
          filePath);
        Exc2 = Exc2.addCause(Exc1);
        throw(Exc2);
      end

      RctRawData = D;
    end



  end    % methods(Static)



end
