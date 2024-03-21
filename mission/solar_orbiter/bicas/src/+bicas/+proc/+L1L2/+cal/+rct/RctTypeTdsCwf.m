%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctTypeTdsCwf < bicas.proc.L1L2.cal.rct.RctType



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant, GetAccess=public)
    filenameRegexpSettingKey = 'PROCESSING.RCT_REGEXP.TDS-LFM-CWF';
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % NOTE: TDS CWF cwfFactorsIvpt are already inverted (can be seen from
    % units).
    function RctData = read_RCT(filePath)

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

        RctData = [];
        RctData.factorsIvpt = factorsIvpt;

      catch Exc1
        Exc2 = MException(...
          'BICAS:FailedToReadInterpretRCT', ...
          ['Error when interpreting calibration file (TDS team''s', ...
          ' LFM CWF RCT for BIAS/BICAS) "%s"'], ...
          filePath);
        Exc2 = Exc2.addCause(Exc1);
        throw(Exc2);
      end
    end



    % NOTE: Contains no TFs. Data is therefore trivial to use as it is in
    % the RCT.
    function RctData2 = modify_RCT_data(RctData1)
      RctData2 = RctData1;
    end



    function log_RCT(RctData, L)

      L.logf(bicas.proc.L1L2.cal.rct.RctType.RCT_DATA_LL, ...
        'TDS CWF calibration factors: %s [ivolt/TM]', ...
        bicas.proc.L1L2.cal.utils.vector_string('%g', RctData.factorsIvpt));
    end



  end    % methods(Static)



end
