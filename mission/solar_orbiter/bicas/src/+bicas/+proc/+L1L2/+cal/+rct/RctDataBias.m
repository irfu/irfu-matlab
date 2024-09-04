%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctDataBias < bicas.proc.L1L2.cal.rct.RctData



  %##########################
  %##########################
  % PUBLIC STATIC PROPERTIES
  %##########################
  %##########################
  % properties(Constant)
  %   RCTTID = 'BIAS'
  % end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant, GetAccess=private)

    % Minimum number of numerator or denominator coefficients in the BIAS RCT.
    N_MIN_TF_NUMER_DENOM_COEFFS = 8;
  end
  properties(SetAccess=immutable)
    epochL
    epochH
    Current
    FtfRctSet
    ItfSet
    dcSingleOffsetsAVolt
    DcDiffOffsets
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = RctDataBias(filePath)
      obj@bicas.proc.L1L2.cal.rct.RctData(filePath)

      RctRawData = bicas.proc.L1L2.cal.rct.RctDataBias.read_RCT(filePath);

      %=============================================
      % Modify file data and store it in the object
      %=============================================
      FtfRctSet = RctRawData.FtfSet;

      % ASSERTIONS
      nTime = irf.assert.sizes(...
        FtfRctSet.DcSingleAvpiv, [-1, 1], ...
        FtfRctSet.DcDiffAvpiv,   [-1, 1], ...
        FtfRctSet.AclgAvpiv,     [-1, 1], ...
        FtfRctSet.AchgAvpiv,     [-1, 1]);

      % NOTE: Derive ITFs.
      ItfSet = [];
      for iTf = 1:nTime
        % INVERT: FTF --> ITF

        % Temporary variables which are stored in the definitions of
        % anonymous functions later.
        % * Might speed up code by eliminating calls to method .inverse()
        % * Reduces size of individual expressions.
        TempItfDcSingleAvpiv = FtfRctSet.DcSingleAvpiv{iTf}.inverse();
        TempItfDcDiffAvpiv   = FtfRctSet.DcDiffAvpiv  {iTf}.inverse();
        TempItfAclgAvpiv     = FtfRctSet.AclgAvpiv    {iTf}.inverse();
        TempItfAchgAvpiv     = FtfRctSet.AchgAvpiv    {iTf}.inverse();

        ItfSet.dcSingleAvpiv{iTf} = @(omegaRps) (TempItfDcSingleAvpiv.eval(omegaRps));
        ItfSet.dcDiffAvpiv  {iTf} = @(omegaRps) (TempItfDcDiffAvpiv  .eval(omegaRps));
        ItfSet.aclgAvpiv    {iTf} = @(omegaRps) (TempItfAclgAvpiv    .eval(omegaRps));
        ItfSet.achgAvpiv    {iTf} = @(omegaRps) (TempItfAchgAvpiv    .eval(omegaRps));
      end

      obj.epochL               = RctRawData.epochL;
      obj.epochH               = RctRawData.epochH;
      obj.Current              = RctRawData.Current;
      obj.FtfRctSet            = FtfRctSet;    % Change name of field (sic!).
      obj.ItfSet               = ItfSet;
      obj.dcSingleOffsetsAVolt = RctRawData.dcSingleOffsetsAVolt;
      obj.DcDiffOffsets        = RctRawData.DcDiffOffsets;
    end



    % Log some indicative value(s) for a BIAS RCT.
    %
    % NOTE: Does not log file path. Caller is assumed to do that.
    function log_RCT(obj, L)

      % Logging parameters
      DC_FREQ_HZ       = [0];   % Single & diffs.
      AC_DIFF_FREQS_HZ = [0, 1000];
      LL               = bicas.proc.L1L2.cal.rct.RctData.RCT_DATA_LL;

      %=====================
      % Iterate over EpochL
      %=====================
      for iEpochL = 1:numel(obj.epochL)

        L.logf(LL, 'Below values are used for data beginning %s:', ...
          bicas.utils.TT2000_to_UTC_str(obj.epochL(iEpochL)))

        % Log bias current calibration
        L.logf(LL, '    BIAS current offsets: %s [aampere]',         ...
          bicas.proc.L1L2.cal.utils.vector_string(...
          '% 10e', obj.Current.offsetsAAmpere(iEpochL, :)))
        L.logf(LL, '    BIAS current gain   : %s [aampere/TM unit]', ...
          bicas.proc.L1L2.cal.utils.vector_string(...
          '% 10e', obj.Current.gainsAapt(     iEpochL, :)))

        % Log transfer functions (frequency domain) at selected
        % frequencies.
        L.logf(LL, ...
          ['    Note: Not logging the exact RCT BIAS TFs', ...
          ' (FTFs; obj.FtfRctSet) since the inversion is trivial.'])
        log_TF('    BIAS ITF DC single',          DC_FREQ_HZ,       obj.ItfSet.dcSingleAvpiv)
        log_TF('    BIAS ITF DC diff',            DC_FREQ_HZ,       obj.ItfSet.dcDiffAvpiv)
        log_TF('    BIAS ITF AC diff, low  gain', AC_DIFF_FREQS_HZ, obj.ItfSet.aclgAvpiv)
        log_TF('    BIAS ITF AC diff, high gain', AC_DIFF_FREQS_HZ, obj.ItfSet.achgAvpiv)
      end

      %=====================
      % Iterate over EpochH
      %=====================
      % NOTE: Must work for multiple CDF records.
      dcDiffOffsetsAVolt = [...
        obj.DcDiffOffsets.E12AVolt, ...
        obj.DcDiffOffsets.E13AVolt, ...
        obj.DcDiffOffsets.E23AVolt];
      irf.assert.sizes(dcDiffOffsetsAVolt, [NaN, 3]);
      for iEpochH = 1:numel(obj.epochH)

        L.logf(LL, 'Below values are used for data beginning %s:', ...
          bicas.utils.TT2000_to_UTC_str(obj.epochH(iEpochH)))

        L.logf(LL, ...
          '    BIAS DC single voltage offsets ( V1, V2, V3): %s [avolt]', ...
          bicas.proc.L1L2.cal.utils.vector_string('%g', ...
          obj.dcSingleOffsetsAVolt(iEpochH, :)))
        L.logf(LL, ...
          '    BIAS DC diff   voltage offsets (E12,E13,E23): %s [avolt]', ...
          bicas.proc.L1L2.cal.utils.vector_string('%g', ...
          dcDiffOffsetsAVolt(iEpochH)))
      end

      %###################################################################
      % Nested utility function.
      % NOTE: Implicitly function of iEpochL, L, LL.
      function log_TF(name, freqArray, ItfList)
        bicas.proc.L1L2.cal.utils.log_TF_function_handle(...
          LL, name, 'avolt/ivolt', freqArray, ...
          ItfList{iEpochL}, L);
      end
      %###################################################################
    end



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function [RctRawData] = read_RCT(filePath)
      % TODO-DEC: How handle time?
      %   PROPOSAL: "Only" access the BIAS values (trans.func and other) through a function instead of selecting
      %             indices in a data struct.
      %       PROPOSAL: (private method) [omegaRps, zVpc] = get_transfer_func(epoch, signalType)
      %           signalType = 'DC_SINGLE' etc

      Do = dataobj(filePath);

      % Constants for interpreting the array indices in the CDF.
      I_NUMERATOR   = 1;
      I_DENOMINATOR = 2;
      %
      I_DC_SINGLE = 1;
      I_DC_DIFF   = 2;
      I_AC_LG     = 3;
      I_AC_HG     = 4;
      %
      I_E12 = 1;
      I_E13 = 2;
      I_E23 = 3;

      try
        % NOTE: Assumes 1 CDF record or many (time-dependent values).
        % ==> Must handle that dataobj assigns differently for these two
        %     cases.
        epochL                    = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.Epoch_L);
        epochH                    = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.Epoch_H);
        biasCurrentOffsetsAAmpere = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.BIAS_CURRENT_OFFSET);      % DEPEND_0 = Epoch_L
        biasCurrentGainsAapt      = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.BIAS_CURRENT_GAIN);        % DEPEND_0 = Epoch_L
        dcSingleOffsetsAVolt      = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.V_OFFSET);                 % DEPEND_0 = Epoch_H
        dcDiffOffsetsAVolt        = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.E_OFFSET);                 % DEPEND_0 = Epoch_H
        ftfCoeffs                 = bicas.proc.L1L2.cal.rct.RctDataBias.normalize_dataobj_ZV(Do.data.TRANSFER_FUNCTION_COEFFS); % DEPEND_0 = Epoch_L

        nEpochL = size(epochL, 1);
        nEpochH = size(epochH, 1);

        % IMPLEMENTATION NOTE: Corrects for what seems to be a bug in
        % dataobj. dataobj permutes/removes indices, and permutes them
        % differently depending on the number of CDF records (but wrong
        % in all cases).
        %
        % 1 CDF record : cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       F/TTT"   # 3=number of dimensions/record
        % 2 CDF records: cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       T/TTT"
        % 1 CDF record:  size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [4 2 8]
        % 2 CDF records: size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [2 4 2 8]
        ftfCoeffs = permute(ftfCoeffs, [1, 4, 3, 2]);



        %=======================================================
        % ASSERTIONS: Size of tfCoeffs/TRANSFER_FUNCTION_COEFFS
        %=======================================================
        % ND = Numerator Denominator
        nNdCoeffs = irf.assert.sizes(ftfCoeffs, [nEpochL, -1, 2, 4]);
        assert(nNdCoeffs >= bicas.proc.L1L2.cal.rct.RctDataBias.N_MIN_TF_NUMER_DENOM_COEFFS)

        %================================
        % Assign struct that is returned
        %================================
        D = struct();
        D.epochL = epochL;
        D.epochH = epochH;

        D.Current.offsetsAAmpere = biasCurrentOffsetsAAmpere;
        D.Current.gainsAapt      = biasCurrentGainsAapt;
        D.dcSingleOffsetsAVolt   = dcSingleOffsetsAVolt;
        D.DcDiffOffsets.E12AVolt = dcDiffOffsetsAVolt(:, I_E12);
        D.DcDiffOffsets.E13AVolt = dcDiffOffsetsAVolt(:, I_E13);
        D.DcDiffOffsets.E23AVolt = dcDiffOffsetsAVolt(:, I_E23);

        % NOTE: Using name "FtfSet" only to avoid "Ftfs" (plural).
        % (List, Table would be wrong? Use "FtfTable"?)
        D.FtfSet.DcSingleAvpiv = bicas.proc.L1L2.cal.rct.RctDataBias.create_TF_sequence(...
          ftfCoeffs(:, :, I_NUMERATOR,   I_DC_SINGLE), ...
          ftfCoeffs(:, :, I_DENOMINATOR, I_DC_SINGLE));

        D.FtfSet.DcDiffAvpiv = bicas.proc.L1L2.cal.rct.RctDataBias.create_TF_sequence(...
          ftfCoeffs(:, :, I_NUMERATOR,   I_DC_DIFF), ...
          ftfCoeffs(:, :, I_DENOMINATOR, I_DC_DIFF));

        D.FtfSet.AclgAvpiv = bicas.proc.L1L2.cal.rct.RctDataBias.create_TF_sequence(...
          ftfCoeffs(:, :, I_NUMERATOR,   I_AC_LG), ...
          ftfCoeffs(:, :, I_DENOMINATOR, I_AC_LG));

        D.FtfSet.AchgAvpiv = bicas.proc.L1L2.cal.rct.RctDataBias.create_TF_sequence(...
          ftfCoeffs(:, :, I_NUMERATOR,   I_AC_HG), ...
          ftfCoeffs(:, :, I_DENOMINATOR, I_AC_HG));

        % ASSERTIONS
        irf.assert.sizes(...
          D.FtfSet.DcSingleAvpiv, [nEpochL, 1], ...
          D.FtfSet.DcDiffAvpiv,   [nEpochL, 1], ...
          D.FtfSet.AclgAvpiv,     [nEpochL, 1], ...
          D.FtfSet.AchgAvpiv,     [nEpochL, 1]);
        for iEpochL = 1:nEpochL
          %assert(Bias.ItfSet.DcSingleAvpiv{iEpochL}.eval(0) > 0, 'BICAS:FailedToReadInterpretRCT', 'DC single inverted transfer function is not positive (and real) at 0 Hz. (Wrong sign?)');
          %assert(Bias.ItfSet.DcDiffAvpiv{iEpochL}.eval(0)   > 0, 'BICAS:FailedToReadInterpretRCT',   'DC diff inverted transfer function is not positive (and real) at 0 Hz. (Wrong sign?)');
          % Unsure if assertion makes sense for AC, or possibly even
          % for DC.
          % 2020-03-10: This criterion is not true for AC high-gain
          % transfer function fit now used (but does for AC diff
          % low-gain).
        end

        %==============================================================
        % ASSERTIONS:
        % All variables NOT based on tfCoeffs/TRANSFER_FUNCTION_COEFFS
        %==============================================================
        bicas.utils.assert_ZV_Epoch(D.epochL)
        bicas.utils.assert_ZV_Epoch(D.epochH)
        validateattributes(D.epochL, {'numeric'}, {'increasing'})
        validateattributes(D.epochH, {'numeric'}, {'increasing'})

        irf.assert.sizes(...
          D.Current.offsetsAAmpere, [nEpochL, 3], ...
          D.Current.gainsAapt,      [nEpochL, 3], ...
          D.dcSingleOffsetsAVolt,   [nEpochH, 3]);

        for fn = fieldnames(D.DcDiffOffsets)'
          irf.assert.sizes(D.DcDiffOffsets.(fn{1}), [nEpochH, 1]);
        end

      catch Exc1
        Exc2 = MException(...
          'BICAS:FailedToReadInterpretRCT', ...
          'Error when interpreting calibration file (BIAS RCT) "%s"', filePath);
        Exc2 = Exc2.addCause(Exc1);
        throw(Exc2)
      end

      RctRawData = D;
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Utility function
    %
    % Function for normalizing the indices of dataobj zVariables. dataobj
    % zVariable arrays have different meanings for their indices depending
    % on whether there are one record or many. If there is one record, then
    % there is not record index. If there are multiple records, then the
    % first index represents the record number. This function inserts a
    % size-one index as the first index.
    %
    % DO   = dataobj(...)
    % data = Do.data.TRANSFER_FUNCTION_COEFFS.data
    %
    % NOTE: Not well tested on different types of zvar array sizes.
    %
    function data = normalize_dataobj_ZV(DataobjZv)
      % PROPOSAL: Move to utils.
      % PROPOSAL: Shorter name:
      %   norm_dataobj_zv
      %   norm_do_zv

      data = DataobjZv.data;

      if DataobjZv.nrec == 1
        data = shiftdim(data, -1);
      end
    end



    % Utility function to simplify read_BIAS_RCT. Arguments correspond to
    % zVariables in BIAS RCT.
    %
    %
    % ARGUMENTS
    % =========
    % ftfNumCoeffs,
    % ftfDenomCoeffs : 2D matrix of numerator/denominator coefficients for
    %                  a sequence of FTFs. (iTime, iCoeff).
    %
    %
    % RETURN VALUE
    % ============
    % FtfArray : 1D column cell array of FTFs
    %            (irf.utils.rational_func_transform).
    %
    function FtfArray = create_TF_sequence(ftfNumCoeffs, ftfDenomCoeffs)

      % ASSERTIONS
      nTime = irf.assert.sizes(...
        ftfNumCoeffs,   [-1, -2], ...
        ftfDenomCoeffs, [-1, -2]);
      %assert(size(ftfNumCoeffs, 1) == size(ftfDenomCoeffs, 1))
      % The last FTF denominator coefficient (highest index, for which the
      % value is non-zero) must be =1.
      assert(...
        ftfDenomCoeffs(find(ftfDenomCoeffs, 1, 'last')) == 1, ...
        'BICAS:FailedToReadInterpretRCT', ...
        ['RCT should contain forward transfer function (FTF)', ...
        ' denominator coefficients,', ...
        ' where the highest-order (non-zero) coefficient', ...
        ' is the number one (1).', ...
        ' The data does not satisfy this criterion.'])

      FtfArray = {};
      for iTime = 1:nTime

        Ftf = irf.utils.rational_func_transform(...
          ftfNumCoeffs(  iTime, :), ...
          ftfDenomCoeffs(iTime, :));

        % ASSERTIONS
        assert(Ftf.has_real_impulse_response())
        % Assert FTF. Can not set proper error message.
        assert(Ftf.zero_in_high_freq_limit(), ...
          'BICAS:FailedToReadInterpretRCT', ...
          ['Transfer function is expected to be "forward",', ...
          ' i.e. in the direction of the physical signal.', ...
          ' It seems not to be.'])

        FtfArray{end+1, 1} = Ftf;    % Force column array.
      end
    end



  end    % methods(Static, Access=private)



end
