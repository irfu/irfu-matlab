%
% Nominal implementation of superclass.
%
%
% IMPLEMENTATION NOTE
% ===================
% The class (instance methods, including constructor) deliberately does
% not itself read the RCTs, nor figure out which ones should be read.
% This is useful since
% ** it completely separates
%       (a) algorithms for determining RCTs to load, and
%       (b) reading RCT,
%    from the class (better modularization, better for automatic test
%    code).
% ** it makes it possible to inspect & modify the RCT content before
%    submitting it to bicas.proc.L1L2.cal.VoltageCalibration
% ** it simplifies the constructor.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VoltageCalibrationDataSupplierImpl < bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierAbstract
% PROPOSAL: Remove useZvcti2 as TDS argument since not supported.
%   NOTE: Constructor must still have it as argument since it is relevant for
%         LFR calibration.
%     PROPOSAL: useZvti2-->useLfrZvcti2.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % RCT calibration data
    Rctdc

    % Non-RCT calibration data
    % ------------------------
    % BIAS scalar (simplified) calibration, not in the RCTs. For
    % debugging/testing purposes.
    BiasScalarGain

    biasOffsetsDisabled
    useBiasTfScalar

    % Whether to use ZVCTI2 for calibration.
    useZvcti2
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = VoltageCalibrationDataSupplierImpl(Rctdc, useZvcti2, Bso)
      assert(isa(Rctdc, 'bicas.proc.L1L2.cal.RctdCollection'))
      assert(islogical(useZvcti2) & isscalar(useZvcti2))

      obj.Rctdc                          = Rctdc;

      obj.BiasScalarGain.alphaIvpav      = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV');
      obj.BiasScalarGain.betaIvpav       = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV');
      obj.BiasScalarGain.gammaIvpav.achg = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN');
      obj.BiasScalarGain.gammaIvpav.aclg = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN');

      obj.biasOffsetsDisabled            = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED');
      obj.useZvcti2                      = useZvcti2;

      %-------------------------
      % Set obj.useBiasTfScalar
      %-------------------------
      settingBiasTf                      = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF');
      switch(settingBiasTf)
        case 'FULL'
          obj.useBiasTfScalar = false;
        case 'SCALAR'
          obj.useBiasTfScalar = true;
        otherwise
          error(...
            'BICAS:Assertion:ConfigurationBug', ...
            ['Illegal value for setting',...
            ' PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF="%s".'], ...
            settingBiasTf)
      end
    end



    % Return BIAS ITF and BIAS offset.
    %
    % NOTE: May return calibration values corresponding to scalar
    % calibration, depending on BSO:
    %
    %
    % ARGUMENTS
    % =========
    % isAchg
    %       NUMERIC value: 0=Off, 1=ON, or NaN=Value not known.
    %       IMPLEMENTATION NOTE: Needs value to represent that isAchg
    %       is unknown. Sometimes, if isAchg is unknown, then it is
    %       useful to process as usual since some of the data can still be
    %       derived/calibrated, so that the caller does not need to handle
    %       the special case.
    %
    function [itfAvpiv, offsetAvolt] = get_BIAS_ITF_and_offset(obj, ...
        ssid, isAchg, iCalibTimeL, iCalibTimeH)

      % PROPOSAL: Separate functions for TF and scalar calibration.

      % ASSERTION
      assert(bicas.proc.L1L2.const.is_SSID(ssid) & isscalar(ssid))
      assert(bicas.proc.L1L2.const.SSID_is_ASR(ssid))
      assert(isscalar(isAchg) && isnumeric(isAchg))
      assert(isscalar(iCalibTimeL))
      assert(isscalar(iCalibTimeH))



      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');
      BiasRctd   = BiasRctdCa{1};

      %###################################################################
      % kFtfIvpav = Multiplication factor "k" that represents/replaces the
      % (forward) transfer function for scalar calibration.
      asid         = bicas.proc.L1L2.const.SSID_ASR_to_ASID( ssid);
      asidCategory = bicas.proc.L1L2.const.get_ASID_category(asid);
      antennas     = bicas.proc.L1L2.const.get_ASID_antennas(asid);
      switch(asidCategory)

        case 'DC_SINGLE'

          % NOTE: List of ITFs for different times.
          itfAvpiv    = BiasRctd.ItfSet.dcSingleAvpiv{iCalibTimeL};
          kFtfIvpav   = obj.BiasScalarGain.alphaIvpav;
          offsetAvolt = BiasRctd.dcSingleOffsetsAVolt(iCalibTimeH, antennas);

        case 'DC_DIFF'

          itfAvpiv = BiasRctd.ItfSet.dcDiffAvpiv{iCalibTimeL};
          kFtfIvpav    = obj.BiasScalarGain.betaIvpav;
          if     isequal(antennas, [1,2]);   offsetAvolt = BiasRctd.DcDiffOffsets.E12AVolt(iCalibTimeH);
          elseif isequal(antennas, [1,3]);   offsetAvolt = BiasRctd.DcDiffOffsets.E13AVolt(iCalibTimeH);
          elseif isequal(antennas, [2,3]);   offsetAvolt = BiasRctd.DcDiffOffsets.E23AVolt(iCalibTimeH);
          else
            error('BICAS:Assertion:IllegalArgument', 'Illegal argument "ssid".');
          end

        case 'AC_DIFF'

          if     isAchg == 0
            itfAvpiv    = BiasRctd.ItfSet.aclgAvpiv{iCalibTimeL};
            kFtfIvpav   = obj.BiasScalarGain.gammaIvpav.aclg;
            offsetAvolt = 0;
          elseif isAchg == 1
            itfAvpiv    = BiasRctd.ItfSet.achgAvpiv{iCalibTimeL};
            kFtfIvpav   = obj.BiasScalarGain.gammaIvpav.achg;
            offsetAvolt = 0;
          elseif isnan(isAchg)
            % CASE: AC GAIN unknown when it is NEEDED (i.e. when AC data).
            itfAvpiv    = bicas.const.NAN_TF;
            kFtfIvpav   = NaN;
            offsetAvolt = NaN;
          else
            error('BICAS:Assertion:IllegalArgument', ...
              'Illegal argument isAchg=%g.', isAchg)
          end

        otherwise
          error('BICAS:Assertion:IllegalArgument', ...
            ['Illegal argument "ssid" implies illegal ASID category="%s".', ...
            ' Can not obtain calibration data for this type of signal.'], ...
            asidCategory)
      end

      %============================================
      % Modify return values for specific settings
      %============================================
      if obj.useBiasTfScalar
        % NOTE: Overwrites "itfAvpiv".
        itfAvpiv = @(omegaRps) (ones(size(omegaRps)) / kFtfIvpav);
      end
      if obj.biasOffsetsDisabled
        % NOTE: Overwrites "offsetAvolt".
        offsetAvolt = 0;
      end
    end



    % Return LFR ITF, but also handle the case that should never happen for
    % actual non-NaN data (LSF F3 + BLTS 4 or 5) and return an ITF that only
    % returns NaN instead. BICAS may still iterate over that combination when
    % calibrating.
    %
    % RETURN VALUE
    % ============
    % itfIvpt
    %       LFR ITF, TM-->ivolt
    %
    function itfIvpt = get_LFR_ITF(obj, iBlts, iLsf, iNonBiasRct, zvcti2)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      bicas.proc.L1L2.cal.utils.assert_iLsf( iLsf)
      assert(isscalar(iNonBiasRct))
      assert(iNonBiasRct >= 1, 'Illegal iNonBiasRct=%g', iNonBiasRct)
      % No assertion on zvcti2 unless used (determined later).

      %==================================================
      % The only place to potentially make use of zvcti2
      %==================================================
      if obj.useZvcti2
        % ASSERTIONS
        assert(isscalar(zvcti2), ...
          'BICAS:IllegalArgument:Assertion', ...
          'Argument zvcti2 is not scalar.')
        assert(zvcti2 >= 0, ...
          'BICAS:IllegalArgument:Assertion', ...
          ['Illegal argument zvcti2=%g', ...
          ' (=zVar CALIBRATION_TABLE_INDEX(iRecord, 2))'], ...
          zvcti2)
        assert(iLsf == zvcti2+1, ...
          'BICAS:IllegalArgument:Assertion', ...
          'zvcti2+1=%i != iLsf=%i (before overwriting iLsf)', ...
          zvcti2+1, iLsf)

        % NOTE: Override earlier iLsf.
        % NOTE: This is the only place zvcti2 is used in this class.
        iLsf = zvcti2 + 1;
      end



      if (iLsf == 4) && ismember(iBlts, [4,5])
        % CASE: F3 and BLTS={4,5}

        % IMPLEMENTATION NOTE: There is no tabulated LFR TF for this case and
        % the h/w does not support it, so an accurate TF can not be returned
        % even in principle. However, the BICAS implementation (2025-07-11)
        % iterates over all 5x BLTS channels no matter the value of LSF, e.g.
        % for F3 (iLsf=4) when there is only real data on BLTS1-3. It can
        % therefore request "calibration" for this case anyway, even if it means
        % calibrating an empty channel (converting NaN values to NaN). For this
        % reason, the code can not raise an exception for this case.
        itfIvpt = bicas.const.NAN_TF;
      else
        LfrRctdCa = obj.Rctdc.get_RCTD_CA('LFR');

        % ASSERTION
        % IMPLEMENTATION NOTE: Anonymous function below will fail at a
        % later stage if these assertions are false. Checking for these
        % criteria here makes it easier to understand these particular
        % types of error.
        assert(iNonBiasRct <= numel(LfrRctdCa), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa is too small for argument iLfrRctd=%g.', ...
          ' This could indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
          ' value is larger than glob. attr. CALIBRATION TABLE allows.'], ...
          iNonBiasRct)
        assert(~isempty(LfrRctdCa{iNonBiasRct}), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa contains no RCT data corresponding', ...
          ' to argument iNonBiasRct=%g. This may indicate that', ...
          ' a zVar CALIBRATION_TABLE_INDEX(:,1) value is wrong or', ...
          ' that BICAS did not try to load the corresponding RCT', ...
          ' in glob. attr. CALIBRATION_TABLE.'], ...
          iNonBiasRct)

        itfIvpt = LfrRctdCa{iNonBiasRct}.ItfModifIvptCaCa{iLsf}{iBlts};
      end
    end    % get_LFR_ITF()



    % Return TDS-CWF nominal ITF as stored in RCT, except also handle BLTS 4-5.
    %
    % ARGUMENTS AND RETURN VALUES
    % ===========================
    % See get_LFR_ITF().
    %
    function itfIvpt = get_TDS_CWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(iNonBiasRct >= 1)
      if obj.useZvcti2
        % TODO? ASSERTION: zvcti2 = 0???
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-CWF calibration never uses ZVCTI2.')
      end



      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        RctdCa        = obj.Rctdc.get_RCTD_CA('TDS-CWF');
        tdsFactorIvpt = RctdCa{iNonBiasRct}.factorsIvpt(iBlts);

        itfIvpt       = @(omegaRps) (tdsFactorIvpt * ones(size(omegaRps)));
      else
        % CASE: BLTS 4-5 which TDS does NOT support (forbidden in h/w).

        itfIvpt = bicas.const.NAN_TF;
      end
    end    % get_TDS_CWF_ITF()



    % Return TDS-RSWF nominal ITF as stored in RCT, except also handle BLTS 4-5.
    %
    % ARGUMENTS AND RETURN VALUES
    % ===========================
    % See get_LFR_ITF().
    %
    function itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, iNonBiasRct, zvcti2)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(iNonBiasRct >= 1)
      if obj.useZvcti2
        % TODO? ASSERTION: zvcti2 = 0???
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-RSWF calibration never uses ZVCTI2.')
      end



      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        %======================================
        % Create combined ITF for TDS and BIAS
        %======================================
        RctdCa  = obj.Rctdc.get_RCTD_CA('TDS-RSWF');
        itfIvpt = RctdCa{iNonBiasRct}.itfModifIvptCa{iBlts};

      else
        % CASE: BLTS 4-5 which TDS does not support (forbidden in h/w).

        itfIvpt = bicas.const.NAN_TF;
      end
    end    % get_TDS_RSWF_ITF()



    function iCalibL = get_BIAS_calibration_time_index_L(obj, Epoch)
      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');

      iCalibL = bicas.proc.L1L2.cal.utils.get_calibration_time_index(...
        Epoch, BiasRctdCa{1}.epochL);
    end



    function iCalibH = get_BIAS_calibration_time_index_H(obj, Epoch)
      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');

      iCalibH = bicas.proc.L1L2.cal.utils.get_calibration_time_index(...
        Epoch, BiasRctdCa{1}.epochH);
    end



  end    % methods(Access=public)



end
