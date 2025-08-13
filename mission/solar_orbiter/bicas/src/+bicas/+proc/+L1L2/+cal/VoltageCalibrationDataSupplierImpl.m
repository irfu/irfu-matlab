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
  % PROPOSAL: Remove useNbci as TDS argument since not supported.
  %   NOTE: Constructor must still have it as argument since it is relevant for
  %         LFR calibration.
  %   PROPOSAL: Move to VCAL.
  %     PRO: This class should not do logic and is not meant to be tested.
  %
  % PROPOSAL: Separate calibration data-retrieving methods
  %     get_BIAS_ITF_and_offset()
  %     get_LFR_ITF()
  %     get_TDS_CWF_ITF()
  %     get_TDS_RSWF_ITF()
  %   into
  %     (1) "core": crash for cases when there is no data,
  %     and
  %     (2) "wrapper": call "core", but return NaN TF etc for cases when there is
  %                    no data.
  %   Move "wrapper" to VCAL.
  %   PRO: Can test wrapper functionality.
  %     Ex: CALIBRATION_TABLE_INDEX(i, 1)=wrong when BW=0/1.



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

    % Whether to use NBCI for calibration.
    useNbci
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = VoltageCalibrationDataSupplierImpl(Rctdc, useNbci, Bso)
      assert(isa(Rctdc, 'bicas.proc.L1L2.cal.RctdCollection'))
      assert(islogical(useNbci) & isscalar(useNbci))

      obj.Rctdc                          = Rctdc;

      % NOTE: Below values are FTF.
      obj.BiasScalarGain.alphaIvpav      = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV');
      obj.BiasScalarGain.betaIvpav       = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV');
      obj.BiasScalarGain.gammaIvpav.achg = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN');
      obj.BiasScalarGain.gammaIvpav.aclg = Bso.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN');

      obj.useNbci                        = useNbci;
    end



    function [itfAvpiv, kItfAvpiv, offsetAvolt] = get_BIAS_ITF_and_offset(obj, ...
        ssid, isAchg, iCalibTimeL, iCalibTimeH)
      % PROPOSAL: Convert isAchg (double!) to FPA.
      %   PRO: Can handle fill/unknown/missing value.

      % ASSERTION
      assert(bicas.proc.L1L2.const.is_SSID(ssid) & isscalar(ssid))
      assert(bicas.proc.L1L2.const.SSID_is_ASR(ssid))
      assert(isscalar(isAchg) && isnumeric(isAchg))
      assert(isscalar(iCalibTimeL))
      assert(isscalar(iCalibTimeH))



      BiasRctdCa = obj.Rctdc.get_RCTD_CA('BIAS');
      BiasRctd   = BiasRctdCa{1};

      asid         = bicas.proc.L1L2.const.SSID_ASR_to_ASID( ssid);
      asidCategory = bicas.proc.L1L2.const.get_ASID_category(asid);
      antennas     = bicas.proc.L1L2.const.get_ASID_antennas(asid);
      switch(asidCategory)

        case 'DC_SINGLE'

          % NOTE: List of ITFs for different times.
          itfAvpiv    = BiasRctd.ItfSet.dcSingleAvpiv{iCalibTimeL};
          kFtfIvpav   = obj.BiasScalarGain.alphaIvpav;
          offsetAvolt = BiasRctd.dcSingleOffsetsAvolt(iCalibTimeH, antennas);

        case 'DC_DIFF'

          itfAvpiv = BiasRctd.ItfSet.dcDiffAvpiv{iCalibTimeL};
          kFtfIvpav    = obj.BiasScalarGain.betaIvpav;
          if     isequal(antennas, [1,2]);   offsetAvolt = BiasRctd.DcDiffOffsets.E12Avolt(iCalibTimeH);
          elseif isequal(antennas, [1,3]);   offsetAvolt = BiasRctd.DcDiffOffsets.E13Avolt(iCalibTimeH);
          elseif isequal(antennas, [2,3]);   offsetAvolt = BiasRctd.DcDiffOffsets.E23Avolt(iCalibTimeH);
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

      kItfAvpiv = 1 / kFtfIvpav;
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
    function itfIvpt = get_LFR_ITF(obj, iBlts, NbriFpa, NbciFpa, iLsf)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      bicas.proc.L1L2.cal.utils.assert_iLsf( iLsf)
      assert(isscalar(NbriFpa))
      assert(isscalar(NbciFpa))

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

        Rctd = obj.get_RCTD(NbriFpa, 'LFR');

        if isempty(Rctd)
          itfIvpt = bicas.const.NAN_TF;
        else
          %==================================================================
          % The only place to potentially make use of NbciFpa
          %==================================================================
          % IMPLEMENTATION NOTE: Executing this after that a RCTD has been
          % identified and first when this value is actually needed to simplify
          % the assertions on NbciFpa.
          if obj.useNbci
            % ASSERTIONS
            assert(isscalar(NbciFpa))
            assert(~NbciFpa.fpAr)

            nbci = NbciFpa.array();

            assert(nbci >= 0)
            % assert(iLsf == nbci+1, ...
            %   'BICAS:IllegalArgument:Assertion', ...
            %   'nbci=%i != iLsf=%i (before overwriting iLsf)', ...
            %   nbci, iLsf)

            % NOTE: Override earlier iLsf.
            % NOTE: This is the only place where nbci is used in
            %       this class.
            iLsf = nbci + 1;
          end

          itfIvpt = Rctd.ItfModifIvptCaCa{iLsf}{iBlts};
        end

      end
    end    % get_LFR_ITF()



    % Utility function for reducing code.  Mostly for asserting that
    % NbriFpa is valid and does not imply a bad
    % CALIBRATION_TABLE_INDEX(i,1) value.
    %
    function Rctd = get_RCTD(obj, NbriFpa, rcttid)
      assert(isscalar(NbriFpa))

      RctdCa = obj.Rctdc.get_RCTD_CA(rcttid);

      if NbriFpa.fpAr
        % CASE: There is no RCTD for the given nbri
        % ------------------------------------------------
        % This should ALWAYS imply that LFR ZV BW=0 (BIAS h/W is off) since
        % preceding code has set it to fill value.
        % Ex: solo_L1R_rpw-lfr-surv-swf-e-cdag_20200225_V10.cdf
        % Ex: solo_L1R_rpw-lfr-surv-swf-e-cdag_20200228_V10.cdf

        Rctd = [];
      else
        % CASE: NbriFpa is scalar AND not a fill position.

        nbri = NbriFpa.array();

        % ASSERTIONS
        % Assert that the value of ZV CALIBRATION_TABLE_INDEX(i, 1) is correct.
        assert(nbri >= 1, 'Illegal nbri=%g', nbri)
        assert(nbri <= numel(RctdCa), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR/TDA RctdCa is too small for argument nbri=%g.', ...
          ' This could indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
          ' value is larger than glob. attr. CALIBRATION TABLE allows.'], ...
          nbri)
        assert(~isempty(RctdCa{nbri}), ...
          'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
          ['LFR LfrRctdCa contains no RCT data corresponding', ...
          ' to argument nbri=%g. This may indicate that', ...
          ' zVariable CALIBRATION_TABLE_INDEX(:,1) value is wrong or', ...
          ' that BICAS did not try to load the corresponding RCT', ...
          ' specified in glob. attr. CALIBRATION_TABLE.'], ...
          nbri)

        Rctd = RctdCa{nbri};
      end

    end



    % Return TDS-CWF nominal ITF as stored in RCT, except also handle BLTS 4-5.
    %
    % ARGUMENTS AND RETURN VALUES
    % ===========================
    % See get_LFR_ITF().
    %
    function itfIvpt = get_TDS_CWF_ITF(obj, iBlts, NbriFpa, NbciFpa)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(~NbriFpa.fpAr, 'Illegal CALIBRATION_TABLE_INDEX(i,1) value.')
      nbri = NbriFpa.array();
      assert(nbri >= 1)

      if obj.useNbci
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-CWF calibration never uses NBCI.')
      end



      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        RctdCa        = obj.Rctdc.get_RCTD_CA('TDS-CWF');
        tdsFactorIvpt = RctdCa{nbri}.factorsIvpt(iBlts);

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
    function itfIvpt = get_TDS_RSWF_ITF(obj, iBlts, NbriFpa, NbciFpa)

      % ASSERTIONS
      bicas.proc.L1L2.cal.utils.assert_iBlts(iBlts)
      assert(~NbriFpa.fpAr, 'Illegal CALIBRATION_TABLE_INDEX(i,1) value.')
      nbri = NbriFpa.array();
      assert(nbri >= 1)
      if obj.useNbci
        error(...
          'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
          'TDS-RSWF calibration never uses NBCI.')
      end



      if ismember(iBlts, [1,2,3])
        % CASE: BLTS 1-3 which TDS does support.

        %======================================
        % Create combined ITF for TDS and BIAS
        %======================================
        RctdCa  = obj.Rctdc.get_RCTD_CA('TDS-RSWF');
        itfIvpt = RctdCa{nbri}.itfModifIvptCa{iBlts};

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
