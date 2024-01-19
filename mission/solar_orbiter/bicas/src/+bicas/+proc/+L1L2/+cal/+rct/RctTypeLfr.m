%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctTypeLfr < bicas.proc.L1L2.cal.rct.RctType



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant, GetAccess=public)
    filenameRegexpSettingKey = 'PROCESSING.RCT_REGEXP.LFR';
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % LfrFtfTpivTable : {iLsf}{iBlts}. Table of LFR FTFs.
    %                   iLsf=1..3 : iBlts=1..5 for BLTS 1-5.
    %                   iLsf=4    : iBlts=1..3 for BIAS 1-3.
    function RctData = read_RCT(filePath)
      Do = dataobj(filePath);

      try
        % ASSUMPTION: Exactly 1 CDF record.
        % IMPLEMENTATION NOTE: Does not want to rely one dataobj special
        % behaviour for 1 record case
        % ==> Remove leading singleton dimensions, much assertions.

        % NOTE: There are separate TFs for each BLTS channel, not just
        % separate LFR sampling frequencies, i.e. there are 5+5+5+3=18
        % TFs (but only 1 frequency table/LSF, since they are recycled).
        % NOTE: The assignment of indices here effectively determines
        % the translation between array index and LFR Sampling Frequency
        % (LSF). This is NOT the same as the values in the LFR zVar
        % FREQ.
        freqTableHz{1}   = shiftdim(Do.data.Freqs_F0.data);    % NOTE: Index {iLsf}.
        freqTableHz{2}   = shiftdim(Do.data.Freqs_F1.data);
        freqTableHz{3}   = shiftdim(Do.data.Freqs_F2.data);
        freqTableHz{4}   = shiftdim(Do.data.Freqs_F3.data);

        amplTableTpiv{1} = shiftdim(Do.data.TF_BIAS_12345_amplitude_F0.data);
        amplTableTpiv{2} = shiftdim(Do.data.TF_BIAS_12345_amplitude_F1.data);
        amplTableTpiv{3} = shiftdim(Do.data.TF_BIAS_12345_amplitude_F2.data);
        amplTableTpiv{4} = shiftdim(Do.data.TF_BIAS_123_amplitude_F3.data);

        phaseTableDeg{1} = shiftdim(Do.data.TF_BIAS_12345_phase_F0.data);
        phaseTableDeg{2} = shiftdim(Do.data.TF_BIAS_12345_phase_F1.data);
        phaseTableDeg{3} = shiftdim(Do.data.TF_BIAS_12345_phase_F2.data);
        phaseTableDeg{4} = shiftdim(Do.data.TF_BIAS_123_phase_F3.data);

        for iLsf = 1:4
          % NOTE: F3 is an exception and has no AC (iBlts={4,5}) TF.
          if iLsf ~= 4   nBlts = 5;
          else           nBlts = 3;
          end

          % NOTE: Values for the specific LSF, hence the prefix.
          lsfFreqTableHz   = freqTableHz{  iLsf};
          lsfAmplTableTpiv = amplTableTpiv{iLsf};
          lsfPhaseTableDeg = phaseTableDeg{iLsf};

          % ASSERTIONS: Check CDF array sizes, and implicitly that the
          % CDF format is the expected one.
          nFreqs = irf.assert.sizes(...
            lsfFreqTableHz,   [-1,       1 ], ...
            lsfAmplTableTpiv, [-1, nBlts], ...
            lsfPhaseTableDeg, [-1, nBlts]);
          assert(nFreqs >= bicas.proc.L1L2.cal.rct.RctType.TF_TABLE_MIN_LENGTH)

          for iBlts = 1:nBlts

            lsfBltsFreqTableHz   = lsfFreqTableHz;
            lsfBltsAmplTableTpiv = lsfAmplTableTpiv(:, iBlts);
            lsfBltsPhaseTableDeg = lsfPhaseTableDeg(:, iBlts);

            FtfTpiv = irf.utils.tabulated_transform(...
              lsfBltsFreqTableHz * 2*pi, ...
              lsfBltsAmplTableTpiv, ...
              deg2rad(lsfBltsPhaseTableDeg));

            % ASSERTION: FTF
            assert(FtfTpiv.toward_zero_at_high_freq())

            FtfTpivTable{iLsf}{iBlts} = FtfTpiv;
          end
        end

        % NOTE: Storing data in struct field to clarify the nature of
        % the content to the caller.
        RctData = [];
        RctData.FtfTpivTable = FtfTpivTable;

      catch Exc1
        Exc2 = MException(...
          'BICAS:FailedToReadInterpretRCT', ...
          ['Error when interpreting calibration file', ...
          ' (LFR team''s RCT for BIAS/BICAS) "%s"'], ...
          filePath);
        Exc2 = Exc2.addCause(Exc1);
        throw(Exc2);
      end
    end



    function RctData2 = modify_RCT_data(RctData1)

      FtfRctTpivCaCa = RctData1.FtfTpivTable;

      % Read LFR FTFs, derive ITFs and modify them.
      itfModifIvptCaCa = {};
      for iLsf = 1:numel(FtfRctTpivCaCa)

        itfModifIvptCaCa{end+1} = {};
        for iBlts = 1:numel(FtfRctTpivCaCa{iLsf})

          % INVERT: tabulated FTF --> tabulated ITF
          ItfIvpt = FtfRctTpivCaCa{iLsf}{iBlts}.inverse();

          % MODIFY tabulated ITF
          ItfModifIvpt = bicas.proc.L1L2.cal.utils.extrapolate_tabulated_TF_to_zero_Hz(ItfIvpt);

          % MODIFY tabulated ITF --> Function TF
          %
          % NOTE: Can not blindly forbid extrapolation (beyond the
          % extrapolation to 0 Hz already done above) by setting value
          % outside table=NaN (which deliberately triggers error
          % elsewhere). LFR's tabulated TFs do in theory cover
          % frequencies up to the Nyquist frequency, but in practice,
          % the actual sampling frequency varies slightly. This means
          % that when the Nyquist frequency also varies slightly and
          % sometimes it exceeds the tabulated frequencies.
          % Ex: solo_L1R_rpw-lfr-surv-cwf-e-cdag_20201102_V01.cd
          %
          % 2020-11-06: LFR tables (RCT):
          % F0=24576 Hz: f={  12--12288} [Hz]
          % F1= 4096 Hz: f={0.01-- 2048} [Hz]
          % F2=  256 Hz: f={0.01--  128} [Hz]
          % F3=   16 Hz: f={0.01--    8} [Hz]
          VALUE_OUTSIDE_TABLE = 0;
          %VALUE_OUTSIDE_TABLE = NaN;   % Does not work. See comments above.
          itfModifIvpt = @(omegaRps) (bicas.proc.L1L2.cal.utils.eval_tabulated_TF(...
            ItfModifIvpt, omegaRps, VALUE_OUTSIDE_TABLE));
          clear ItfModifIvpt

          itfModifIvptCaCa{iLsf}{iBlts} = itfModifIvpt;
        end
      end

      RctData2 = [];
      % NOTE: RctData.FtfRctTpivCaCa is still kept (for debugging).
      RctData2.FtfRctTpivCaCa   = FtfRctTpivCaCa;    % Just copied.
      RctData2.ItfModifIvptCaCa = itfModifIvptCaCa;
    end



    function log_RCT(RctData, L)
      % NOTE: Frequencies may go outside of tabulated data.
      %FREQ_HZ = [0, 1, 5];
      FREQ_HZ = [0, 1, 100];

      % CASE: This index corresponds to an actually loaded RCT (some are
      % intentionally empty).
      for iLsf = 1:4
        if iLsf ~= 4   nBltsMax = 5;
        else           nBltsMax = 3;
        end

        for iBlts = 1:nBltsMax

          itfNamePrefix = sprintf('LFR, F%i, BLTS/BIAS_%i', iLsf-1, iBlts);

          itfName = sprintf('%s FTF (as in RCT)', itfNamePrefix);
          bicas.proc.L1L2.cal.utils.log_TF_tabulated(...
            bicas.proc.L1L2.cal.rct.RctType.RCT_DATA_LL, ...
            itfName, ...
            RctData.FtfRctTpivCaCa{iLsf}{iBlts}, ...
            L);

          itfIvpt          = RctData.ItfModifIvptCaCa{iLsf}{iBlts};
          itfName = sprintf('%s ITF (modif., interp.)', itfNamePrefix);

          bicas.proc.L1L2.cal.utils.log_TF_function_handle(...
            bicas.proc.L1L2.cal.rct.RctType.RCT_DATA_LL, ...
            itfName, ...
            'ivolt/TM unit', FREQ_HZ, itfIvpt, L)

        end
      end    % for

    end



  end    % methods(Static)



end
