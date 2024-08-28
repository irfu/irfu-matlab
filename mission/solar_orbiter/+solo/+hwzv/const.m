%
% SolO-related constants.
%
% HWZV = HardWare, zVariables
%
% NOTE: Not to be confused with BICAS's constants.
%
% sampere = set ampere
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-26
%
classdef const



  properties(Constant)

    % LFR sampling frequencies.
    LFR_F0_HZ = 24576;
    LFR_F1_HZ =  4096;
    LFR_F2_HZ =   256;
    LFR_F3_HZ =    16;

    % LSF = LFR Sampling Frequencies: F0, F1, F2, F3
    % The string names (F[0-3]) follow LFR's naming scheme.
    LSF_HZ = [...
      solo.hwzv.const.LFR_F0_HZ, ...
      solo.hwzv.const.LFR_F1_HZ, ...
      solo.hwzv.const.LFR_F2_HZ, ...
      solo.hwzv.const.LFR_F3_HZ];
    LSF_NAME_ARRAY = {'F0', 'F1', 'F2', 'F3'};

    % Should at least refer to the "normal" LFR snapshot length that BICAS
    % uses. Notes imply that there may be other ones (calibration? LFR-HF?
    % LFR-SCM?).
    LFR_SWF_SNAPSHOT_LENGTH = 2048;

    TDS_RSWF_SNAPSHOT_LENGTH_MIN = 2^10;
    TDS_RSWF_SNAPSHOT_LENGTH_MAX = 2^15;

    % Number of samples reserved for a snapshot in TDS (LFM) RSWF datasets.
    %
    % NOTE: This the max length of normal snapshots. The exception is
    % FULL_BAND mode snapshots which are 2^18=262144 samples/snapshot ("This
    % mode is however meant for calibration and testing, not for science
    % operations.").
    %
    % IMPORTANT NOTE
    % ==============
    % Empirically, this value varies between L1R and L1 datasets. Since
    % BICAS does not officially support reading L1 datasets, the code should
    % use the value for L1R. The code is also (probably) not adapted to
    % having another number of samples/record. This means that BICAS
    % currently (2023-10-09) can not read L1 SOLO_L1_RPW-TDS-LFM-RSWF
    % datasets.
    % Ex:
    % solo_L1R_rpw-tds-lfm-rswf-e-cdag_20200409_V12.cdf: 32768 samples/record
    % solo_L1_rpw-tds-lfm-rswf-cdag_20200409_V09.cdf   : 16384 samples/record
    TDS_RSWF_L1R_SAMPLES_PER_RECORD = 32768;

    % Max absolute value of set current.
    %
    % NOTE: Does not take into consideration that the actual min & max might
    % be slightly different due to that TM = -2^15 ... (2^15-1)., i.e.
    % max=(2^15-1)/2^15 * 60e-9 sampere.
    MAX_ABS_SAMPERE = 60e-6;

    TM_PER_SAMPERE = 32768 / solo.hwzv.const.MAX_ABS_SAMPERE;
  end



end
