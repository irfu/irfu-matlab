%
% Class that collects "processing functions" as public static methods. Only
% covers processing L2-->L3.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like"
%   data, use the first MATLAB array index to represent CDF records.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% See readme.txt.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef L2L3
%
% PROPOSAL: Split up processing between (a) density, and (b) E-field & SCPOT.
%   PRO: Faster
%       CON: Not very heavy operation.
%   PRO: Leads to better organization of code.
%       PRO: process_L2_to_L3() is too large and should be split up anyway.
%   CON: DENSITY is a function EFIELD+SCPOT, and thus has to be processed after
%        the latter.
%
% NOTE: Class only has one function.
%   PROPOSAL: Convert to function file.



    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################

    
    
        % Processing function for processing L2-->L3 (not VHT).
        %
        % NOTE: Function assumes that (some) fill values for integer-valued
        % zVariables are identical in input and output datasets.
        %
        % NOTE: Function does not discard data with QUALITY_FLAG==fill value, as
        % opposed to QUALITY_FLAG < threshold.
        %
        % NOTE: Sets QUALITY_FLAG==fill value when ALL data in record is NaN.
        % Both OSR and DSR. The same is not(?) enforced in L2 processing, but
        % should maybe be. /EJ 2021-05-12
        %
        function [OutEfieldOsr,  OutEfieldDsr, ...
                  OutScpotOsr,   OutScpotDsr, ...
                  OutDensityOsr, OutDensityDsr] ...
                = process_L2_to_L3(InLfrCwf, Ec, SETTINGS, L)

            % PROPOSAL: Split up into different parts for EFIELD, SCPOT, DENSITY
            %           (still combine non-downsampled and downsampled).
            %   CON: Slows down overall processing.
            %       PRO: Must read same L2 dataset multiple times.
            %       PRO: Must read L3 SCPOT dataset to produce L3 DENSITY dataset.
            %   CON: There is much shared functionality for 3 quality ZVs.
            %       PRO: Same ~constants
            %           Ex: INPUT_DSI, BIN_LENGTH_WOLS_NS, BIN_TIMESTAMP_POS_WOLS_NS
            %       PRO: Read setting QUALITY_FLAG_MIN_FOR_USE
            %       PRO: Normalizing CWF zVar names.
            %       PRO: Preparations for downsampled.
            %           Bin locations, bundling of records,
            %           Downsampling of quality variables
            %               (QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK).
            %           DELTA_PLUS_MINUS_dsr
            %
            % NOTE: ROC BUG: https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/48
            %         L1 QUALITY_BITMASK seems to use the wrong value (255) as
            %         fill value (FILLVAL=65535). ==> A BICAS bug fix would not
            %         fix the entire issue.
            %   PROPOSAL: Use double also for CDF integer variables so NaN can
            %             represent fill value also for these.
            %   PROPOSAL: Implement MATLAB equivalent of the JUICE pipeline's
            %             FPA class.
            %
            % NOTE: L2 LFR-CWF-E skt previously had zVar
            %   QUALITY_BITMASK=CDF_UINT1, fill value=255 (wrong)
            % until skt V12 when it was changed to
            %   QUALITY_BITMASK=CDF_UINT2, fill value 65535 (correct).

            tTicToc = tic();
            assert(isa(Ec, 'bicas.proc.L2L3.ExternalCodeAbstract'))



            %===========
            % Constants
            %===========
            % The only acceptable input DSI.
            INPUT_DSI                 = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
            % NS = Nanoseconds
            BIN_LENGTH_WOLS_NS        = int64(10e9);
            BIN_TIMESTAMP_POS_WOLS_NS = int64(BIN_LENGTH_WOLS_NS / 2);
            


            %======================
            % Normalize zVar names
            %======================
            [InLfrCwf.Zv, fnChangeList] = ...
                irf.ds.normalize_struct_fieldnames(InLfrCwf.Zv, ...
                {{{'VDC', 'V'}, 'VDC'}}, 'Assert one matching candidate');

            bicas.proc.utils.handle_ZV_name_change(...
                fnChangeList, INPUT_DSI, SETTINGS, L, 'VDC', ...
                'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')
            
            
            
            %=======================================
            % Call BICAS-external code to calculate
            % (1) EFIELD, SCPOT, and from that
            % (2) DENSITY.
            %=======================================
            LfrCwfZv = [];
            LfrCwfZv.Epoch            = InLfrCwf.Zv.Epoch;
            LfrCwfZv.VDC              = InLfrCwf.Zv.VDC;
            LfrCwfZv.EDC              = InLfrCwf.Zv.EDC;
            LfrCwfZv.QUALITY_FLAG_Fpa = InLfrCwf.ZvFpa.QUALITY_FLAG;
            R = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT_DENSITY(LfrCwfZv, Ec, SETTINGS);


            
            %===================================================================
            % ~HACK: MODIFY INPUT ARGUMENT InLfrCwf
            % -------------------------------------
            % IMPLEMENTATION NOTE: This is to modify QUALITY_FLAG for both OSR
            % and DSR datasets. In principle, this is for keeping the interface
            % to bicas.proc.dsr.init_shared_DSR_ZVs() simple.
            %===================================================================
            InLfrCwf.ZvFpa.QUALITY_FLAG(R.bNotUsed) = bicas.utils.FPArray.FP_UINT8;



            %====================================================================
            % Derive values for CDF global attribute "Misc_calibration_versions"
            %====================================================================
            % NOTE: Should not add BICAS version to GA
            % "Misc_calibration_versions" since this is already encoded in GA
            % "Software_version" (together with "Software_name").
            %
            % NOTE: Density "Misc_calibration_versions" contains all three
            % versions, since density is derived from PSP.
            vdccalStr    = ['solo.vdccal() code version ',     R.vdccalCodeVerStr];
            vdccalMatStr = ['solo.vdccal() calibration file ', R.vdccalMatVerStr];
            psp2neStr    = ['solo.psp2ne() code version ',     R.psp2neCodeVerStr];
            gaEfieldScpot_Misc_calibration_versions = {vdccalStr, vdccalMatStr};
            gaDensity_Misc_calibration_versions     = ...
                [gaEfieldScpot_Misc_calibration_versions, {psp2neStr}];

            
            
            %=========================================
            % Misc. variables shared between datasets
            %=========================================
            % Global attributes -- shared between all OSR+DSR datasets.
            InitialGa = struct();
            InitialGa.OBS_ID                = InLfrCwf.Ga.OBS_ID;
            InitialGa.SOOP_TYPE             = InLfrCwf.Ga.SOOP_TYPE;
            % zVariables -- shared between all OSR datasets.
            InitialOsrZv = struct();
            InitialOsrZv.Epoch              = InLfrCwf.Zv.Epoch;
            InitialOsrZv.QUALITY_FLAG       = InLfrCwf.ZvFpa.QUALITY_FLAG;
            InitialOsrZv.QUALITY_BITMASK    = InLfrCwf.ZvFpa.QUALITY_BITMASK;
            InitialOsrZv.L2_QUALITY_BITMASK = InLfrCwf.ZvFpa.L2_QUALITY_BITMASK;
            InitialOsrZv.DELTA_PLUS_MINUS   = InLfrCwf.Zv.DELTA_PLUS_MINUS;
            %
            InitialOsr = struct(...
                'Ga', InitialGa, ...
                'Zv', InitialOsrZv);
            %
            [InitialDsrZv, iRecordsInBinCa] = bicas.proc.dsr.init_shared_DSR_ZVs(...
                InLfrCwf, ...
                BIN_LENGTH_WOLS_NS, ...
                BIN_TIMESTAMP_POS_WOLS_NS, ...
                L);
            % NOTE: Not setting DSR ".Ga"/global attributes here, since DSR
            % datasets later copy ".Ga" from the respective OSR datasets.
            InitialDsr = struct('Zv', InitialDsrZv);
            % nRecordsDsr = numel(iRecordsInBinCa);
            


            %====================
            % ZVs for EFIELD OSR
            %====================
            OutEfieldOsr = InitialOsr;
            OutEfieldOsr.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutEfieldOsr.Zv.EDC_SRF                   = R.EdcSrfMvpmFpa.cast('single');
            %
            b = all(OutEfieldOsr.Zv.EDC_SRF.fpAr, 2);    % Rows which are only FPs.
            OutEfieldOsr.Zv.QUALITY_FLAG(b) = bicas.utils.FPArray.FP_UINT8;
        


            %===================
            % ZVs for SCPOT OSR
            %===================
            OutScpotOsr = InitialOsr;
            OutScpotOsr.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutScpotOsr.Zv.SCPOT                     = R.ScpotVoltFpa.cast('single');
            OutScpotOsr.Zv.PSP                       = R.PspVoltFpa.cast('single');
            %
            b = OutScpotOsr.Zv.SCPOT.fpAr & ...
                OutScpotOsr.Zv.PSP.fpAr;
            OutScpotOsr.Zv.QUALITY_FLAG(b) = bicas.utils.FPArray.FP_UINT8;



            %=====================
            % ZVs for DENSITY OSR
            %=====================
            OutDensityOsr = InitialOsr;
            OutDensityOsr.Ga.Misc_calibration_versions = gaDensity_Misc_calibration_versions;
            %
            OutDensityOsr.Zv.DENSITY                   = R.NeScpCm3Fpa.cast('single');
            %
            b = OutDensityOsr.Zv.DENSITY.fpAr;
            OutDensityOsr.Zv.QUALITY_FLAG(b)           = bicas.utils.FPArray.FP_UINT8;



            %====================
            % ZVs for EFIELD DSR
            %====================
            OutEfieldDsr    = InitialDsr;
            OutEfieldDsr.Ga = OutEfieldOsr.Ga;
            %            
            [EdcSrfDsrFpa, EdcstdSrfDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                R.EdcSrfMvpmFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            OutEfieldDsr.Zv.EDC_SRF    = EdcSrfDsrFpa.cast('single');
            OutEfieldDsr.Zv.EDCSTD_SRF = EdcstdSrfDsrFpa.cast('single');
            %
            b = all(OutEfieldDsr.Zv.EDC_SRF.fpAr, 2);    % Rows which are only FPs.
            OutEfieldDsr.Zv.QUALITY_FLAG(b) = bicas.utils.FPArray.FP_UINT8;

            

            %===================
            % ZVs for SCPOT DSR
            %===================
            OutScpotDsr    = InitialDsr;
            OutScpotDsr.Ga = OutScpotOsr.Ga;
            %
            [ScpotDsrFpa, ScpotstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                R.ScpotVoltFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            OutScpotDsr.Zv.SCPOT    = ScpotDsrFpa.cast('single');
            OutScpotDsr.Zv.SCPOTSTD = ScpotstdDsrFpa.cast('single');
            %
            [PspDsrFpa, PspstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                R.PspVoltFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            OutScpotDsr.Zv.PSP    = PspDsrFpa.cast('single');
            OutScpotDsr.Zv.PSPSTD = PspstdDsrFpa.cast('single');
            %
            b = OutScpotDsr.Zv.SCPOT.fpAr & ...
                OutScpotDsr.Zv.PSP.fpAr;
            OutScpotDsr.Zv.QUALITY_FLAG(b) = bicas.utils.FPArray.FP_UINT8;



            %=====================
            % ZVs for DENSITY DSR
            %=====================
            OutDensityDsr    = InitialDsr;
            OutDensityDsr.Ga = OutDensityOsr.Ga;
            %
            [DensityDsrFpa, DensitystdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                R.NeScpCm3Fpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            OutDensityDsr.Zv.DENSITY    = DensityDsrFpa.cast('single');
            OutDensityDsr.Zv.DENSITYSTD = DensitystdDsrFpa.cast('single');
            %
            b = OutDensityDsr.Zv.DENSITY.fpAr;
            OutDensityDsr.Zv.QUALITY_FLAG(b) = bicas.utils.FPArray.FP_UINT8;



            nRecordsOsr = size(InLfrCwf.Zv.Epoch,   1);
            nRecordsDsr = size(InitialDsr.Zv.Epoch, 1);
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L3.process_L2_to_L3', tTicToc, ...
                nRecordsOsr, 'OSR record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L3.process_L2_to_L3', tTicToc, ...
                nRecordsDsr, 'DSR record')

        end    % process_L2_to_L3
        
        
        
    end    % methods(Static, Access=public)



end
