%
% SWMP for processing L2 LFR CWF to L3 OSR+DSR density, E field, and ScPot.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like"
%   data, use the first MATLAB array index to represent CDF records.
%
%
% BUG?
% ====
% bNotUsed is used when calling
% bicas.proc.L2L3.ext.calc_EFIELD_SCPOT() but not when calling
% bicas.proc.L2L3.ext.calc_DENSITY() (should maybe be used for both).
% At the same time, bNotUsed is used for creating the OSR template which is used
% for all datasets (DENSITY, EFIELD, SCPOT).
%   InLfrCwf.ZvFpa.QUALITY_FLAG(R.bNotUsed) = bicas.utils.FPArray.FP_UINT8;
%   TemplateOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.get_OSR_template(InLfrCwf);
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef L3OsrDsrSwmProcessing < bicas.proc.SwmProcessing
% PROPOSAL: Automatic test code.
%   NOTE: There are limited tests.
%
% PROPOSAL: Better name.
%   OSR, DSR
%   L3
%   Density, Efield, ScPot = DES
%
% PROPOSAL: Split up processing between (a) density, and (b) E-field & SCPOT.
%   PRO: Faster
%       CON: Not very heavy operation.
%   PRO: Leads to better organization of code.
%       PRO: process_L2_to_L3() is too large and should be split up anyway.
%   CON: DENSITY is a function EFIELD+SCPOT, and thus has to be processed after
%        the latter.
%
% PROPOSAL: Instead of sharing initial "template variables", have function for
%           generating those variables.



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        % OVERRIDE
        function OutputDatasetsMap = production_function(obj, ...
            InputDatasetsMap, rctDir, NsoTable, Bso, L)
            
            InputLfrCwfCdf = InputDatasetsMap('LFR-SURV-CWF-E_cdf');

            Ec = bicas.proc.L2L3.ExternalCodeImplementation();

            %==============
            % Process data
            %==============
            [EfieldOsrCdf,  EfieldDsrCdf, ...
             ScpotOsrCdf,   ScpotDsrCdf, ...
             DensityOsrCdf, DensityDsrCdf] = ...
                bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3(InputLfrCwfCdf, Ec, Bso, L);

            OutputDatasetsMap = containers.Map();
            OutputDatasetsMap('EFIELD_OSR_cdf')  = EfieldOsrCdf;
            OutputDatasetsMap('EFIELD_DSR_cdf')  = EfieldDsrCdf;
            OutputDatasetsMap('SCPOT_OSR_cdf')   = ScpotOsrCdf;
            OutputDatasetsMap('SCPOT_DSR_cdf')   = ScpotDsrCdf;
            OutputDatasetsMap('DENSITY_OSR_cdf') = DensityOsrCdf;
            OutputDatasetsMap('DENSITY_DSR_cdf') = DensityDsrCdf;
        end



    end    % methods(Access=public)



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % ~Process L2-->L3 (not VHT).
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
                = process_L2_to_L3(InLfrCwf, Ec, Bso, L)

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

            Tmk = bicas.utils.Timekeeper('bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3', L);
            assert(isa(Ec, 'bicas.proc.L2L3.ExternalCodeAbstract'))



            %===========
            % Constants
            %===========
            % The only acceptable input DSI.
            %INPUT_DSI                 = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
            % NS = Nanoseconds
            BIN_LENGTH_WOLS_NS        = int64(10e9);
            BIN_TIMESTAMP_POS_WOLS_NS = int64(BIN_LENGTH_WOLS_NS / 2);
            


            %=======================================
            % Call BICAS-external code to calculate
            % (1) EFIELD, SCPOT, and from that
            % (2) DENSITY.
            %=======================================
            LfrCwfZv = [];
            LfrCwfZv.Epoch            = InLfrCwf.Zv.Epoch;
            LfrCwfZv.VDC_Fpa          = InLfrCwf.ZvFpa.VDC;
            LfrCwfZv.EDC_Fpa          = InLfrCwf.ZvFpa.EDC;
            LfrCwfZv.QUALITY_FLAG_Fpa = InLfrCwf.ZvFpa.QUALITY_FLAG;
            R = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT_DENSITY(LfrCwfZv, Ec, Bso);



            %===================================================================
            % ~HACK: MODIFY INPUT ARGUMENT InLfrCwf
            % -------------------------------------
            % IMPLEMENTATION NOTE: This is to modify QUALITY_FLAG for both OSR
            % and DSR datasets. In principle, this is for keeping the interface
            % to bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template() simple.
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

            
            
            %================================================================
            % Misc. variables shared between datasets and later modified for
            % specific datasets
            %================================================================
            TemplateOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.get_OSR_template(InLfrCwf);
            [TemplateDsrZv, iRecordsInBinCa] = bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template(...
                InLfrCwf, ...
                BIN_LENGTH_WOLS_NS, ...
                BIN_TIMESTAMP_POS_WOLS_NS, ...
                L);
            % NOTE: Not setting DSR ".Ga"/global attributes here, since DSR
            % datasets later copy ".Ga" from the respective OSR datasets.
            TemplateDsr = struct('Zv', TemplateDsrZv);



            %=======================================
            % Generate data structures for datasets
            %=======================================            
            OutEfieldOsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_efield( TemplateOsr, R.EdcSrfMvpmFpa,                       gaEfieldScpot_Misc_calibration_versions);
            OutScpotOsr   = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_scpot(  TemplateOsr, R.ScpotVoltFpa,  R.PspVoltFpa,         gaEfieldScpot_Misc_calibration_versions);
            OutDensityOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_density(TemplateOsr, R.NeScpCm3Fpa,   R.NeScpQualityBitFpa, gaDensity_Misc_calibration_versions);
            %
            OutEfieldDsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_efield( TemplateDsr, OutEfieldOsr.Ga,  R.EdcSrfMvpmFpa,                                    iRecordsInBinCa, L);
            OutScpotDsr   = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_scpot(  TemplateDsr, OutScpotOsr.Ga,   R.ScpotVoltFpa, R.PspVoltFpa,                       iRecordsInBinCa, L);
            OutDensityDsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_density(TemplateDsr, OutDensityOsr.Ga, R.NeScpCm3Fpa, OutDensityOsr.Zv.L3_QUALITY_BITMASK, iRecordsInBinCa, L);



            nRecordsOsr = size(InLfrCwf.Zv.Epoch,    1);
            nRecordsDsr = size(TemplateDsr.Zv.Epoch, 1);
            Tmk.stop_log(nRecordsOsr, 'OSR record', nRecordsDsr, 'DSR record')
        end    % process_L2_to_L3



    end    % methods(Static)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Starting template for OSR datasets. Return value is modified 
        function TemplateOsr = get_OSR_template(InLfrCwf)
            Ga = struct();
            Ga.OBS_ID             = InLfrCwf.Ga.OBS_ID;
            Ga.SOOP_TYPE          = InLfrCwf.Ga.SOOP_TYPE;

            Zv = struct();
            Zv.Epoch              = InLfrCwf.Zv.Epoch;
            Zv.QUALITY_FLAG       = InLfrCwf.ZvFpa.QUALITY_FLAG;
            Zv.QUALITY_BITMASK    = InLfrCwf.ZvFpa.QUALITY_BITMASK;
            Zv.L2_QUALITY_BITMASK = InLfrCwf.ZvFpa.L2_QUALITY_BITMASK;
            Zv.DELTA_PLUS_MINUS   = InLfrCwf.ZvFpa.DELTA_PLUS_MINUS;

            TemplateOsr = struct('Ga', Ga, 'Zv', Zv);
        end
        
        
        
        function Out = OSR_efield(TemplateOsr, EdcSrfMvpmFpa, gaMisc_calibration_versions)
            Out = TemplateOsr;
            Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;
            
            Out.Zv.EDC_SRF                   = EdcSrfMvpmFpa.cast('single');
            
            bFp = all(Out.Zv.EDC_SRF.fpAr, 2);    % Rows which are only FPs.
            Out.Zv.QUALITY_FLAG(bFp)         = bicas.utils.FPArray.FP_UINT8;
        end



        function Out = OSR_scpot(TemplateOsr, ScpotVoltFpa, PspVoltFpa, gaMisc_calibration_versions)
            Out = TemplateOsr;
            Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;

            Out.Zv.SCPOT                     = ScpotVoltFpa.cast('single');
            Out.Zv.PSP                       = PspVoltFpa.  cast('single');

            bFp = Out.Zv.SCPOT.fpAr & ...
                  Out.Zv.PSP.fpAr;
            Out.Zv.QUALITY_FLAG(bFp)         = bicas.utils.FPArray.FP_UINT8;
        end



        function Out = OSR_density(TemplateOsr, NeScpCm3Fpa, NeScpQualityBitFpa, gaMisc_calibration_versions)
            Out = TemplateOsr;
            Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;

            Out.Zv.DENSITY                   = NeScpCm3Fpa.cast('single');

            % NOTE: Behaviour w.r.t. FPs:
            %   Density FP     ==> L3_QUALITY_BITMASK FP
            %                      QUALITY_FLAG       FP
            %   Density bit FP ==> L3_QUALITY_BITMASK density bit=false
            %                      (since there is no FP for individual quality bits).
            [QUALITY_FLAG, L3_QUALITY_FLAG] = bicas.proc.L2L3.qual.get_quality_ZVs_density(NeScpQualityBitFpa.array(false));
            Out.Zv.QUALITY_FLAG              = Out.Zv.QUALITY_FLAG.min(QUALITY_FLAG);
            Out.Zv.L3_QUALITY_BITMASK        = bicas.utils.FPArray(L3_QUALITY_FLAG);
            
            bFp = Out.Zv.DENSITY.fpAr;
            Out.Zv.QUALITY_FLAG(bFp)         = bicas.utils.FPArray.FP_UINT8;
            Out.Zv.L3_QUALITY_BITMASK(bFp)   = bicas.utils.FPArray.FP_UINT16;
        end



        function Out = DSR_efield(TemplateDsr, OutEfieldOsrGa, EdcSrfMvpmOsrFpa, iRecordsInBinCa, L)
            Out    = TemplateDsr;
            Out.Ga = OutEfieldOsrGa;

            [EdcSrfDsrFpa, EdcstdSrfDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                EdcSrfMvpmOsrFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, L);
            Out.Zv.EDC_SRF    = EdcSrfDsrFpa.   cast('single');
            Out.Zv.EDCSTD_SRF = EdcstdSrfDsrFpa.cast('single');

            bFp = all(Out.Zv.EDC_SRF.fpAr, 2);    % Rows which are only FPs.
            Out.Zv.QUALITY_FLAG(bFp) = bicas.utils.FPArray.FP_UINT8;
        end



        function Out = DSR_scpot(TemplateDsr, OutScpotOsrGa, ScpotVoltOsrFpa, PspVoltOsrFpa, iRecordsInBinCa, L)
            Out    = TemplateDsr;
            Out.Ga = OutScpotOsrGa;

            % Downsample SCPOT
            [ScpotDsrFpa, ScpotstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                ScpotVoltOsrFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, ...
                L);
            Out.Zv.SCPOT    = ScpotDsrFpa.   cast('single');
            Out.Zv.SCPOTSTD = ScpotstdDsrFpa.cast('single');

            % Downsample PSP
            [PspDsrFpa, PspstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                PspVoltOsrFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, L);
            Out.Zv.PSP    = PspDsrFpa.   cast('single');
            Out.Zv.PSPSTD = PspstdDsrFpa.cast('single');

            bFp = Out.Zv.SCPOT.fpAr & ...
                  Out.Zv.PSP.fpAr;
            Out.Zv.QUALITY_FLAG(bFp) = bicas.utils.FPArray.FP_UINT8;
        end



        function Out = DSR_density(TemplateDsr, OutDensityOsrGa, NeScpCm3OsrFpa, osr_L3_QUALITY_BITMASK, iRecordsInBinCa, L)
            Out    = TemplateDsr;
            Out.Ga = OutDensityOsrGa;

            [DensityDsrFpa, DensitystdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
                NeScpCm3OsrFpa, ...
                bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
                iRecordsInBinCa, L);
            Out.Zv.DENSITY    = DensityDsrFpa.   cast('single');
            Out.Zv.DENSITYSTD = DensitystdDsrFpa.cast('single');

            bFp = Out.Zv.DENSITY.fpAr;
            Out.Zv.QUALITY_FLAG(bFp)       = bicas.utils.FPArray.FP_UINT8;
            Out.Zv.L3_QUALITY_BITMASK      = bicas.proc.dsr.downsample_ZV_bitmask(...
                osr_L3_QUALITY_BITMASK, iRecordsInBinCa);
            Out.Zv.L3_QUALITY_BITMASK(bFp) = bicas.utils.FPArray.FP_UINT16;   % ?!
        end



    end    % methods(Static, Access=private)



end
