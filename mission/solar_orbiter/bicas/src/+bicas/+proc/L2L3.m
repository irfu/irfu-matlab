%
% Class that collects "processing functions" as public static methods. Only
% covers processing L2-->L3.
%
% This class is not meant to be instantiated.
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
                = process_L2_to_L3(InLfrCwf, SETTINGS, L)

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
            % Regular expression for the format of version strings from
            % BICAS-external code.
            % Equivalent to: yyyy-mm-ddThh:mm:ss
            CODE_VER_STR_REGEXP = ...
                '[0-9]{4}-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]';
            


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
            % solo.vdccal()
            R = bicas.proc.L2L3.calc_EFIELD_SCPOT(InLfrCwf.Zv, SETTINGS);
            % solo.psp2ne()
            [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr] = bicas.proc.L2L3.calc_DENSITY(R.PspTs);
            clear NeScpQualityBitTs
            % NOTE: Ignoring return value NeScpQualityBit(Ts) for now. Value is
            %       expected to be used by BICAS later.


            
            %===================================================================
            % ~HACK: MODIFY INPUT ARGUMENT InLfrCwf
            % -------------------------------------
            % IMPLEMENTATION NOTE: This is to modify QUALITY_FLAG for both OSR
            % and DSR datasets. In principle, this is for keeping the interface
            % to bicas.proc.dsr.init_shared_DSR_ZVs() simple.
            %===================================================================
            InLfrCwf.Zv.QUALITY_FLAG(R.bNotUsed) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %====================================================================
            % Derive values for CDF global attribute "Misc_calibration_versions"
            %====================================================================
            assert(~isempty(R.vdccalMatVerStr), ...
                ['solo.vdccal() returns an empty vdccalMatVerStr', ...
                ' (string representing the version of the corresponding', ...
                ' .mat file). BICAS therefore needs to be updated.'])
            irf.assert.castring_regexp(R.vdccalCodeVerStr, CODE_VER_STR_REGEXP)
            irf.assert.castring_regexp(psp2neCodeVerStr,   CODE_VER_STR_REGEXP)
            %
            % NOTE: Should not add BICAS version to glob.attr.
            % "Misc_calibration_versions" since this is already encoded in
            % global attribute "Software_version" (together with
            % "Software_name").
            %
            % NOTE: Density "Misc_calibration_versions" contains all three
            % versions, since density is derived from PSP.
            vdccalStr    = ['solo.vdccal() code version ',     R.vdccalCodeVerStr];
            vdccalMatStr = ['solo.vdccal() calibration file ', R.vdccalMatVerStr];
            psp2neStr    = ['solo.psp2ne() code version ',     psp2neCodeVerStr];
            gaEfieldScpot_Misc_calibration_versions = {vdccalStr, vdccalMatStr};
            gaDensity_Misc_calibration_versions     = ...
                [gaEfieldScpot_Misc_calibration_versions, {psp2neStr}];

            
            
            %=========================================
            % Misc. variables shared between datasets
            %=========================================
            % Global attributes -- shared between all OSR+DSR datasets.
            InitialGa = struct();
            InitialGa.OBS_ID                 = InLfrCwf.Ga.OBS_ID;
            InitialGa.SOOP_TYPE              = InLfrCwf.Ga.SOOP_TYPE;
            % zVariables -- shared between all OSR datasets.
            InitialOsrZv = struct();
            InitialOsrZv.Epoch              = InLfrCwf.Zv.Epoch;
            InitialOsrZv.QUALITY_BITMASK    = InLfrCwf.Zv.QUALITY_BITMASK;
            InitialOsrZv.L2_QUALITY_BITMASK = InLfrCwf.Zv.L2_QUALITY_BITMASK;
            InitialOsrZv.QUALITY_FLAG       = InLfrCwf.Zv.QUALITY_FLAG;
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
            OutEfieldOsr.Zv.EDC_SRF                   = R.zvEdcMvpm;
            %
            b = all(isnan(OutEfieldOsr.Zv.EDC_SRF), 2);
            OutEfieldOsr.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;
        


            %===================
            % ZVs for SCPOT OSR
            %===================
            OutScpotOsr = InitialOsr;
            OutScpotOsr.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutScpotOsr.Zv.SCPOT                     = R.ScpotTs.data;
            OutScpotOsr.Zv.PSP                       = R.PspTs.data;
            %
            b = isnan(OutScpotOsr.Zv.SCPOT) & ...
                isnan(OutScpotOsr.Zv.PSP);
            OutScpotOsr.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %=====================
            % ZVs for DENSITY OSR
            %=====================
            OutDensityOsr = InitialOsr;
            OutDensityOsr.Ga.Misc_calibration_versions = gaDensity_Misc_calibration_versions;
            %
            OutDensityOsr.Zv.DENSITY                   = NeScpTs.data;
            %
            b = isnan(OutDensityOsr.Zv.DENSITY);
            OutDensityOsr.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %====================
            % ZVs for EFIELD DSR
            %====================
            OutEfieldDsr    = InitialDsr;
            OutEfieldDsr.Ga = OutEfieldOsr.Ga;
            %
            [OutEfieldDsr.Zv.EDC_SRF, ...
             OutEfieldDsr.Zv.EDCSTD_SRF] = bicas.proc.dsr.downsample_sci_ZV(...
                OutEfieldOsr.Zv.EDC_SRF, ...
                bicas.const.N_MIN_SAMPLES_PER_DSR_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            % NOTE: Merge across samples in same record.
            b = all(isnan(OutEfieldDsr.Zv.EDC_SRF), 2);
            OutEfieldDsr.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;

            

            %===================
            % ZVs for SCPOT DSR
            %===================
            OutScpotDsr    = InitialDsr;
            OutScpotDsr.Ga = OutScpotOsr.Ga;
            %
            [OutScpotDsr.Zv.SCPOT, ...
             OutScpotDsr.Zv.SCPOTSTD] = bicas.proc.dsr.downsample_sci_ZV(...
                OutScpotOsr.Zv.SCPOT, ...
                bicas.const.N_MIN_SAMPLES_PER_DSR_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            [OutScpotDsr.Zv.PSP, ...
             OutScpotDsr.Zv.PSPSTD] = bicas.proc.dsr.downsample_sci_ZV(...
                OutScpotOsr.Zv.PSP, ...
                bicas.const.N_MIN_SAMPLES_PER_DSR_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            b = isnan(OutScpotDsr.Zv.SCPOT) & ...
                isnan(OutScpotDsr.Zv.PSP);
            OutScpotDsr.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %=====================
            % ZVs for DENSITY DSR
            %=====================
            OutDensityDsr    = InitialDsr;
            OutDensityDsr.Ga = OutDensityOsr.Ga;
            %
            [OutDensityDsr.Zv.DENSITY, ...
             OutDensityDsr.Zv.DENSITYSTD] = bicas.proc.dsr.downsample_sci_ZV(...
                OutDensityOsr.Zv.DENSITY, ...
                bicas.const.N_MIN_SAMPLES_PER_DSR_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            b = isnan(OutDensityDsr.Zv.DENSITY);
            OutDensityDsr.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



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



    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################
    
    
    
        % Calculate both
        %   (1) ELECTRIC FIELD, and
        %   (2) SPACECRAFT POTENTIALS
        % via the same BICAS-external code solo.vdccal() (still inside
        % irfu-matlab).
        %
        % Largely a wrapper around solo.vdccal().
        %
        % NOTE: Needs to be careful with the units, and incompatible updates to
        % solo.vdccal() without the knowledge of the BICAS author. Therefore
        % uses extra assertions to detect such changes.
        %
        % RETURN VALUE
        % ============
        % R : Struct with multiple variables.
        %     NOTE: Return values are packaged as a struct to provide named
        %     return values and avoid confusing similar return results with each
        %     other.
        %
        function R = calc_EFIELD_SCPOT(...
                InLfrCwfZv, SETTINGS)
            
            QUALITY_FLAG_minForUse = SETTINGS.get_fv(...
                'PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN');


            
            % Shorten recurring variables.
            zv_VDC   = InLfrCwfZv.VDC;
            zv_EDC   = InLfrCwfZv.EDC;
            zv_Epoch = InLfrCwfZv.Epoch;            
            
            %======================================================
            % Create input variables for solo.vdccal()
            % ----------------------------------------
            % Set records to NaN for QUALITY_FLAG below threshold.
            %======================================================
            % NOTE: Comparison will technically fail for QUALITY_FLAG fill
            % value, but that is acceptable (ideal result is ambiguous anyway).
            bNotUsed = InLfrCwfZv.QUALITY_FLAG < QUALITY_FLAG_minForUse;
            clear InLfrCwfZv
            zv_VDC(bNotUsed, :) = NaN;
            zv_EDC(bNotUsed, :) = NaN;
            %
            % NOTE: Should TSeries objects really use TensorOrder=1 and
            % repres={x,y,z}?!! VDC and EDC are not time series of vectors, but
            % fo three scalars. Probably does not matter. solo.vdccal() does
            % indeed use VDC.x, EDC.x etc.
            VdcTs = TSeries(...
                EpochTT(zv_Epoch), zv_VDC, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});
            EdcTs = TSeries(...
                EpochTT(zv_Epoch), zv_EDC, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});
            
            
            
            %==========================
            % CALL BICAS-EXTERNAL CODE
            %==========================
            % NOTE: Not specifying calibration file.
            % ==> Use current official calibration file, hardcoded in
            %     solo.vdccal(), that should be used for official datasets.
            [EdcSrfTs, PspTs, ScpotTs, vdccalCodeVerStr, vdccalMatVerStr] ...
                = solo.vdccal(VdcTs, EdcTs, []);
            clear VdcTs EdcTs
            %==========================
            
            
            
            % ASSERTIONS: Check solo.vdccal() return values.
            irf.assert.sizes(...
                zv_Epoch,      [-1, 1], ...
                EdcSrfTs.data, [-1, 3], ...
                PspTs.data,    [-1, 1], ...
                ScpotTs.data,  [-1, 1]);
            assert(strcmp(EdcSrfTs.units,            'mV/m'))
            assert(strcmp(EdcSrfTs.coordinateSystem, 'SRF'))
            assert(strcmp(PspTs.units,               'V'))
            assert(strcmp(ScpotTs.units,             'V'))
            assert(all(zv_Epoch == EdcSrfTs.time.ttns))
            assert(all(zv_Epoch ==    PspTs.time.ttns))
            assert(all(zv_Epoch ==  ScpotTs.time.ttns))

            
            
            %===================================================================
            % Normalize the representation of E-field X-component
            % (EdcSrfTs --> zvEdcMvpm)
            % ---------------------------------------------------
            % Set E_x = NaN, but ONLY if assertion deems that the corresponding
            % information is missing.
            %
            % IMPLEMENTATION NOTE: solo.vdccal() sets antenna 1 values to be
            % zero, if its input data is non-fill value/NaN, but NaN if fill
            % value. Must therefore check for both zero and NaN.
            % Ex: Dataset 2020-08-01
            %===================================================================
            zvEdcMvpm = EdcSrfTs.data;    % MVPM = mV/m
            clear EdcSrfTs
            % IMPLEMENTATION NOTE: ismember() does not work for NaN.
            assert(all(zvEdcMvpm(:, 1) == 0 | isnan(zvEdcMvpm(:, 1))), ...
                ['EDC for antenna 1 returned from', ...
                ' solo.vdccal() is neither zero nor NaN and can therefore', ...
                ' not be assumed to be unknown anymore.', ...
                ' Verify that this is correct solo.vdccal() behaviour and', ...
                ' (if correct) then update BICAS to handle this.'])
            zvEdcMvpm(:, 1) = NaN;
            
            
            
            % Prepare return struct.
            R = [];
            R.PspTs            = PspTs;
            R.ScpotTs          = ScpotTs;
            R.zvEdcMvpm        = zvEdcMvpm;
            R.vdccalCodeVerStr = vdccalCodeVerStr;
            R.vdccalMatVerStr  = vdccalMatVerStr;
            R.bNotUsed         = bNotUsed;
            
        end



        % Calculate DENSITY via a BICAS-external code solo.psp2ne() (still
        % inside irfu-matlab).
        %
        % Essentially a wrapper around solo.psp2ne().
        % 
        % NOTE: One needs to be careful with units and incompatible updates to
        % solo.vdccal() without the knowledge of the BICAS author. Therefore
        % uses extra assertions to detect such changes.
        %
        % NOTE: Empirically, some return values are NaN.
        % NOTE: Shortening "SCP" = SCPOT comes from the return variable name in
        % solo.psp2ne().
        %
        % IMPLEMENTATION NOTE: Does not need to check QUALITY_FLAG limit since
        % relies on PSP values for which this has already been done.
        %
        function [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr] = calc_DENSITY(PspTs)
            
            %==========================
            % CALL BICAS-EXTERNAL CODE
            %==========================
            [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr] = solo.psp2ne(PspTs);
            %==========================
            


            % ASSERTIONS: Check solo.psp2ne() return values.
            irf.assert.sizes(...
                PspTs.data,             [-1, 1], ...
                NeScpTs.data,           [-1, 1], ...
                NeScpQualityBitTs.data, [-1, 1] ...
            );
            assert(strcmp(NeScpTs.units, 'cm^-3'))
            assert(all(PspTs.time == NeScpTs.time          ))
            assert(all(PspTs.time == NeScpQualityBitTs.time))
            assert(all( (NeScpTs.data > 0) | isnan(NeScpTs.data)), ...
                'solo.psp2ne() returned non-positive (non-NaN) plasma density.')
            assert(strcmp(NeScpTs.units, 'cm^-3'))
            % NOTE: Not permitting NaN quality bit. Unsure if that is the
            %       best behaviour.
            assert(...
                all(ismember(NeScpQualityBitTs.data, [0, 1])), ...
                'solo.psp2ne() returned illegal NeScpTsQualityBitTs.')
        end



    end    % methods(Static, Access=private)

end
