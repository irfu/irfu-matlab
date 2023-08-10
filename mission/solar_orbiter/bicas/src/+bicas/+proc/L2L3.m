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
        % Both ORIS and DWNS. The same is not(?) enforced in L2 processing, but
        % should maybe be. /EJ 2021-05-12
        %
        function [OutEfieldOris,  OutEfieldDwns, ...
                  OutScpotOris,   OutScpotDwns, ...
                  OutDensityOris, OutDensityDwns] ...
                = process_L2_to_L3(InLfrCwf, SETTINGS, L)

            % PROPOSAL: Split up into different parts for EFIELD, SCPOT, DENSITY
            %           (still combine non-downsampled and downsampled).
            %   CON: Slows down overall processing.
            %       PRO: Must read same L2 dataset multiple times.
            %       PRO: Must read L3 SCPOT dataset to produce L3 DENSITY dataset.
            %   CON: There is much shared functionality for 3 quality zVars.
            %       PRO: Same ~constants
            %           Ex: INPUT_DATASET_ID, BIN_LENGTH_WOLS_NS, BIN_TIMESTAMP_POS_WOLS_NS
            %       PRO: Read setting QUALITY_FLAG_MIN_FOR_USE
            %       PRO: Normalizing CWF zVar names.
            %       PRO: Preparations for downsampled.
            %           Bin locations, bundling of records,
            %           Downsampling of quality variables
            %               (QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK).
            %           DELTA_PLUS_MINUS_dwns
            %
            % NOTE: ROC BUG:
            %               https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/48
            %         L1 QUALITY_BITMASK seems to use the wrong value (255) as
            %         fill value (FILLVAL=65535). ==> A bug fix would not fix
            %         the entire issue.
            %   PROPOSAL: Use double also for CDF integer variables so NaN can
            %             represent fill value also for these.
            %
            % NOTE: L2 LFR-CWF-E skt previously had zVar
            %   QUALITY_BITMASK=CDF_UINT1, fill value=255 (wrong)
            % until skt V12 when it was changed to
            %   QUALITY_BITMASK=CDF_UINT2, fill value 65535 (correct).

            tTicToc = tic();



            %===========
            % Constants
            %===========
            % The only acceptable input DATASET_ID.
            INPUT_DATASET_ID          = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
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

            bicas.proc.utils.handle_zv_name_change(...
                fnChangeList, INPUT_DATASET_ID, SETTINGS, L, 'VDC', ...
                'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')
            
            

            %=======================================
            % Call BICAS-external code to calculate
            % (1) EFIELD, SCPOT, and from that
            % (2) DENSITY.
            %=======================================
            % solo.vdccal()
            R = bicas.proc.L2L3.calc_EFIELD_SCPOT(InLfrCwf.Zv, SETTINGS);
            % solo.psp2ne()
            [NeScpTs, psp2neCodeVerStr] = bicas.proc.L2L3.calc_DENSITY(R.PspTs);


            
            %===================================================================
            % ~HACK: MODIFY INPUT ARGUMENT InLfrCwf
            % -------------------------------------
            % IMPLEMENTATION NOTE: This is to modify QUALITY_FLAG for both ORIS
            % and DWNS datasets. In principle, this is to keep the interface to
            % init_shared_downsampled() simple.
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
            % Global attributes -- shared between all ORIS+DWNS datasets.
            InitialGa = struct();
            InitialGa.OBS_ID                 = InLfrCwf.Ga.OBS_ID;
            InitialGa.SOOP_TYPE              = InLfrCwf.Ga.SOOP_TYPE;
            % zVariables -- shared between all ORIS datasets.
            InitialOrisZv = struct();
            InitialOrisZv.Epoch              = InLfrCwf.Zv.Epoch;
            InitialOrisZv.QUALITY_BITMASK    = InLfrCwf.Zv.QUALITY_BITMASK;
            InitialOrisZv.L2_QUALITY_BITMASK = InLfrCwf.Zv.L2_QUALITY_BITMASK;
            InitialOrisZv.QUALITY_FLAG       = InLfrCwf.Zv.QUALITY_FLAG;
            InitialOrisZv.DELTA_PLUS_MINUS   = InLfrCwf.Zv.DELTA_PLUS_MINUS;
            %
            InitialOris = struct(...
                'Ga', InitialGa, ...
                'Zv', InitialOrisZv);
            %
            [InitialDwnsZv, iRecordsInBinCa] = bicas.proc.dwns.init_shared_downsampled(...
                InLfrCwf, ...
                BIN_LENGTH_WOLS_NS, ...
                BIN_TIMESTAMP_POS_WOLS_NS, ...
                L);
            % NOTE: Not setting DWNS ".Ga"/global attributes here, since DWNS
            % datasets later copy ".Ga" from the respective ORIS datasets.
            InitialDwns = struct('Zv', InitialDwnsZv);
            


            %=======================
            % zVars for EFIELD ORIS
            %=======================
            OutEfieldOris = InitialOris;
            OutEfieldOris.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutEfieldOris.Zv.EDC_SRF                   = R.zvEdcMvpm;
            %
            b = all(isnan(OutEfieldOris.Zv.EDC_SRF), 2);
            OutEfieldOris.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;
        


            %======================
            % zVars for SCPOT ORIS
            %======================
            OutScpotOris = InitialOris;
            OutScpotOris.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            %
            OutScpotOris.Zv.SCPOT                     = R.ScpotTs.data;
            OutScpotOris.Zv.PSP                       = R.PspTs.data;
            %
            b = isnan(OutScpotOris.Zv.SCPOT) & ...
                isnan(OutScpotOris.Zv.PSP);
            OutScpotOris.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %========================
            % zVars for DENSITY ORIS
            %========================
            OutDensityOris = InitialOris;
            OutDensityOris.Ga.Misc_calibration_versions = gaDensity_Misc_calibration_versions;
            %
            OutDensityOris.Zv.DENSITY                   = NeScpTs.data;
            %
            b = isnan(OutDensityOris.Zv.DENSITY);
            OutDensityOris.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %==============================
            % zVars for EFIELD DOWNSAMPLED
            %==============================
            OutEfieldDwns    = InitialDwns;
            OutEfieldDwns.Ga = OutEfieldOris.Ga;
            %
            [OutEfieldDwns.Zv.EDC_SRF, ...
             OutEfieldDwns.Zv.EDCSTD_SRF] = bicas.proc.dwns.downsample_sci_zVar(...
                OutEfieldOris.Zv.EDC_SRF, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            % NOTE: Merge across samples in same record.
            b = all(isnan(OutEfieldDwns.Zv.EDC_SRF), 2);
            OutEfieldDwns.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;

            

            %=============================
            % zVars for SCPOT DOWNSAMPLED
            %=============================
            OutScpotDwns    = InitialDwns;
            OutScpotDwns.Ga = OutScpotOris.Ga;
            %
            [OutScpotDwns.Zv.SCPOT, ...
             OutScpotDwns.Zv.SCPOTSTD] = bicas.proc.dwns.downsample_sci_zVar(...
                OutScpotOris.Zv.SCPOT, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            [OutScpotDwns.Zv.PSP, ...
             OutScpotDwns.Zv.PSPSTD] = bicas.proc.dwns.downsample_sci_zVar(...
                OutScpotOris.Zv.PSP, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            b = isnan(OutScpotDwns.Zv.SCPOT) & ...
                isnan(OutScpotDwns.Zv.PSP);
            OutScpotDwns.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            %===============================
            % zVars for DENSITY DOWNSAMPLED
            %===============================
            OutDensityDwns    = InitialDwns;
            OutDensityDwns.Ga = OutDensityOris.Ga;
            %
            [OutDensityDwns.Zv.DENSITY, ...
             OutDensityDwns.Zv.DENSITYSTD] = bicas.proc.dwns.downsample_sci_zVar(...
                OutDensityOris.Zv.DENSITY, ...
                bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN, ...
                iRecordsInBinCa, ...
                L);
            %
            b = isnan(OutDensityDwns.Zv.DENSITY);
            OutDensityDwns.Zv.QUALITY_FLAG(b) = InLfrCwf.ZvFv.QUALITY_FLAG;



            nRecordsOris = size(InLfrCwf.Zv.Epoch,    1);
            nRecordsDwns = size(InitialDwns.Zv.Epoch, 1);
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L3.process_L2_to_L3', tTicToc, ...
                nRecordsOris, 'ORIS record')
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L2L3.process_L2_to_L3', tTicToc, ...
                nRecordsDwns, 'DWNS record')

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
            
            
            
            % ASSERTIONS
            irf.assert.sizes(...
                zv_Epoch,      [-1, 1], ...
                EdcSrfTs.data, [-1, 3], ...
                PspTs.data,    [-1, 1], ...
                ScpotTs.data,  [-1, 1]);
            assert(strcmp(EdcSrfTs.units,            'mV/m'))
            assert(strcmp(EdcSrfTs.coordinateSystem, 'SRF'))
            assert(strcmp(PspTs.units,               'V'))
            assert(strcmp(ScpotTs.units,             'V'))
            
            
            
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
        % NOTE: Shortening "SCP" comes from the return variable name in
        % solo.psp2ne(). Do not know what it means (SpaceCraft Potential?).
        %
        % IMPLEMENTATION NOTE: Does not need to check QUALITY_FLAG limit since
        % relies on PSP values for which this has already been done.
        %
        function [NeScpTs, psp2neCodeVerStr] = calc_DENSITY(PspTs)
            
            %==========================
            % CALL BICAS-EXTERNAL CODE
            %==========================
            [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr] = solo.psp2ne(PspTs);
            %==========================
            % NOTE: Ignoring return value NeScpQualityBit(Ts) for now except for
            %       assertions on it. Value is expected to be used by BICAS
            %       later.
            
            % ASSERTIONS
            irf.assert.sizes(...
                PspTs.data,             [-1, 1], ...
                NeScpTs.data,           [-1, 1], ...
                NeScpQualityBitTs.data, [-1, 1] ...
            );
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
