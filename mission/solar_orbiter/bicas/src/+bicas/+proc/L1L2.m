%
% Class that collects "processing functions" as public static methods. Only
% covers processing L1/L1R-->L2.
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
% SOME INTERMEDIATE PROCESSING DATA FORMATS
% =========================================
% - PreDC = Pre-(Demuxing & Calibration) Data
%       Generic data format that can represent all forms of input datasets
%       before demuxing and calibration. Can use an arbitrary number of samples
%       per record. Some variables are therefore not used in CWF output
%       datasets.
% - PostDC = Post-(Demuxing & Calibration) Data
%       Data format that includes calibrated currents & calibrated & demuxed
%       voltages.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
classdef L1L2
%#######################################################################################################################
%
% PROPOSAL: Move normalize_CALIBRATION_TABLE_INDEX() to some collection of utils.
%   PROPOSAL: bicas.proc.utils
%       CON: Function is too specific. Has inputDsi as argument.
%           CON: Could be less bad than this file.
%
% PROPOSAL: Submit zVar variable attributes.
%   PRO: Can interpret fill values.
%       Ex: Can doublecheck TDS RSWF snapshot length using fill values and compare with zVar SAMPS_PER_CH (which seems
%           to be bad).
%
% PROPOSAL: Return (to execute_SWM), global attributes.
%   PRO: Needed for output datasets: CALIBRATION_TABLE, CALIBRATION_VERSION
%       ~CON: CALIBRATION_VERSION refers to algorithm and should maybe be a SETTING.
%
% PROPOSAL: Class for HkSciTime.
%
%#######################################################################################################################



    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################



        % Processing function
        %
        % NOTE: Only converts relevant HK ZVs to be on SCI Epoch. Later
        % (other) code decides whether to use it (mux mode).
        function HkSciTime = process_HK_CDF_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS, L)
            % PROPOSAL: Separate function for the actual interpolation of data
            %           (changing time array HK-->SCI).
            
            

            % ASSERTIONS
            assert(isa(InSci, 'bicas.InputDataset'))
            assert(isa(InHk,  'bicas.InputDataset'))

            HkSciTime = [];



            %===================================================================
            % Select whether HK should use
            %   (1) (HK) Epoch, or
            %   (2) (HK) ACQUISITION_TIME (not always available).
            % ---------------------------------------------------
            % IMPLEMENTATION NOTE: Historically, there have been datasets where
            % Epoch contains errors, but ACQUISITION_TIME seems OK. This should
            % be phased out eventually.
            %===================================================================
            ACQUISITION_TIME_EPOCH_UTC = SETTINGS.get_fv('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC');
            USE_ZV_ACQUISITION_TIME_HK = SETTINGS.get_fv('PROCESSING.HK.USE_ZV_ACQUISITION_TIME');
            if USE_ZV_ACQUISITION_TIME_HK
                hkEpoch = bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
                    InHk.Zv.ACQUISITION_TIME, ...
                    ACQUISITION_TIME_EPOCH_UTC);

                L.logf('warning', 'Using HK zVar ACQUISITION_TIME instead of Epoch.')
            else
                hkEpoch = InHk.Zv.Epoch;
            end



            %==================================================================
            % Log time intervals to enable comparing available SCI and HK data
            %==================================================================
            TimeVars = [];    % Temporary struct only used for logging.
            TimeVars.HK_Epoch  = InHk.Zv.Epoch;
            TimeVars.SCI_Epoch = InSci.Zv.Epoch;
            if isfield(InHk.Zv, 'ACQUISITION_TIME')
                TimeVars.HK_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
                        InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            if isfield(InSci.Zv, 'ACQUISITION_TIME') && ~isempty(InSci.Zv.ACQUISITION_TIME)
                TimeVars.SCI_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
                    InSci.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            bicas.utils.log_ZVs(TimeVars, SETTINGS, L);



            %===================
            % WARNINGS / ERRORS
            %===================
            if ~issorted(hkEpoch, 'strictascend')
                % Ex: zVar ACQUISITION_TIME in test file
                % TDS___TESTDATA_RGTS_TDS_CALBA_V0.8.6/
                % solo_HK_rpw-bia_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                % is not monotonically increasing (in fact, it is completely
                % strange).
                error(...
                    ['HK timestamps do not increase monotonically', ...
                    ' (USE_ZV_ACQUISITION_TIME_HK=%g).'], ...
                    USE_ZV_ACQUISITION_TIME_HK)
                
            end
            if ~irf.utils.ranges_intersect(InSci.Zv.Epoch, hkEpoch)
                %---------------------------------------
                % CASE: SCI does not overlap HK in time
                %---------------------------------------

                % NOTE: "WARNING" (rather than error) only makes sense if it is
                % possible to later meaningfully permit non-intersection.
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY');
                bicas.default_anomaly_handling(L, ...
                    settingValue, settingKey, 'E+W+illegal', ...
                    'SCI and HK time ranges do not overlap in time.', ...
                    'BICAS:SWMProcessing')
                
            elseif ~irf.utils.is_range_subset(InSci.Zv.Epoch, hkEpoch)
                %-------------------------------------------------
                % CASE: SCI does not cover a subset of HK in time
                %-------------------------------------------------
                % NOTE: This anomaly is obviously implied by the anomaly above
                % (SCI, HK do not overlap). It is therefore only meaningful to
                % detect it if the above anomaly is not detected.
                hk1RelativeSec = 1e-9 * (min(hkEpoch) - min(InSci.Zv.Epoch));
                hk2RelativeSec = 1e-9 * (max(hkEpoch) - max(InSci.Zv.Epoch));

                anomalyDescrMsg = sprintf(...
                    ['HK time range is not a superset of SCI time range.', ...
                    ' Can not reliably interpolate HK data for all of SCI.', ...
                    ' HK begins %g s AFTER SCI begins. HK ends %g s BEFORE SCI ends.'], ...
                    hk1RelativeSec, ...
                    -hk2RelativeSec);

                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.HK.TIME_NOT_SUPERSET_OF_SCI_POLICY');
                bicas.default_anomaly_handling(L, ...
                    settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:DatasetFormat:SWMProcessing')
                
            end



            % Derive time margin within which the nearest HK value will be used.
            % NOTE: Requires >=2 records. 0 or 1 records ==> NaN (and MATLAB
            %       warning).
            %
            % BUG?: This looks as if it could yield bad results. If the margin
            %       is small, then interpolation would yield too many NaN.
            %   NOTE: HK nominally has a much lower "sampling rate" than science
            %         data. ==> The derived margin tends to be higher.
            %   NOTE: Time difference between HK samples varies, which means
            %         that mode() will not necessarily identify the de facto
            %         most common time difference.
            hkEpochExtrapMargin = mode(diff(hkEpoch)) / 2;
            
            

            %=============================================================
            % Derive MUX_SET
            % --------------
            % NOTE: Only obtains one MUX_SET per record
            %       ==> Can not change MUX_SET in the middle of a record.
            % NOTE: Can potentially obtain MUX_SET from LFR SCI.
            %=============================================================
            HkSciTime.MUX_SET = bicas.utils.interpolate_nearest(...
                hkEpochExtrapMargin, ...
                hkEpoch, ...
                InHk.Zv.HK_BIA_MODE_MUX_SET, ...
                InSci.Zv.Epoch);



            %==================================================================
            % Derive DIFF_GAIN
            % ----------------
            % NOTE: Not perfect handling of time when 1 snapshot/record, since
            % one should ideally use time stamps for every LFR _sample_.
            %==================================================================
            HkSciTime.DIFF_GAIN = bicas.utils.interpolate_nearest(...
                hkEpochExtrapMargin, ...
                hkEpoch, ...
                InHk.Zv.HK_BIA_DIFF_GAIN, ...
                InSci.Zv.Epoch);



            %===============================
            % Derive HK_BIA_MODE_DIFF_PROBE
            %===============================
            % NOTE: FPA uint8 --> float-NaN --> FPA logical
            HkSciTime.Fpa_HK_BIA_MODE_DIFF_PROBE = bicas.utils.FillPositionsArray.floatNan2logical(...
                bicas.utils.interpolate_nearest(...
                    hkEpochExtrapMargin, ...
                    hkEpoch, ...
                    InHk.ZvFpa.HK_BIA_MODE_DIFF_PROBE.int2doubleNan(), ...
                    InSci.Zv.Epoch));



            % ASSERTIONS
            %irf.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN', 'HK_BIA_MODE_DIFF_PROBE'}, {})
            irf.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN', 'Fpa_HK_BIA_MODE_DIFF_PROBE'}, {})
        end



        % Utility function to shorten code.
        %
        % NOTE: Operates on entire ZvStruct since CALIBRATION_TABLE_INDEX exists
        % for L1R, but not L1, and the corresponding field may thus be or not be
        % present.
        function CALIBRATION_TABLE_INDEX = normalize_CALIBRATION_TABLE_INDEX(...
                ZvStruct, nRecords, inputDsi)

            C = bicas.classify_BICAS_L1_L1R_to_L2_DSI(inputDsi);

            if C.isL1r
                CALIBRATION_TABLE_INDEX = ZvStruct.CALIBRATION_TABLE_INDEX;
            elseif C.isL1
                CALIBRATION_TABLE_INDEX = nan(nRecords, 2);
            else
                error(...
                    ['Can not normalize CALIBRATION_TABLE_INDEX', ...
                    ' for this DSI classification.'])
            end

            irf.assert.sizes(CALIBRATION_TABLE_INDEX, [nRecords, 2])
        end



    end    % methods(Static, Access=public)



end
