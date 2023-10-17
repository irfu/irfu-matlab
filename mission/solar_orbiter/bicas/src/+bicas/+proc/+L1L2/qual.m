%
% Collection of code relating to quality variables for L1/L1R to L2 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual
    % PROPOSAL: Redefine as applying to quality variables for all processing
    %           (between any archiving levels).



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Set quality zVariables.
        % Overwrite selected data in selected CDF records with fill values/NaN.
        %
        % NOTE: Does not handle PROCESSING.ZV_QUALITY_FLAG_MAX. That is handled
        %       by bicas.write_dataset_CDF().
        % NOTE: Does not seem able to ever set zv_L2_QUALITY_BITMASK to fill
        %       value.
        function [zvUfv, zv_QUALITY_FLAG_Fpa, zv_L2_QUALITY_BITMASK] = ...
                modify_quality_filter(ZvIn, isLfr, NsoTable, SETTINGS, L)

            irf.assert.struct(ZvIn, {'Epoch', 'ufv', 'bdmFpa', 'QUALITY_FLAG_Fpa'}, {})
            zv_Epoch            = ZvIn.Epoch;
            zvUfv               = ZvIn.ufv;
            zvBdmFpa            = ZvIn.bdmFpa;
            zv_QUALITY_FLAG_Fpa = ZvIn.QUALITY_FLAG_Fpa;
            clear ZvIn

            % ASSERTIONS
            assert(isscalar(isLfr) && islogical(isLfr))
            assert(isa(zv_QUALITY_FLAG_Fpa, 'bicas.utils.FPArray'))
            nRecords = irf.assert.sizes( ...
                zv_Epoch,            [-1], ...
                zvBdmFpa,            [-1], ...
                zv_QUALITY_FLAG_Fpa, [-1], ...
                zvUfv,               [-1]);

            % Pre-allocate
            zv_L2_QUALITY_BITMASK = zeros(nRecords, 1, 'uint16');


            %============================================================
            % Find CDF records to remove due to settings and LFR ZV "BW"
            %============================================================
            zvUfvSettings = bicas.proc.L1L2.qual.get_UFV_records_from_settings(...
                zv_Epoch, zvBdmFpa, isLfr, SETTINGS, L);

            zvUfv = zvUfv | zvUfvSettings;



            %========================================
            % Take actions based on NSO events table
            %========================================
            [zv_QUALITY_FLAG_CapFpa, zv_L2_QUALITY_BITMASK_new] = bicas.proc.L1L2.qual.get_quality_by_NSOs(...
                zv_Epoch, NsoTable, L);
            zv_QUALITY_FLAG_Fpa   = zv_QUALITY_FLAG_Fpa.min(zv_QUALITY_FLAG_CapFpa);
            zv_L2_QUALITY_BITMASK = bitor(zv_L2_QUALITY_BITMASK, zv_L2_QUALITY_BITMASK_new);



            assert(isa(zv_L2_QUALITY_BITMASK, 'uint16'))
        end    % modify_quality_filter
        
        
        
        % Derive QUALITY_FLAG *cap* and L2_QUALITY_FLAG *bits* set due to NSO
        % table by itself. Return values are then supposed to be used for
        % creating global versions of the actual ZVs.
        %
        %
        % RETURN VALUES
        % =============
        % QUALITY_FLAG_CapFpa
        %       FPA. Cap (highest allowed value) for ZV QUALITY_FLAG.
        % L2_QUALITY_BITMASK
        %       Array. L2_QUALITY_BITMASK bits set based on NSOs only. Should be
        %       merged (OR:ed) with global L2_QUALITY_BITMASK.
        %
        function [QUALITY_FLAG_CapFpa, L2_QUALITY_BITMASK] = get_quality_by_NSOs(...
                zv_Epoch, NsoTable, L)

            % PROPOSAL: Arguments for the QUALITY_FLAG cap and
            %           L2_QUALITY_BITMASK values to set for the respective NSOIDs.
            %   PRO: Better for testing.
            %   PRO: Can better document the consequences of different NSOIDs.
            %       PROPOSAL: Can document as ~constants.
            %   PROPOSAL: Map argument.

            % Variable naming conventions:
            % CE = CDF Event    = NSO event that overlaps with CDF records.
            % GE = Global Event = NSO event in global NSO event table.
            % NA                = Numeric Array

            % NOTE: iCdfEventNa = CDF events as indices to global events.
            [bCeRecordsCa, ceNsoidCa, iCeNa] = ...
                NsoTable.get_NSO_timestamps(zv_Epoch);
            nCe = numel(ceNsoidCa);
            nGe = numel(NsoTable.evtNsoidCa);
            L.logf('info', ...
                ['Searched non-standard operations (NSO) table.', ...
                ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
                nCe, nGe);



            % Pre-allocate
            QUALITY_FLAG_CapFpa = bicas.utils.FPArray(...
                bicas.const.QUALITY_FLAG_MAX * ones(size(zv_Epoch), 'uint8'), ...
                'NO_FILL_POSITIONS');
            L2_QUALITY_BITMASK = zeros(size(zv_Epoch), 'uint16');

            % Iterate over index into LOCAL/CDF NSO events table.
            for kCe = 1:nCe

                % Index into GLOBAL NSO events table.
                iGe = iCeNa(kCe);
                eventNsoid = ceNsoidCa{kCe};
                % Indices into ZVs.
                bCeRecords = bCeRecordsCa{kCe};

                %===========================================================
                % Log the relevant NSO event in the GLOBAL NSO events table
                %===========================================================
                L.logf('info', '    %s -- %s %s', ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGe)), ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGe)), ...
                    eventNsoid);

                %================================
                % Take action depending on NSOID
                %================================
                % Temporary shorter variable name.
                L2_QUALITY_BITMASK_ce = L2_QUALITY_BITMASK(bCeRecords);

                switch(eventNsoid)

                    case bicas.const.NSOID.PARTIAL_SATURATION

                        QUALITY_FLAG_capCe = bicas.const.QUALITY_FLAG_CAP_PARTIAL_SATURATION;

                        L2_QUALITY_BITMASK_ce = bitor(...
                            L2_QUALITY_BITMASK_ce, ...
                            bicas.const.L2QBM_PARTIAL_SATURATION);

                    case bicas.const.NSOID.FULL_SATURATION

                        QUALITY_FLAG_capCe = bicas.const.QUALITY_FLAG_CAP_FULL_SATURATION;

                        % NOTE: Also set PARTIAL saturation bit when FULL
                        % saturation. /YK 2020-10-02.
                        L2_QUALITY_BITMASK_ce = bitor(...
                            L2_QUALITY_BITMASK_ce, ...
                            bicas.const.L2QBM_FULL_SATURATION + ...
                            bicas.const.L2QBM_PARTIAL_SATURATION);

                    case bicas.const.NSOID.THRUSTER_FIRING

                        QUALITY_FLAG_capCe = bicas.const.QUALITY_FLAG_CAP_THRUSTER_FIRING;

                        % NOTE: There will be an L1 QUALITY_BITMASK bit for
                        % thruster firings eventually according to
                        % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
                        % Therefore(?) not setting any bit in
                        % L2_QUALITY_BITMASK. (YK 2020-11-03 did not ask for any
                        % to be set.)

                    otherwise
                        % ASSERTION
                        % NOTE: Not perfect assertion on legal NSOIDs since
                        % code only checks those relevant for the data (time
                        % interval) currently processed. (Therefore also checks
                        % all NSOIDs when reads NSO table.)
                        error('Can not interpret RCS NSOID "%s".', ...
                            ceNsoidCa{kCe})

                end
                
                % Old implementation which should have bug which tests do not
                % test: Can not handle time-overlapping NSO events (for
                % different NSOIDs).
%                 QUALITY_FLAG_CapFpa(bCdfEvent, 1) = bicas.utils.FPArray(QUALITY_FLAG_capCdfEvent, 'NO_FILL_POSITIONS');
%                 L2_QUALITY_BITMASK( bCdfEvent, 1) = L2_QUALITY_BITMASK_cdfEvent;
                % New implementation which should not have abovementioned bug.
                QUALITY_FLAG_CapFpaTemp            = QUALITY_FLAG_CapFpa(bCeRecords, 1);
                QUALITY_FLAG_CapFpa(bCeRecords, 1) = QUALITY_FLAG_CapFpaTemp.min(bicas.utils.FPArray(QUALITY_FLAG_capCe, 'NO_FILL_POSITIONS'));
                L2_QUALITY_BITMASK( bCeRecords, 1) = bitor(L2_QUALITY_BITMASK(bCeRecords, 1), L2_QUALITY_BITMASK_ce);
                
            end    % for

        end



        % Overwrite selected records of voltage & current with FVs.
        %
        % ARGUMENTS
        % =========
        % zv_Epoch
        %       NOTE: Only needed for logging.
        % zvAsrSamplesAVoltSrm
        %       ASR samples.
        %       NOTE: Handle object which is MODIFIED.
        function zvCurrentAAmpere = set_voltage_current_FV(...
                zv_Epoch, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, zvUfv, L)
            % PROPOSAL: Separate functions for ASR samples and bias currents.
            
            assert(islogical(zvUfv))
            assert(isa(zvAsrSamplesAVoltSrm, 'bicas.utils.SameRowsMap'))

            % Log
            logHeaderStr = sprintf(...
                ['Interval(s) of CDF records for which data should be set', ...
                ' to fill values (i.e. removed), regardless of reason.\n']);
            bicas.proc.L1L2.qual.log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)

            % Set current values to fill value/NaN.
            zvCurrentAAmpere(zvUfv, :) = NaN;

            % ====================================
            % Set voltage values to fill value/NaN
            % ====================================
            % NOTE: Should really use future bicas.utils.SameSizeTypeMap here
            %       which contains size on other dimensions.
            keysCa = zvAsrSamplesAVoltSrm.keys;
            nSpr   = size(zvAsrSamplesAVoltSrm(keysCa{1}), 2);

            % IMPLEMENTATION NOTE: bicas.utils.SameRowsMap.setRows() can not
            % handle logical indexing.
            iUfv = find(zvUfv);
            nanArray = NaN(size(iUfv, 1), nSpr);
            tempSrm = bicas.utils.SameRowsMap(...
                'char', size(nanArray, 1), ...
                'CONSTANT', nanArray, zvAsrSamplesAVoltSrm.keys);
            zvAsrSamplesAVoltSrm.setRows(tempSrm, iUfv);
        end



        % Find CDF records to remove based on settings (not data itself, almost,
        % since MUX mode is data).
        %
        % Ex: Sweeps
        %
        % ARGUMENTS
        % ---------
        % zvBdm
        %       Demultiplexer data, from BIAS HK or LFR.
        %       Fill positions are not recognized as BDMs for setting UFV.
        %
        function zvUfv = get_UFV_records_from_settings(...
                zv_Epoch, zvBdmFpa, isLfr, SETTINGS, L)
            % PROPOSAL: Only derive UFV records based on settings. Not take
            %           previously found UFV records (BW) into account. Merging UFV
            %           records from settings and BW respectively can be done
            %           outside (trivial).
            % PROPOSAL: Separate function for logging which records that should be removed.

            bicas.utils.assert_ZV_Epoch(zv_Epoch)
            assert(islogical(isLfr));
            assert(isa(zvBdmFpa, 'bicas.utils.FPArray'))

            %===============
            % Read settings
            %===============
            [bdmRemoveArray, settingBdmRemoveKey] = SETTINGS.get_fv(...
                'PROCESSING.L2.REMOVE_DATA.MUX_MODES');
            bdmRemoveArray = bdmRemoveArray(:);
            if     isLfr   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S';    % LFR
            else           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S';    % TDS
            end
            [removeMarginSec, settingMarginKey] = SETTINGS.get_fv(settingMarginKey);

            %==========================================
            % Find exact indices/CDF records to remove
            %==========================================
            zvUfv = irf.utils.true_with_margin(...
                zv_Epoch, ...
                ismember(zvBdmFpa.int2doubleNan(), bdmRemoveArray), ...
                removeMarginSec * 1e9);

            %=====
            % Log
            %=====
            logHeaderStr = sprintf(...
                ['Found interval(s) of CDF records for which data should be set to', ...
                ' fill values (i.e. removed) based on settings.\n', ...
                '    NOTE: This may not be all CDF records which will be removed.\n', ...
                '    Setting %s = [%s]\n', ...
                '    Setting %s = %f\n'], ...
                settingBdmRemoveKey, ...
                strjoin(irf.str.sprintf_many('%g', bdmRemoveArray), ', '), ...
                settingMarginKey, ...
                removeMarginSec);
            bicas.proc.L1L2.qual.log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)
        end



    end    % methods(Static)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Log UFV records
        %
        % NOTE: Only logs (including header) if there are records to remove.
        function log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)
            LL = 'info';    % LL = Log Level

            [i1Array, i2Array] = irf.utils.split_by_false(zvUfv);
            nUfvIntervals = numel(i1Array);
            if nUfvIntervals > 0

                %==============
                % Log settings
                %==============
                L.logf(LL, logHeaderStr)

                %===============
                % Log intervals
                %===============
                for iRi = 1:nUfvIntervals
                    iCdfRecord1 = i1Array(iRi);
                    iCdfRecord2 = i2Array(iRi);
                    utc1  = irf.cdf.TT2000_to_UTC_str(zv_Epoch(iCdfRecord1));
                    utc2  = irf.cdf.TT2000_to_UTC_str(zv_Epoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s -- %s', ...
                        iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
            end

        end



    end    % methods(Static, Access=private)

end
