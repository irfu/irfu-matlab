%
% Collection of code relating to quality variables.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual



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
        function ZvOut = modify_quality_filter(ZvIn, isLfr, NsoTable, SETTINGS, L)
            % PROPOSAL: Separate function for handling UFV.
            %   CON: Other quality variable processing might want to read or
            %        modify the UFV.
            % PROPOSAL: Separate function for modifying voltage & current using
            %           UFV, not for deriving it. Function called outside. This
            %           function (modify_quality_filter) only returns ZV UFV.
            %   PRO: Good for testing.
            %       PRO: Fewer variables in & out of function when testing.
            %
            % PROPOSAL: Structs for arguments & return values. -- IMPLEMENTED
            %   PRO: Safer w.r.t. confusing variables.
            %   CON: Can not as easily see in function which are the arguments & return values.
            %       CON: Easier to see arguments & return values when calling
            %            function, assuming that caller unpacks the return struct.
            %           CON: Caller may forget to unpack field in return value.
            %       CON-PROPOSAL: Explicitly convert between struct fields and
            %                     variables at beginning and end of function.
            %           CON: Longer code.
            %           PRO: Implementation is clear on what goes in and out of function.

            zv_Epoch         = ZvIn.Epoch;
            zv_MUX_SET       = ZvIn.MUX_SET;
            zv_QUALITY_FLAG  = ZvIn.QUALITY_FLAG;
            zvDemuxerOutput  = ZvIn.DemuxerOutput;
            zvCurrentAAmpere = ZvIn.currentAAmpere;
            zvUfv            = ZvIn.ufv;
            clear ZvIn

            % ASSERTIONS
            assert(isscalar(isLfr) && islogical(isLfr))
            nRecords = irf.assert.sizes( ...
                zv_Epoch,        [-1], ...
                zv_MUX_SET,      [-1], ...
                zv_QUALITY_FLAG, [-1], ...
                zvUfv,           [-1]);



            %============================================================
            % Find CDF records to remove due to settings and LFR ZV "BW"
            %============================================================
            zvUfvSettings = bicas.proc.L1L2.qual.get_UFV_records_from_settings(...
                zv_Epoch, zv_MUX_SET, isLfr, SETTINGS, L);

            zvUfv = zvUfv | zvUfvSettings;



            %========================================
            % Take actions based on NSO events table
            %========================================
            % Variable naming conventions:
            % CDF event    = NSO event that overlaps with CDF records.
            % Global event = NSO event in global NSO event table.
            % NA           = Numeric Array
            
            % NOTE: iCdfEventNa = CDF events as indices to global events.
            [bCdfEventRecordsCa, cdfEventNsoIdCa, iCdfEventNa] = ...
                NsoTable.get_NSO_timestamps(zv_Epoch);
            nCdfEvents    = numel(cdfEventNsoIdCa);
            nGlobalEvents = numel(NsoTable.evtNsoIdCa);
            L.logf('info', ...
                ['Searched non-standard operations (NSO) table.', ...
                ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
                nCdfEvents, nGlobalEvents);

            zv_L2_QUALITY_BITMASK = zeros(nRecords, 1, 'uint16');

            % Iterate over index into LOCAL/CDF NSO events table.
            for kCdfEvent = 1:nCdfEvents

                % Index into GLOBAL NSO events table.
                iGlobalEvent = iCdfEventNa(kCdfEvent);
                eventNsoId   = cdfEventNsoIdCa{kCdfEvent};

                %===========================================================
                % Log the relevant NSO event in the GLOBAL NSO events table
                %===========================================================
                L.logf('info', '    %s -- %s %s', ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGlobalEvent)), ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGlobalEvent)), ...
                    eventNsoId);



                %=================================
                % Take action depending on NSO ID
                %=================================
                % Temporary shorter variable name.
                zv_QUALITY_FLAG_cdfEvent       = zv_QUALITY_FLAG      (bCdfEventRecordsCa{kCdfEvent});
                zv_L2_QUALITY_BITMASK_cdfEvent = zv_L2_QUALITY_BITMASK(bCdfEventRecordsCa{kCdfEvent});

                switch(eventNsoId)

                    case bicas.constants.NSOID.PARTIAL_SATURATION
                        zv_QUALITY_FLAG_cdfEvent       = min(zv_QUALITY_FLAG_cdfEvent, 1, 'includeNaN');
                        zv_L2_QUALITY_BITMASK_cdfEvent = bitor(...
                            zv_L2_QUALITY_BITMASK_cdfEvent, ...
                            bicas.constants.L2QBM_PARTIAL_SATURATION);

                    case bicas.constants.NSOID.FULL_SATURATION
                        zv_QUALITY_FLAG_cdfEvent       = min(zv_QUALITY_FLAG_cdfEvent, 0, 'includeNaN');
                        zv_L2_QUALITY_BITMASK_cdfEvent = bitor(...
                            zv_L2_QUALITY_BITMASK_cdfEvent, ...
                            bicas.constants.L2QBM_FULL_SATURATION);
                        zv_L2_QUALITY_BITMASK_cdfEvent = bitor(...
                            zv_L2_QUALITY_BITMASK_cdfEvent, ...
                            bicas.constants.L2QBM_PARTIAL_SATURATION);
                        % NOTE: Also set PARTIAL saturation bit when FULL
                        % saturation. /YK 2020-10-02.

                    case bicas.constants.NSOID.THRUSTER_FIRING
                        zv_QUALITY_FLAG_cdfEvent = min(zv_QUALITY_FLAG_cdfEvent, 1, 'includeNaN');
                        % NOTE: There will be an L1 QUALITY_BITMASK bit for
                        % thruster firings eventually according to
                        % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
                        % Therefore(?) not setting any bit in
                        % L2_QUALITY_BITMASK. (YK 2020-11-03 did not ask for any
                        % to be set.)

                    otherwise
                        % ASSERTION
                        % NOTE: Not perfect assertion on legal NSO IDs since
                        % code only checks those relevant for the data (time
                        % interval) currently processed. (Therefore also checks
                        % all NSO IDs when reads NSO table.)
                        error('Can not interpret RCS NSO ID "%s".', ...
                            cdfEventNsoIdCa{kCdfEvent})

                end
                zv_QUALITY_FLAG      (bCdfEventRecordsCa{kCdfEvent}) = zv_QUALITY_FLAG_cdfEvent;
                zv_L2_QUALITY_BITMASK(bCdfEventRecordsCa{kCdfEvent}) = zv_L2_QUALITY_BITMASK_cdfEvent;

            end    % for



            %=================================================================
            % Set zVariables for CURRENTS and VOLTAGES to NaN based on zvUfv.
            %=================================================================
            [zvDemuxerOutput, zvCurrentAAmpere] = bicas.proc.L1L2.qual.set_voltage_current_fill_value(...
                zv_Epoch, zvDemuxerOutput, zvCurrentAAmpere, zvUfv, L);

            ZvOut = struct();
            ZvOut.DemuxerOutput      = zvDemuxerOutput;
            ZvOut.currentAAmpere     = zvCurrentAAmpere;
            ZvOut.QUALITY_FLAG       = zv_QUALITY_FLAG;
            ZvOut.L2_QUALITY_BITMASK = zv_L2_QUALITY_BITMASK;

        end    % modify_quality_filter



        % Overwrite selected records of voltage & current with fill values.
        %
        % ARGUMENTS
        % =========
        % zv_Epoch
        %       NOTE: Only needed for logging.
        function [zvDemuxerOutput, zvCurrentAAmpere] = set_voltage_current_fill_value(...
                zv_Epoch, zvDemuxerOutput, zvCurrentAAmpere, zvUfv, L)
            assert(islogical(zvUfv))

            % Log
            logHeaderStr = sprintf(...
                ['Interval(s) of CDF records for which data should be set', ...
                ' to fill values (i.e. removed), regardless of reason.\n']);
            bicas.proc.L1L2.qual.log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)

            % Set current values to fill value/NaN.
            zvCurrentAAmpere(zvUfv, :) = NaN;

            % Set voltage values to fill value/NaN.
            fnCa = fieldnames(zvDemuxerOutput);
            for iFn = 1:numel(fnCa)
                zvDemuxerOutput.(fnCa{iFn})(zvUfv, :, :) = NaN;
            end
        end



        % Find CDF records to remove based on settings (not data itself, almost,
        % since MUX mode is data).
        %
        % Ex: Sweeps
        %
        % ARGUMENTS
        % ---------
        % zv_MUX_SET
        %   Demultiplexer data, from BIAS HK or LFR.
        %
        function zvUfv = get_UFV_records_from_settings(...
                zv_Epoch, zv_MUX_SET, isLfr, SETTINGS, L)
            % PROPOSAL: Only derive UFV records based on settings. Not take
            %           previously found UFV records (BW) into account. Merging UFV
            %           records from settings and BW respectively can be done
            %           outside (trivial).
            % PROPOSAL: Separate function for logging which records that should be removed.

            bicas.utils.assert_ZV_Epoch(zv_Epoch)
            assert(islogical(isLfr));

            %===============
            % Read settings
            %===============
            [muxModesRemove, settingMuxModesKey] = SETTINGS.get_fv(...
                'PROCESSING.L2.REMOVE_DATA.MUX_MODES');
            muxModesRemove = muxModesRemove(:);
            if     isLfr   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S';    % LFR
            else           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S';    % TDS
            end
            [removeMarginSec, settingMarginKey] = SETTINGS.get_fv(settingMarginKey);

            %==========================================
            % Find exact indices/CDF records to remove
            %==========================================
            zvUfv = irf.utils.true_with_margin(...
                zv_Epoch, ...
                ismember(zv_MUX_SET, muxModesRemove), ...
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
                settingMuxModesKey, ...
                strjoin(irf.str.sprintf_many('%g', muxModesRemove), ', '), ...
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
