%
% Collection of code relating to quality variables.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual   % < handle
    % PROPOSAL: Automatic test code.



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)
        % Overwrite selected data in selected CDF records with fill values/NaN.
        % Modify quality zVariables.
        %
        % NOTE: Almost does not modify PreDc.
        %   Exception: Modifies PreDc.Zv.QUALITY_FLAG
        %
        % Sets
        %   PreDc.Zv.QUALITY_FLAG (modifies)
        %   PostDc.Zv.L2_QUALITY_BITMASK
        %   PostDc.Zv.DemuxerOutput
        %   PostDc.Zv.currentAAmpere
        %
        %
        % RATIONALE
        % =========
        % Does NOT want to operate on structs that mimic the input or output
        % datasets, but on struct that are as similiar as possible for all forms
        % of L1R-->L2 processing.
        %
        function [PreDc, PostDc] = modify_quality_filter(...
                PreDc, PostDc, NsoTable, SETTINGS, L)

            % NOTE: Adds zVar L2_QUALITY_FLAG to PostDc, technically altering the PostDc format.
            %   NOTE: Also overwrites voltage with fill values.
            %   PROPOSAL: Treat output PostDc as another format?
            %   PROPOSAL: Initialize empty L2_QUALITY_FLAG when PostDc first created.
            %
            % PROPOSAL: Generalize function to be used in L3.
            %   CON: Can not be done since this function is meant to have access
            %        to arbitrary L1/L1R and L2 data to make decisions, although
            %        this is not much used yet.
            %
            % PROPOSAL: Automated test code.
            % PROPOSAL: Only return the modified zVariables, not PreDc & PostDc.
            %   Current (2023-08-11) effective input:
            %       PreDc.isLfr
            %       PreDc.Zv.Epoch
            %       PreDc.Zv.MUX_SET
            %       PreDc.Zv.ufv
            %       PreDc.Zv.QUALITY_FLAG
            %       PostDc.Zv.L2_QUALITY_BITMASK
            %   Current (2023-08-11) effective output:
            %       PreDc.Zv.QUALITY_FLAG
            %       PostDc.Zv.L2_QUALITY_BITMASK
            %       PostDc.Zv.DemuxerOutput
            %       PostDc.Zv.currentAAmpere
            %   --
            %   PRO: Makes modifications/output clearer.
            %   PRO: Simpler test code.
            %   PROPOSAL: Only have arguments for the required variables.
            %       PRO: Makes dependence clearer.
            %   PROPOSAL: Split up into two functions: (1) Set QUALITY_FLAG,
            %             L2_QUALITY_BITMASK, and (2) UFV.
            %
            % PROPOSAL: Call function from
            %       bicas.proc.L1L2.dc.process_calibrate_demux() (and redefine that
            %       function to include quality calculations.
            %   PRO: Can return modified values to caller which sets PostDc.
            %       PRO: Can avoid modifying PreDc.
            %       PRO: Can avoid modifying PostDc.
            %
            % PROPOSAL: Separate function for handling UFV.

            % ASSERTION
            assert(isa(PreDc,  'bicas.proc.L1L2.PreDc'))
            assert(isa(PostDc, 'bicas.proc.L1L2.PostDc'))
            nRecords = irf.assert.sizes(PreDc.Zv.Epoch, [-1]);



            % NOTE: Preallocates and ADDS zVar to PostDc.
            PostDc.Zv.L2_QUALITY_BITMASK = zeros(nRecords, 1, 'uint16');



            %============================================
            % Find CDF records to remove due to settings
            %============================================
            zvUfvSettings = bicas.proc.L1L2.qual.get_UFV_records_from_settings(...
                PreDc.Zv.Epoch, PreDc.Zv.MUX_SET, PreDc.isLfr, SETTINGS, L);

            zvUfv = PreDc.Zv.ufv | zvUfvSettings;



            %========================================
            % Take actions based on NSO events table
            %========================================
            % Variable naming convention:
            % CDF event    = NSO event that overlaps with CDF records.
            % Global event = NSO event in global NSO event table.

            % NOTE: iCdfEventNa = CDF events as indices to global events.
            [bCdfEventRecordsCa, cdfEventNsoIdCa, iCdfEventNa] = ...
                NsoTable.get_NSO_timestamps(PreDc.Zv.Epoch);
            nCdfEvents    = numel(cdfEventNsoIdCa);
            nGlobalEvents = numel(NsoTable.evtNsoIdCa);
            L.logf('info', ...
                ['Searched non-standard operations (NSO) table.', ...
                ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
                nCdfEvents, nGlobalEvents);

            % Index into LOCAL/CDF NSO events table.
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
                zv_QUALITY_FLAG       = PreDc.Zv.QUALITY_FLAG       (bCdfEventRecordsCa{kCdfEvent});
                zv_L2_QUALITY_BITMASK = PostDc.Zv.L2_QUALITY_BITMASK(bCdfEventRecordsCa{kCdfEvent});

                switch(eventNsoId)

                    case bicas.constants.NSOID.PARTIAL_SATURATION
                        zv_QUALITY_FLAG       = min(zv_QUALITY_FLAG, 1, 'includeNaN');
                        zv_L2_QUALITY_BITMASK = bitor(...
                            zv_L2_QUALITY_BITMASK, ...
                            bicas.constants.L2QBM_PARTIAL_SATURATION);

                    case bicas.constants.NSOID.FULL_SATURATION
                        zv_QUALITY_FLAG       = min(zv_QUALITY_FLAG, 0, 'includeNaN');
                        zv_L2_QUALITY_BITMASK = bitor(...
                            zv_L2_QUALITY_BITMASK, ...
                            bicas.constants.L2QBM_FULL_SATURATION);
                        zv_L2_QUALITY_BITMASK = bitor(...
                            zv_L2_QUALITY_BITMASK, ...
                            bicas.constants.L2QBM_PARTIAL_SATURATION);
                        % NOTE: Also set PARTIAL saturation bit when FULL
                        % saturation. /YK 2020-10-02.

                    case bicas.constants.NSOID.THRUSTER_FIRING
                        zv_QUALITY_FLAG = min(zv_QUALITY_FLAG, 1, 'includeNaN');
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
                PreDc.Zv.QUALITY_FLAG       (bCdfEventRecordsCa{kCdfEvent}) = zv_QUALITY_FLAG;
                PostDc.Zv.L2_QUALITY_BITMASK(bCdfEventRecordsCa{kCdfEvent}) = zv_L2_QUALITY_BITMASK;

            end    % for
            
            % IMPLEMENTATION NOTE: Reminder that bicas.proc.L1L2.PreDc can be
            % modified by code (depending on NSO table), despite that
            % bicas.proc.L1L2.PreDc should ideally be immutable but can
            % currently not be. Since modification only happens for NSO events,
            % this modification might not be run, depending on the time
            % interval. Therefore always running this "null modification" to
            % make sure that a mistakenly immutable bicas.proc.L1L2.PreDc always
            % triggers error.
            PreDc.Zv.QUALITY_FLAG = PreDc.Zv.QUALITY_FLAG;



            %=================================================================
            % Set zVariables for CURRENTS and VOLTAGES to NaN based on zvUfv.
            %=================================================================
            % Log
            logHeaderStr = sprintf(...
                ['All interval(s) of CDF records for which data should be set', ...
                ' to fill values (i.e. removed), regardless of reason.\n']);
            bicas.proc.L1L2.qual.log_UFV_records(PreDc.Zv.Epoch, zvUfv, logHeaderStr, L)
            %
            PostDc.Zv.currentAAmpere(zvUfv, :) = NaN;
            %
            fnCa = fieldnames(PostDc.Zv.DemuxerOutput);
            for iFn = 1:numel(fnCa)
                PostDc.Zv.DemuxerOutput.(fnCa{iFn})(zvUfv, :, :) = NaN;
            end



            % ASSERTION
            assert(isa(PreDc,  'bicas.proc.L1L2.PreDc'))
            assert(isa(PostDc, 'bicas.proc.L1L2.PostDc'))

        end    % modify_quality_filter



    end    % methods(Static)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



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
                zvEpoch, zv_MUX_SET, isLfr, SETTINGS, L)
            % PROPOSAL: Only derive UFV records based on settings. Not take
            %           previously found UFV records (BW) into account. Merging UFV
            %           records from settings and BW respectively can be done
            %           outside (trivial).
            % PROPOSAL: Separate function for logging which records that should be removed.

            bicas.utils.assert_ZV_Epoch(zvEpoch)
            assert(islogical(isLfr));

            %===============
            % Read settings
            %===============
            [muxModesRemove, settingMuxModesKey] = SETTINGS.get_fv(...
                'PROCESSING.L2.REMOVE_DATA.MUX_MODES');
            if     isLfr   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S';    % LFR
            else           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S';    % TDS
            end
            [removeMarginSec, settingMarginKey] = SETTINGS.get_fv(settingMarginKey);

            %==========================================
            % Find exact indices/CDF records to remove
            %==========================================
            zvUfv = irf.utils.true_with_margin(...
                zvEpoch, ...
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
            bicas.proc.L1L2.qual.log_UFV_records(zvEpoch, zvUfv, logHeaderStr, L)
        end



        % Log UFV records
        %
        % NOTE: Only logs (including header) if there are records to remove.
        function log_UFV_records(zvEpoch, zvUfv, logHeaderStr, L)
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
                    utc1  = irf.cdf.TT2000_to_UTC_str(zvEpoch(iCdfRecord1));
                    utc2  = irf.cdf.TT2000_to_UTC_str(zvEpoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s -- %s', ...
                        iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
            end

        end



    end    % methods(Static, Access=private)

end
