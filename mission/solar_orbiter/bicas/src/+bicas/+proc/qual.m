%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual
%
% Data structure for storing QRC-->quality ZV modification: QrdidSetting:
%   PROPOSAL: Replace QUALITY_FLAG and *_QUALITY_BITMASK with containers.Map:
%             key--> QUALITY_FLAG and *_QUALITY_BITMASK value.
%       PRO: Can simultaneously specify different behaviour different types of output datasets.
%       TODO-DEC: What should map key be?
%           PROPOSAL: DSI
%               CON: Too many identical entries.
%           PROPOSAL: Informal group of output datasets.
%               Ex: L2, L3_density
%       PROPOSAL: Class method with DSI as argument returns only relevant
%                 values.
%           CON: Can not handle (hypothetically) qualitatively different types
%                 of quality variables for different groups of output datasets,
%                 if method output should always be on the same format.
%               Ex: Density has L3_QUALITY_BITMASK but EFIELD and SCPOT do not
%                   (I think).
%               Ex: Hypothetical: Varying number of bits in *_QUALITY_BITMASK in
%                   different output datasets.
%               CON: There is no (real) such case.
%
% PROBLEM: How represent all translations from QRC to quality ZV modifications?
%   NOTE: QUALITY_FLAG has the same definition for all output datasets.
%   NOTE: L2_QUALITY_BITMASK and L3_QUALITY_BITMASK have different bit
%         definitions.
%   NOTE: L3_QUALITY_BITMASK *might* have different definitions for
%         different datasets in the future.
%   NOTE: Not all QRCs are defined for all datasets.
%       Ex: Density quality bit
%   NOTE: Not all datasets set quality bits at all.
%       Ex: EFIELD, SCPOT.
%       Ex: L2 CWF DSR
%   NOTE: Quality bits which are set for one dataset may effectively be
%         inherited from another output dataset in the processing.
%       Ex: DENSITY OSR --> DENSITY DSR
%       Ex: LFR CWF OSR --> LFR CWF DSR
%   NOTE: Quality bits are inherited from L2 to L3 (and should remain
%         unchanged).
%   PROPOSAL: In principle: One "quality function" per group of datasets
%             which sets quality variables.
%               [QUALITY_FLAG_array, L2/L3_QUALITY_BITMASK_array] = ...
%                   get_quality_ZVs(qrcArray1, ..., qrcArray2, [NsoTable])
%       PROBLEM: How reuse code that translates QRC arrays into different bitmask array?
%                (L2_QUALITY_BITMASK vs L3_QUALITY_BITMASK)
%           NOTE: Not a problem for QUALITY_FLAG



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % IMPLEMENTATION NOTE: Without allQrcidCa, the function can not create a
        % return value map that contains keys for all QRCIDs, in case the
        % NsoTable does not contain all QRCIDs.
        %
        % IMPLEMENTATION NOTE: allQrcidCa is an argument due to automated tests.
        % Could otherwise be derived from constants.
        %
        % ARGUMENTS
        % =========
        % allQrcidCa
        %       1D cell array of all QRCIDs.
        %
        % RETURN VALUE
        % ============
        % QrcFlagsMap
        %       containers.Map. QRCID->logical array
        %       Contains keys for all QRCIDs specified in allQrcidCa, not just
        %       those present in NsoTable.
        %
        function QrcFlagsMap = NSO_table_to_QRC_flag_arrays(...
                allQrcidCa, NsoTable, Epoch, L)

            % Local variable naming conventions:
            % ----------------------------------
            % GE = Global Event = NSO event in global NSO event table.
            % CE = CDF Event    = NSO event that overlaps with CDF records.
            % Ar                = (Non-cell) Array

            % NOTE: iCeAr = CDF events as indices to global events.
            [bCeRecordsCa, ceQrcidCa, iCeAr] = NsoTable.get_NSO_timestamps(Epoch);
            nCe = numel(ceQrcidCa);
            nGe = numel(NsoTable.evtQrcidCa);
            L.logf('info', ...
                ['Searched non-standard operations (NSO) table.', ...
                ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
                nCe, nGe);

            % Initialize "empty" QrcFlagsMap (all QRCs set to false)
            % ------------------------------------------------------
            % IMPLEMENTATION NOTE: valueType=logical implies scalar (sic!).
            QrcFlagsMap = containers.Map('keyType', 'char', 'valueType', 'any');
            for i = 1:numel(allQrcidCa)
                QrcFlagsMap(allQrcidCa{i}) = false(size(Epoch));
            end

            % Iterate over events in NSO events table which apply to the
            % specified timestamps (Epoch).
            for kCe = 1:nCe

                % Index into GLOBAL NSO events table.
                iGe        = iCeAr(kCe);
                eventQrcid = ceQrcidCa{kCe};
                % Indices into ZVs.
                bCeRecords = bCeRecordsCa{kCe};

                %===========================================================
                % Log the relevant NSO event in the GLOBAL NSO events table
                %===========================================================
                L.logf('info', '    %s -- %s %s', ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGe)), ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGe)), ...
                    eventQrcid);

                % ASSERTION
                % NOTE: Not perfect assertion on legal QRCIDs since code only
                % checks those relevant for the data (time interval) currently
                % processed. (Therefore also checks all QRCIDs when reads NSO
                % table.)
                assert(ismember(eventQrcid, allQrcidCa), 'Can not interpret QRCID "%s".', eventQrcid)

                %======================================
                % Set corresponding QRC array elements
                %======================================
                bQrc                    = QrcFlagsMap(eventQrcid);
                bQrc(bCeRecords)        = true;
                QrcFlagsMap(eventQrcid) = bQrc;
            end    % for
        end



        % NOTE: Does not return FPA, since internal algorithm can not produce
        % unknown values.
        %
        % ARGUMENTS
        % =========
        % nRec
        %       Number of CDF records (rows).
        %       IMPLEMENTATION NOTE: Needed for handling the case of zero
        %       QRCIDs.
        %
        function [QUALITY_FLAG, L2_QUALITY_BITMASK] = QRC_flag_arrays_to_quality_ZVs(...
                nRec, QrcFlagsMap, QrcidSettingsMap)

            % Create "empty" arrays
            QUALITY_FLAG               = ones( nRec, 1, 'uint8' ) * bicas.const.QUALITY_FLAG_MAX;
            L2_QUALITY_BITMASK         = zeros(nRec, 1, 'uint16');

            qrcidCa = QrcFlagsMap.keys();
            for i = 1:numel(qrcidCa)
                qrcid        = qrcidCa{i};
                QrcidSetting = QrcidSettingsMap(qrcid);
                bQrcid       = QrcFlagsMap(qrcid);

                assert(isequal( size(bQrcid), [nRec, 1] ))

                % Set QUALITY_FLAG
                % ----------------
                % IMPLEMENTATION NOTE: Only adjusts relevant indices since the
                % operation is more natural (simpler) that way.
                QUALITY_FLAG(bQrcid) = min(...
                    QUALITY_FLAG(bQrcid), ...
                    QrcidSetting.QUALITY_FLAG);

                % Set L2_QUALITY_BITMASK
                L2_QUALITY_BITMASK = bitor(...
                    L2_QUALITY_BITMASK, ...
                    QrcidSetting.L2_QUALITY_BITMASK * uint16(bQrcid));
            end
        end



    end    % methods(Static)



end
