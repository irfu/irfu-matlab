%
% UNFINISHED / EXPERIMENTAL
%
% Store SolO non-standard operations (NSO) table and supply basic functionality.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-09-22
%
classdef nso_table   % < handle
    % PROPOSAL: Include file reading code.
    %
    % PROPOSAL: Should somehow support logging. Which RCS NSO codes are used for
    % which time interval.
    
    properties(SetAccess = immutable)
        % NOTE: Same RCS NSO code may occur multiple times. Not unique.
        startTt2000Array
        stopTt2000Array
        rcsNsoCodeCa         % NOTE: Not unique strings.
    end
    
    
    
    methods(Access = public)

        
        
        function obj = nso_table(xmlFilePath)
            NsoTable = bicas.read_ns_ops(xmlFilePath);
            
            EJ_library.assert.sizes(...
                NsoTable.startTt2000Array, [-1], ...
                NsoTable.stopTt2000Array,  [-1], ...
                NsoTable.rcsNsoCodeCa,     [-1]);

            obj.startTt2000Array = NsoTable.startTt2000Array;
            obj.stopTt2000Array  = NsoTable.stopTt2000Array;
            obj.rcsNsoCodeCa     = NsoTable.rcsNsoCodeCa;
        end
        
        
        
        % Determine which RCS NSO codes apply to which timestamps.
        %
        % RETURN VALUES
        % =============
        % bArraysCa    : Cell array of arrays of logical indices into tt2000Array.
        %                {iRcsNsoCode}(iTimestamp)
        % rcsNsoCodeCa : List of unique RCS NSO codes. {iRcsNsoCode}
        % iEventArray  : 1D array of indices into NSO event list. To be used for
        %                logging events that affect the submitted timestamps.
        function [bArraysCa, rcsNsoCodeCa, iEventArray] = get_NSO_timestamps(obj, tt2000Array)
            % PROPOSAL: Return b instead i.
            % PROPOSAL: Automatic tests.
            
            bEvents = EJ_library.utils.intervals_intersect(...         
                obj.startTt2000Array, ...
                obj.stopTt2000Array, ...            
                min(tt2000Array), ...
                max(tt2000Array));
            
            % IMPLEMENTATION NOTE: Reduce size of arrays to possibly improve
            % speed if there are many NSO events (cf. Cluster). Maybe pointless.
            % Also removes irrelevant RCS NSO events for the code after.
            startTt2000Array = obj.startTt2000Array(bEvents);
            stopTt2000Array  = obj.stopTt2000Array(bEvents);
            rcsNsoCodeCa     = obj.rcsNsoCodeCa(bEvents);
            
            iEventArray      = find(bEvents);
            
            % IMPLEMENTATION NOTE: rcsNsoCodeCa is NOT a list of unique NSO RCS
            % codes. The return value must be a list of unique NSO RCS codes.
            % Must distinguish between these two.
            uniqCodesCa = unique(rcsNsoCodeCa);
            nUniqCodes  = numel(uniqCodesCa);
            bArraysCa   = cell(nUniqCodes, 1);
            for iUniqCode = 1:nUniqCodes
                iCodeEventArray = find(strcmp(uniqCodesCa{iUniqCode}, rcsNsoCodeCa));
                
                b = false(size(tt2000Array));
                for iEvent = iCodeEventArray(:)'
                    b = b || EJ_library.utils.intervals_intersect(...
                        startTt2000Array(iEvent), ...
                        stopTt2000Array(iEvent), ...
                        tt2000Array, ...
                        tt2000Array);
                end
                bArraysCa{iUniqCode} = b;
            end
            
            rcsNsoCodeCa = uniqCodesCa;   % NOTE: Overwrites old variable.

        end
        
        
        
    end
    
    
        
%     methods(Access = private)
%     end
    
    
        
%     methods(Static)
%     end
    
end