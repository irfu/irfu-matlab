%
% Store SolO non-standard operations (NSO) table and supply basic functionality.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-09-22
%
classdef NSO_table   % < handle
    % PROPOSAL: Include read_ns_ops as static method.
    %   CON: Must include read_ns_ops's internal helper functions.
    %
    % PROPOSAL: Should somehow support logging. Which RCS NSO codes are used for
    % which time interval.
    %
    % PROPOSAL: Terminology similar to MMS, mission/mms/mms_ns_ops.xml, i.e.
    %           event = one time interval in NSO table.
    %   PRO: MMS table more recent design than Cluster table.
    % PROPOSAL: Terminology similar to Cluster, mission/cluster/caa/caa-control/ns_ops.xml, i.e.
    %           operation = one time interval in NSO table.
    %   PRO: Cluster table has been used more than MMS table.
    %
    % PROPOSAL: Some kind of feature for being able to trigger real NSO ID behaviour just for
    % testing, but only with proper setting.
    %   PROPOSAL: process_quality_filter_L2: Translate test NSO ID to proper NSO
    %   ID, but only with proper setting.
    
    properties(SetAccess=immutable, GetAccess=public)
        % NOTE: Same RCS NSO ID may occur multiple times. Not unique.
        startTt2000Array
        stopTt2000Array
        nsoIdCa
    end
    
    
    
    %#####################
    %#####################
    methods(Access=public)
    %#####################
    %#####################
    
        
        
        function obj = NSO_table(xmlFilePath)
            NsoTable = bicas.NSO_table.read_file(xmlFilePath);
            
            EJ_library.assert.sizes(...
                NsoTable.startTt2000Array, [-1], ...
                NsoTable.stopTt2000Array,  [-1], ...
                NsoTable.nsoIdCa,          [-1]);

            obj.startTt2000Array = NsoTable.startTt2000Array;
            obj.stopTt2000Array  = NsoTable.stopTt2000Array;
            obj.nsoIdCa          = NsoTable.nsoIdCa;
        end
        
        
        
        % Determine which RCS NSO IDs apply to which timestamps. Given e.g. a
        % zVar Epoch, obtain lists of indices to CDF records.
        %
        %
        % ARGUMENTS
        % =========
        % tt2000Array : Column array of TT2000 timestamps. Meant to be zVar
        %               Epoch.
        %
        %
        % RETURN VALUES
        % =============
        % bArraysCa   : Cell array of arrays of logical indices into tt2000Array.
        %               {iRcsNsoId}(iTimestamp)
        % nsoIdCa     : List of unique RCS NSO codes. {iRcsNsoId}
        % iEventArray : 1D array of indices into NSO event list. To be used for
        %               logging events that affect the submitted timestamps.
        %
        function [bArraysCa, nsoIdCa, iNsoArray] = get_NSO_timestamps(obj, tt2000Array)
            % PROPOSAL: Automatic tests.
            % PROPOSAL: Static method.
            %   PRO: More natural to have explicit table argument for automatic testing.
            % PROPOSAL: Return list to every event, not unique NSO IDs.
            %   PRO: Can return iStart, iStop instad of bArray.
            %       PRO: Uses less memory than in particular logical indexing.
            %       CON: Less practical. Want to use logical expressions separately record-by-record.
            %           CON-PROPOSAL: Should have function for turning i1:i2
            %                         into bArray.
            %   PRO: Can log every instance in NSO table.
            
            bEvents = EJ_library.utils.intervals_intersect(...
                obj.startTt2000Array, ...
                obj.stopTt2000Array, ...            
                min(tt2000Array), ...
                max(tt2000Array));
            
            % IMPLEMENTATION NOTE: Obtain subset of NSO table corresponding to
            % tt2000Array. Removes irrelevant RCS NSO events for the code after.
            startTt2000Array = obj.startTt2000Array(bEvents);
            stopTt2000Array  = obj.stopTt2000Array(bEvents);
            nsoIdCa          = obj.nsoIdCa(bEvents);
            
            iNsoArray        = find(bEvents);
            
            

            % IMPLEMENTATION NOTE: nsoIdCa is NOT a list of unique NSO IDs.
            % The return value must be a list of unique NSO IDs.
            % Must distinguish between these two.
            nNsoId    = numel(nsoIdCa);
            bArraysCa = cell(nNsoId, 1);
            for iNsoId = 1:nNsoId
                
                tt2000_1 = startTt2000Array(iNsoId);
                tt2000_2 = stopTt2000Array(iNsoId);
                
                b = false(size(tt2000Array));
                b((tt2000_1 <= tt2000Array) & (tt2000Array <= tt2000_2)) = true;
                
                bArraysCa{iNsoId} = b;
            end

        end    % get_NSO_timestamps
        
        
        
    end    % methods(Access = public)


    
    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################
    
    
    
        % Read SolO non-standard operations (NSO) XML file.
        %
        % IMPLEMENTATION NOTE: Separate static method only for testing purposes.
        %
        %
        % ARGUMENTS
        % =========
        % filePath
        %
        %
        % RETURN VALUES
        % =============
        % NsoTable : Struct of arrays representing file content. Not class.
        %
        %
        % Author: Erik P G Johansson, Uppsala, Sweden
        % First created 2020-09-21.
        %
        function NsoTable = read_file(filePath)
            % TODO-DEC: Time format? String? TT2000? Numeric?
            %   NOTE: Want to assert t1<t2.
            % PROPOSAL: Permit multiple forms of XML time input: (t1, t2), (t1,dt), (dt, t2)
            %   CON: Might not be able to convert XML to HTML using CSS.
            
            LEGAL_NSOID_CA = struct2cell(bicas.constants.NSOID);
            
            
            
            RootElem      = xmlread(filePath);
            TablesElem    = bicas.NSO_table.getXmlUniqChildElem(RootElem, 'table');
            EventElemList = TablesElem.getElementsByTagName('event');
            
            nEvents = EventElemList.getLength;
            
            startTt2000Array = int64(zeros(nEvents, 1));
            stopTt2000Array  = int64(zeros(nEvents, 1));
            nsoIdCa          = cell(nEvents, 1);
            
            for i = 1:nEvents
                EventElem = EventElemList.item(i-1);    % NOTE: Subtract by one.
                
                startUtc = bicas.NSO_table.getXmlChildElemStr(EventElem, 'startTimeUtc');
                stopUtc  = bicas.NSO_table.getXmlChildElemStr(EventElem, 'stopTimeUtc');
                nsoId    = bicas.NSO_table.getXmlChildElemStr(EventElem, 'rcsNsoId');
                
                startTt2000 = spdfparsett2000(startUtc);
                stopTt2000  = spdfparsett2000(stopUtc);
                
                % ASSERTION
                assert(ismember(nsoId, LEGAL_NSOID_CA), 'NSO table file contains illegal NSO ID="%s".', nsoId)
                % IMPLEMENTATION NOTE: This assertion requires converting the UTC
                % strings to a numerical format.
                assert(startTt2000 < stopTt2000, ...
                    'BICAS:read_ns_ops:FailedToReadInterpretNsOps', ...
                    'Start time does not preceed stop time for event stated to begin at UTC "%s".', ...
                    startUtc)
                
                startTt2000Array(i, 1) = startTt2000;
                stopTt2000Array(i, 1)  = stopTt2000;
                nsoIdCa{i, 1}          = nsoId;
            end
            
            NsoTable = struct(...
                'startTt2000Array', {startTt2000Array}, ...
                'stopTt2000Array',  {stopTt2000Array}, ...
                'nsoIdCa',          {nsoIdCa});
            
        end
        
    

    end    % methods(Static, Access=public)

    
    
    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################
        
        
        
        % Elem : Element that has exactly one child in the form of an element with
        %        specified tag name.
        function ChildElem = getXmlUniqChildElem(Elem, childTagName)
            ChildElemList = Elem.getElementsByTagName(childTagName);
            if ~(ChildElemList.getLength() == 1)
                error( ...
                    'BICAS:read_ns_ops:FailedToReadInterpretNsOps', ...
                    'XML element (tag name "%s") does not have exactly one child element with tag name "%s".', ...
                    Elem.getNodeName(), childTagName)
            end
            
            ChildElem = ChildElemList.item(0);
        end
        
        
        
        % Elem : Element that only has one child in the form of a text.
        %
        % NOTE: Probably does not really assert enough to ensure that the one element is
        % a text.
        function s = getXmlElemStr(Elem)
            ChildNodesList = Elem.getChildNodes();
            assert(ChildNodesList.getLength == 1, ...
                'BICAS:read_ns_ops:FailedToReadInterpretNsOps', ...
                'XML element does not have exactly one child node as expected.')
            
            s = char(ChildNodesList.item(0).getTextContent);
        end
        
        
        
        function s = getXmlChildElemStr(Elem, childTagName)
            ChildElem = bicas.NSO_table.getXmlUniqChildElem(Elem, childTagName);
            s         = bicas.NSO_table.getXmlElemStr(ChildElem);
        end
        
        
        
    end    % methods(Static, Access=private)



end
