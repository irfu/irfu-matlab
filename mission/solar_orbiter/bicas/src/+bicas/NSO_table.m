%
% Class for storing one SolO NSO table and supply basic functionality.
%
%
% TERMINOLOGY
% ===========
% NSO       : Non-Standard Operations
% NSO event : One continuous time interval with exactly one NSO ID attached to
%             it.
% NSO ID    : Unique string that identifies the actions BICAS should take.
%             Multiple NSO events may have the same NSO ID.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-22
%
classdef NSO_table    
    % PROPOSAL: New name
    %   PROPOSAL: ~events
    %   PROPOSAL: ~list
    %       PRO: Table implies multiple dimensions more than list.
    %            List implies ~1D (or even ~set).
    %       NOTE: solo_ns_ops.xml.html:
    %             Title "Non-Standard Operations (NSO) List"
    %
    % PROPOSAL: Terminology similar to MMS, mission/mms/mms_ns_ops.xml, i.e.
    %           event = one time interval in NSO table.
    %   PRO: MMS table more recent design than Cluster table.
    % PROPOSAL: Terminology similar to Cluster, mission/cluster/caa/caa-control/ns_ops.xml, i.e.
    %           operation = one time interval in NSO table.
    %   PRO: Cluster table has been used more than MMS table.
    %
    % PROPOSAL: Some kind of feature for being able to trigger real NSO ID
    %           behaviour just for testing, but only with proper setting.
    %   PROPOSAL: process_quality_filter_L2: Translate test NSO ID to proper NSO
    %   ID, but only with proper setting.
    %
    % PROPOSAL: Automatic test code for get_NSO_timestamps().
    %   NOTE: Implies separating file reading from initializing object.
    %
    % PROPOSAL: Ability to add time margins for selected NSOIDs.
    %   PRO: Easier to change pre-existing margins than to modify the XML file.
    %       CON: Harder for outsiders to interpret (& edit) the XML file.
    %   Ex: thruster_firings.
    
    
    
    properties(SetAccess=immutable, GetAccess=public)
        % NOTE: Same RCS NSO ID may occur multiple times. Not unique.
        evtStartTt2000Array
        evtStopTt2000Array
        evtNsoIdCa
    end
    
    
    
    %#####################
    %#####################
    methods(Access=public)
    %#####################
    %#####################



        function obj = NSO_table(xmlFilePath)
            NsoTable = bicas.NSO_table.read_file(xmlFilePath);



            %============
            % ASSERTIONS
            %============
            % PROPOSAL: Move ~all assertions to bicas.NSO_table.read_file ?
            % PROPOSAL: Collect ~all assertions, in constructor (here) and in
            %           bicas.NSO_table.read_file ?
            % PROPOSAL: Check that FULL_SATURATION and PARTIAL SATURATION do not overlap.
            
            irf.assert.sizes(...
                NsoTable.evtStartTt2000Array, [-1], ...
                NsoTable.evtStopTt2000Array,  [-1], ...
                NsoTable.evtNsoIdCa,          [-1]);
            
            % IMPLEMENTATION NOTE: Can not assume that both start & stop
            % timestaps are sorted. One event may entirely contain another
            % event (with different NSO ID).
            % IMPLEMENTATION NOTE: Can not assume "strictly ascending" values,
            % since events with separate NSO IDs may begin at the exact same
            % instant.
            if ~issorted(NsoTable.evtStartTt2000Array)
                iEvt = find(diff(NsoTable.evtStartTt2000Array) < 0) + 1;
                assert(~isempty(iEvt));
                
                utcCa = irf.cdf.TT2000_to_UTC_str_many(...
                    NsoTable.evtStartTt2000Array(iEvt));
                
                sCa = irf.str.sprintf_many('    %s\n', utcCa);
                timestampsListStr = strjoin(sCa);
                
                error('BICAS:FailedToReadInterpretNsOps', ...
                    ['NsoTable.evtStartTt2000Array is not sorted. Events', ...
                    ' beginning at the following timestamps begin earlier', ...
                    ' than the precedeing events:\n%s'], ...
                    timestampsListStr);
            end
            
            %----------------------------------------------------------------
            % ASSERTION: Events with the same NSO ID do not overlap (and are
            % time sorted)
            %----------------------------------------------------------------
            uniqueEvtNsoIdCa = unique(NsoTable.evtNsoIdCa);
            for i = 1:numel(uniqueEvtNsoIdCa)
                nsoId = uniqueEvtNsoIdCa{i};
                b = strcmp(nsoId, NsoTable.evtNsoIdCa);
                
                % Sorted (earlier assertion), but not monotonically.
                evtStartTt2000Array = NsoTable.evtStartTt2000Array(b);
                % Can not be assumed to be sorted (no earlier assertion).
                evtStopTt2000Array  = NsoTable.evtStopTt2000Array(b);
                
                % NOTE: This will catch some overlaps, but not all.
                assert(issorted(evtStartTt2000Array, 'strictascend'), ...
                    ['evtStartTt2000Array for nsoId="%s" is not time-sorted. ', ...
                    'At least two events with that NSO ID overlap.'], nsoId)
                assert(issorted(evtStopTt2000Array, 'strictascend'), ...
                    ['evtStopTt2000Array for nsoId="%s" is not time-sorted. ', ...
                    'At least two events with that NSO ID overlap.'], nsoId)
            
                % ASSERTION: Events do not overlap
                % --------------------------------
                % NOTE: ASSUMPTION: Start timestamps are already time-sorted.
                % NOTE: Transposing before 2D-->1D vector.
                % NOTE: 'strictascend' excludes ~adjacent events.
                temp = [...
                    NsoTable.evtStartTt2000Array(b), ...
                    NsoTable.evtStopTt2000Array(b)]';
                tt2000Array = temp(:);
                assert(issorted(tt2000Array, 'strictascend'), ...
                    ['At least two events for nsoId="%s"', ...
                    ' seem to overlap with each other.'], nsoId)
            end
            
            

            obj.evtStartTt2000Array = NsoTable.evtStartTt2000Array;
            obj.evtStopTt2000Array  = NsoTable.evtStopTt2000Array;
            obj.evtNsoIdCa          = NsoTable.evtNsoIdCa;
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
        % bEvtArraysCa
        %       Cell array of arrays of logical indices into tt2000Array.
        %       {iEvent}(iTimestamp)
        % evtNsoIdCa
        %       List of NSO IDs. {iEvent}
        % iGlobalEventsArray
        %       1D array of indices into NSO event list. Can be used for logging
        %       the tabulated (global) NSO events that affect e.g. a particular
        %       CDF.
        %
        function [bEvtArraysCa, evtNsoIdCa, iGlobalEventsArray] = get_NSO_timestamps(obj, tt2000Array)
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
            
            bEvents = irf.utils.intervals_intersect(...
                obj.evtStartTt2000Array, ...
                obj.evtStopTt2000Array, ...            
                min(tt2000Array), ...
                max(tt2000Array), ...
                'closed intervals');

            % IMPLEMENTATION NOTE: Obtain SUBSET of NSO table EVENTS which
            % overlap with timestamps in tt2000Array. Indirectly also removes
            % irrelevant RCS NSO IDs (not just events) for the code after.
            evtStartTt2000Array = obj.evtStartTt2000Array(bEvents);
            evtStopTt2000Array  = obj.evtStopTt2000Array(bEvents);
            evtNsoIdCa          = obj.evtNsoIdCa(bEvents);
            
            iGlobalEventsArray  = find(bEvents);



            % IMPLEMENTATION NOTE: evtNsoIdCa is NOT a list of unique NSO IDs.
            % The return value must be a list of unique NSO IDs.
            % Must distinguish between these two.
            nEvents      = numel(evtNsoIdCa);
            bEvtArraysCa = cell(nEvents, 1);
            for iEvent = 1:nEvents
                
                tt2000_1 = evtStartTt2000Array(iEvent);
                tt2000_2 = evtStopTt2000Array(iEvent);
                
                b = false(size(tt2000Array));
                b((tt2000_1 <= tt2000Array) & (tt2000Array <= tt2000_2)) = true;
                
                bEvtArraysCa{iEvent, 1} = b;
            end

            
            
            % ASSERTIONS
            irf.assert.sizes(...
                bEvtArraysCa,       [-1], ...
                evtNsoIdCa,         [-1], ...
                iGlobalEventsArray, [-1])
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
        % Author: Erik P G Johansson, IRF, Uppsala, Sweden
        % First created 2020-09-21.
        %
        function NsoTable = read_file(filePath)
            % TODO-DEC: Time format? String? TT2000? Numeric?
            %   NOTE: Want to assert t1<t2.
            % PROPOSAL: Permit multiple forms of XML time input: (t1, t2), (t1,dt), (dt, t2)
            %   CON: Might not be able to convert XML to HTML using CSS.
            
            % List of all legal NSO IDs.
            LEGAL_NSOID_CA = struct2cell(bicas.constants.NSOID);
            irf.assert.castring_set(LEGAL_NSOID_CA)
            
            
            
            RootXmlElem      = xmlread(filePath);
            MainXmlElem      = bicas.NSO_table.getXmlUniqChildElem(RootXmlElem, 'main');
            TablesXmlElem    = bicas.NSO_table.getXmlUniqChildElem(MainXmlElem, 'eventsTable');
            EventXmlElemList = TablesXmlElem.getElementsByTagName('event');
            
            nEvents = EventXmlElemList.getLength;
            
            evtStartTt2000Array = int64(zeros(nEvents, 1));
            evtStopTt2000Array  = int64(zeros(nEvents, 1));
            evtNsoIdCa          = cell(nEvents, 1);
            
            for i = 1:nEvents
                % NOTE: Subtract by one.
                EventXmlElem = EventXmlElemList.item(i-1);
                
                startUtc = bicas.NSO_table.getXmlChildElemStr(EventXmlElem, 'startTimeUtc');
                stopUtc  = bicas.NSO_table.getXmlChildElemStr(EventXmlElem, 'stopTimeUtc');
                nsoId    = bicas.NSO_table.getXmlChildElemStr(EventXmlElem, 'rcsNsoId');
                
                
                
                startTt2000 = spdfparsett2000(startUtc);
                stopTt2000  = spdfparsett2000(stopUtc);
                
                % ASSERTIONS
                assert(ismember(nsoId, LEGAL_NSOID_CA), ...
                    'NSO table file contains illegal NSO ID="%s".', nsoId)
                % IMPLEMENTATION NOTE: This assertion requires converting the
                % UTC strings to a numerical format.
                assert(startTt2000 < stopTt2000, ...
                    'BICAS:FailedToReadInterpretNsOps', ...
                    ['Start time does not precede stop time for', ...
                    ' NSO table event stated to begin at UTC "%s".'], ...
                    startUtc)
                
                evtStartTt2000Array(i, 1) = startTt2000;
                evtStopTt2000Array(i, 1)  = stopTt2000;
                evtNsoIdCa{i, 1}          = nsoId;
            end
            
            NsoTable = struct(...
                'evtStartTt2000Array', {evtStartTt2000Array}, ...
                'evtStopTt2000Array',  {evtStopTt2000Array}, ...
                'evtNsoIdCa',          {evtNsoIdCa});
            
        end
        
    

    end    % methods(Static, Access=public)

    
    
    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################
        
        
        
        % Elem
        %   Element that has exactly one child in the form of an element
        %   with the specified tag name.
        function ChildXmlElem = getXmlUniqChildElem(XmlElem, childTagName)
            
            ChildXmlElemList = XmlElem.getElementsByTagName(childTagName);
            if ~(ChildXmlElemList.getLength() == 1)
                error( ...
                    'BICAS:FailedToReadInterpretNsOps', ...
                    ['XML element (tag name "%s") does not have exactly, ', ...
                    ' one child element with tag name "%s" as expected.'], ...
                    XmlElem.getNodeName(), childTagName)
            end
            
            ChildXmlElem = ChildXmlElemList.item(0);
        end
        
        
        
        % XmlElem : Element that only has one child in the form of a text.
        %
        % NOTE: Probably does not really assert enough to ensure that the one
        % element is a text.
        function s = getXmlElemStr(XmlElem)
            
            ChildXmlNodesList = XmlElem.getChildNodes();
            assert(ChildXmlNodesList.getLength == 1, ...
                'BICAS:FailedToReadInterpretNsOps', ...
                'XML element does not have exactly one child node as expected.')
            
            s = char(ChildXmlNodesList.item(0).getTextContent);
        end
        
        
        
        function s = getXmlChildElemStr(XmlElem, childTagName)
            ChildXmlElem = bicas.NSO_table.getXmlUniqChildElem(XmlElem, childTagName);
            s            = bicas.NSO_table.getXmlElemStr(ChildXmlElem);
        end
        
        
        
    end    % methods(Static, Access=private)



end
