%
% Class for storing one SolO NSO table and supply basic functionality.
%
%
% TERMINOLOGY
% ===========
% NSO event : One continuous time interval with exactly one QRCID attached to
%             it.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-22
%
classdef NsoTable
  % PROPOSAL: New name
  %   PROPOSAL: ~table
  %   PROPOSAL: ~events
  %   PROPOSAL: ~list
  %       CON: Table implies multiple dimensions more than list.
  %            List implies ~1D (or even ~set).
  %       NOTE: solo_ns_ops.xml.html:
  %             Title "Non-Standard Operations (NSO) List"
  %   PROPOSAL: NsoList
  %
  % PROPOSAL: Terminology similar to MMS, mission/mms/mms_ns_ops.xml, i.e.
  %           event = one time interval in NSO table.
  %   PRO: MMS table more recent design than Cluster table.
  % PROPOSAL: Terminology similar to Cluster, mission/cluster/caa/caa-control/ns_ops.xml, i.e.
  %           operation = one time interval in NSO table.
  %   PRO: Cluster table has been used more than MMS table.
  %
  % PROPOSAL: Ability to add time margins for selected QRCIDs.
  %   PRO: Easier to change pre-existing margins than to modify the XML file.
  %       CON: Harder for outsiders to interpret (& edit) the XML file.
  %   Ex: thruster_firings.
  %
  % TODO-DEC: Where put requirements/assertions on NSO table content?
  %   PROPOSAL: Distinguish between assertions which are specific for
  %       (1) BICAS (not tests),
  %       (2) file format,
  %       (3) NsoTable
  %   NOTE: Assertions:
  %       t_start <= t_stop           : NsoTable
  %       t_start sorted globally.    : NsoTable? File format?
  %       Events do not overlap.      : NsoTable
  %       Only using official QRCIDs. : BICAS
  %   PROPOSAL:
  %       (1) Strictly necessary assertions in constructor.
  %       (2) Assertions on file content only in read_file_raw(), but after
  %           invoking constructor.
  %       (3) No assertions (besides valid XML file) in read_file_raw().
  %
  % PROPOSAL: Change NSO table XML format: rcsNsoId --> Nsoid
  %   PRO: Consistent with code and readme.txt.
  %   PROBLEM: Must update BICAS version simultaneously as NSO table file,
  %            also on brain. Can not update NSO table file only with e.g. new
  %            thruster firings.



  properties(SetAccess=immutable, GetAccess=public)
    % See constructor.
    % NOTE: All fields are Nx1 vectors.

    evtStartTt2000Array
    evtStopTt2000Array
    evtQrcidCa
  end



  %#####################
  %#####################
  methods(Access=public)
    %#####################
    %#####################


    % ARGUMENTS
    % =========
    % evtStartTt2000Array
    %       Column array of timestamps that represent the beginning of
    %       events.
    %       NOTE: Must increment.
    % evtStopTt2000Array
    %       Column array of timestamps that represent the end of events.
    % evtQrcidCa
    %       Column cell array of QRCIDs.
    %       NOTE: Same RCS QRCID may occur multiple times. Not unique.
    function obj = NsoTable(evtStartTt2000Array, evtStopTt2000Array, evtQrcidCa)

      %============
      % ASSERTIONS
      %============
      % PROPOSAL: Move ~all assertions to bicas.NsoTable.read_file_BICAS ?
      % PROPOSAL: Collect ~all assertions, in constructor (here) and in
      %           bicas.NsoTable.read_file_BICAS ?
      % PROPOSAL: Check that FULL_SATURATION and PARTIAL SATURATION do not overlap.
      irf.assert.sizes(...
        evtStartTt2000Array, [-1], ...
        evtStopTt2000Array,  [-1], ...
        evtQrcidCa,          [-1]);

      assert(isa(evtStartTt2000Array, 'int64'))
      assert(isa(evtStopTt2000Array,  'int64'))
      assert(isa(evtQrcidCa,          'cell' ))

      % ASSERTION: All events have non-negative length.
      assert(all(evtStartTt2000Array <= evtStopTt2000Array), ...
        'Not all events have non-negative length.')

      %--------------------------------------------------
      % ASSERTION: Event start times are sorted globally
      %--------------------------------------------------
      % IMPLEMENTATION NOTE: Can not assume that both start & stop
      % timestaps are sorted. One event may entirely contain another
      % event (with different QRCID) in time. Therefore enforcing only
      % sorted start values, but not sorted stop values.
      % IMPLEMENTATION NOTE: Can not assume "strictly ascending" values,
      % since events with separate QRCIDs may begin at the exact same
      % instant.
      if ~issorted(evtStartTt2000Array)
        iEvt = find(diff(evtStartTt2000Array) < 0) + 1;
        assert(~isempty(iEvt));

        utcCa = irf.cdf.TT2000_to_UTC_str_many(...
          evtStartTt2000Array(iEvt), 9);

        sCa = irf.str.sprintf_many('    %s\n', utcCa);
        timestampsListStr = strjoin(sCa);

        error('BICAS:FailedToReadInterpretNsOps', ...
          ['NsoTable.evtStartTt2000Array is not sorted. Events', ...
          ' beginning at the following timestamps begin earlier', ...
          ' than the precedeing events:\n%s'], ...
          timestampsListStr);
      end

      %----------------------------------------------------------------
      % ASSERTION: Events with the same QRCID do not overlap (and are
      % time sorted).
      %----------------------------------------------------------------
      uniqueEvtQrcidCa = unique(evtQrcidCa);
      for i = 1:numel(uniqueEvtQrcidCa)
        qrcid = uniqueEvtQrcidCa{i};
        b = strcmp(qrcid, evtQrcidCa);

        % NOTE: ASSUMPTION: Start timestamps are already time-sorted.
        % NOTE: Transposing before 2D-->1D vector.
        % NOTE: 'strictascend' excludes ~adjacent events.
        temp = [...
          evtStartTt2000Array(b), ...
          evtStopTt2000Array(b)]';
        tt2000Array = temp(:);
        assert(issorted(tt2000Array, 'strictascend'), ...
          ['At least two events for qrcid="%s"', ...
          ' seem to overlap with each other.'], qrcid)
      end

      % CASE: Data seems OK.

      %=====================
      % Store data in class
      %=====================
      obj.evtStartTt2000Array = evtStartTt2000Array;
      obj.evtStopTt2000Array  = evtStopTt2000Array;
      obj.evtQrcidCa          = evtQrcidCa;
    end



    % Determine which RCS QRCIDs apply to which timestamps. Given e.g. a
    % zVar Epoch, obtain lists of indices to CDF records.
    %
    %
    % ARGUMENTS
    % =========
    % tt2000Array
    %       Column array of TT2000 timestamps. Intended to be zVar Epoch.
    %
    %
    % RETURN VALUES
    % =============
    % bEvtArraysCa
    %       Cell array of arrays of logical indices into tt2000Array.
    %       {iMatchingEvent}(iTimestamp) = logical
    % evtQrcidCa
    %       List of QRCIDs. {iMatchingEvent} = qrcid
    % iGlobalEventsArray
    %       1D array of indices into NSO event list. Can be used for logging
    %       the tabulated (global) NSO events that affect e.g. a particular
    %       CDF.
    %       (iMatchingEvent) = iGlobalEvt
    %       NOTE: Useful for identifying event in NSO table, and hence the
    %             begin & end timestamps for logging.
    %
    function [bEvtArraysCa, evtQrcidCa, iGlobalEventsArray] = get_NSO_timestamps(obj, tt2000Array)
      % PROPOSAL: Sort return data by QRCID, i.e. all return values have a
      %           top-level index iQrcid.
      %   PRO: Quality algorithms should only operate on a per-sample and
      %   per-QRCID basis basis and do not need to be aware of separate.
      %   CON: Too complicated data structure.
      %       CON-PROPOSAL: Use struct array.
      %           MatchArray(iQrcid)
      %               .bMatch
      %               .qrcid
      %               .iNso
      %   PROPOSAL:
      %       bEvtArraysCa{iQrcid}(iTimestamp)
      %       evtQrcidCa{iQrcid}
      %       iGlobalEventsCa{iQrcid}(iEvt)
      %
      % PROPOSAL: Return class to simplify return values.

      assert(isa(tt2000Array, 'int64') && iscolumn(tt2000Array), ...
        'tt2000Array is not an int64 column vector.')

      if isempty(tt2000Array)
        bEvents = false(0, 1);
      else
        bEvents = irf.utils.intervals_intersect(...
          obj.evtStartTt2000Array, ...
          obj.evtStopTt2000Array, ...
          min(tt2000Array), ...
          max(tt2000Array), ...
          'closed intervals');
      end

      % =====================================
      % Assign evtQrcidCa, iGlobalEventsArray
      % =====================================
      % IMPLEMENTATION NOTE: Obtain SUBSET of NSO table EVENTS which
      % overlap with timestamps in tt2000Array. Indirectly also removes
      % irrelevant RCS QRCIDs (not just events) for the code after.
      evtStartTt2000Array = obj.evtStartTt2000Array(bEvents);
      evtStopTt2000Array  = obj.evtStopTt2000Array(bEvents);
      evtQrcidCa          = obj.evtQrcidCa(bEvents);
      iGlobalEventsArray  = find(bEvents);

      % Normalize 0x0 to 0x1
      % --------------------
      % IMPLEMENTATION NOTE: Must normalize empty vectors due to
      % inconsistent MATLAB behaviour. Otherwise column vectors become
      % non-column vectors.
      % Ex:
      %     a = [3, 4, 5]';  size(a(false(3,1)))  == [0, 1]
      %     a = [3];         size(a(false))       == [0, 0]    # NOTE!
      %     a = zeros(0, 1); size(a(false(0, 1))) == [0, 1]
      evtStartTt2000Array = evtStartTt2000Array(:);
      evtStopTt2000Array  = evtStopTt2000Array(:);
      evtQrcidCa          = evtQrcidCa(:);
      % Ex:
      %     size(find(false(0, 1))) == [0, 1]
      %     size(find(false(1, 1))) == [0, 0]    # NOTE!
      %     size(find(false(3, 1))) == [0, 1]
      iGlobalEventsArray  = iGlobalEventsArray(:);



      % ===================
      % Assign bEvtArraysCa
      % ===================
      % IMPLEMENTATION NOTE: obj.evtQrcidCa is NOT a list of unique NSO
      % IDs, but the return value "evtQrcidCa" must be a list of unique
      % QRCIDs. One must distinguish between these two.
      nEvents      = numel(evtQrcidCa);
      bEvtArraysCa = cell(nEvents, 1);
      for iEvent = 1:nEvents    % Matching events (not global).

        tt2000_1 = evtStartTt2000Array(iEvent);
        tt2000_2 = evtStopTt2000Array(iEvent);

        bMatch = (tt2000_1 <= tt2000Array) & (tt2000Array <= tt2000_2);
        bEvtArray = false(size(tt2000Array));
        bEvtArray(bMatch) = true;

        bEvtArraysCa{iEvent, 1} = bEvtArray;
      end



      % ASSERTIONS
      irf.assert.sizes(...
        bEvtArraysCa,       [-1], ...
        evtQrcidCa,         [-1], ...
        iGlobalEventsArray, [-1])
    end    % get_NSO_timestamps



  end    % methods(Access = public)



  %#############################
  %#############################
  methods(Static, Access=public)
    %#############################
    %#############################



    % Read SolO non-standard operations (NSO) XML file for *BICAS* and
    % return the content as an instance of bicas.NsoTable.
    function NsoTable = read_file_BICAS(filePath)
      [evtStartTt2000Array, evtStopTt2000Array, evtQrcidCa] = ...
        bicas.NsoTable.read_file_raw(filePath);

      % ASSERTION: No non-BICAS QRCIDs
      % ------------------------------
      % List of all legal QRCIDs.
      LEGAL_QRCID_CA = struct2cell(bicas.const.QRCID);
      irf.assert.castring_set(LEGAL_QRCID_CA)
      illegalEvtQrcidSet = setdiff(evtQrcidCa, LEGAL_QRCID_CA);
      assert(isempty(illegalEvtQrcidSet), ...
        'NSO table file contains illegal QRCID(s): %s.',  ...
        ['"', strjoin(illegalEvtQrcidSet, '", "'), '"'])

      NsoTable = bicas.NsoTable(...
        evtStartTt2000Array, evtStopTt2000Array, evtQrcidCa);
    end



    % Read SolO non-standard operations (NSO) XML file and return "raw
    % content" (without all checks) as variables.
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
    % Same fields as in class bicas.NsoTable.
    % evtStartTt2000Array
    % evtStopTt2000Array
    % evtQrcidCa
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    % First created 2020-09-21.
    %
    function [evtStartTt2000Array, evtStopTt2000Array, evtQrcidCa] = ...
        read_file_raw(filePath)

      RootXmlElem      = xmlread(filePath);
      MainXmlElem      = bicas.NsoTable.getXmlUniqChildElem(RootXmlElem, 'main');
      TablesXmlElem    = bicas.NsoTable.getXmlUniqChildElem(MainXmlElem, 'eventsTable');
      EventXmlElemList = TablesXmlElem.getElementsByTagName('event');

      nEvents = EventXmlElemList.getLength;

      evtStartTt2000Array = int64(zeros(nEvents, 1));
      evtStopTt2000Array  = int64(zeros(nEvents, 1));
      evtQrcidCa          = cell(nEvents, 1);

      for i = 1:nEvents
        % NOTE: Subtract by one.
        EventXmlElem = EventXmlElemList.item(i-1);

        startUtc = bicas.NsoTable.getXmlChildElemStr(EventXmlElem, 'startTimeUtc');
        stopUtc  = bicas.NsoTable.getXmlChildElemStr(EventXmlElem, 'stopTimeUtc');
        qrcid    = bicas.NsoTable.getXmlChildElemStr(EventXmlElem, 'rcsNsoId');
        % NOTE: NSO XML file contains string "rcsNsoId" which is
        % technically against the naming convention (w.r.t.
        % capitalization) used in the source code. Keeping the old
        % format in the XML file for compatibility.

        startTt2000 = spdfparsett2000(startUtc);
        stopTt2000  = spdfparsett2000(stopUtc);

        evtStartTt2000Array(i, 1) = startTt2000;
        evtStopTt2000Array(i, 1)  = stopTt2000;
        evtQrcidCa{i, 1}          = qrcid;
      end
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
      ChildXmlElem = bicas.NsoTable.getXmlUniqChildElem(XmlElem, childTagName);
      s            = bicas.NsoTable.getXmlElemStr(ChildXmlElem);
    end



  end    % methods(Static, Access=private)



end
