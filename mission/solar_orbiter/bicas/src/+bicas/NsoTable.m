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
    evtQrcidAr
  end
  properties(Dependent)
    nEvents
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods



    function nRows = get.nEvents(obj)
      nRows = numel(obj.evtStartTt2000Array);
    end



  end
  methods(Access=public)



    % Constructor for instantiating object from in-memory variables (i.e. not
    % file).
    % NOTE: There is a separate function for instantiating from file.
    %
    % ARGUMENTS
    % =========
    % evtStartTt2000Array
    %       Column array of timestamps that represent the beginning of
    %       events.
    %       NOTE: Must increment.
    % evtStopTt2000Array
    %       Column array of timestamps that represent the end of events.
    % evtQrcidAr
    %       Column array of QRCIDs for events.
    %       NOTE: The same QRCID may occur multiple times. Not unique.
    %
    function obj = NsoTable(evtStartTt2000Array, evtStopTt2000Array, evtQrcidAr)

      %============
      % ASSERTIONS
      %============
      % PROPOSAL: Collect ~all assertions, in constructor (here) and in
      %           bicas.NsoTable.read_file_validated ?
      % PROPOSAL: Move ~all file format assertions to
      %           bicas.NsoTable.read_file_validated()?
      irf.assert.sizes(...
        evtStartTt2000Array, [-1], ...
        evtStopTt2000Array,  [-1], ...
        evtQrcidAr,          [-1]);

      assert(isa(evtStartTt2000Array, 'int64'))
      assert(isa(evtStopTt2000Array,  'int64'))
      assert(isa(evtQrcidAr,          'string'))

      % ASSERTION: All events have non-negative length.
      assert(all(evtStartTt2000Array <= evtStopTt2000Array), ...
        'Not all events have non-negative length.')

      %--------------------------------------------------
      % ASSERTION: Event start times are sorted globally
      %--------------------------------------------------
      % IMPLEMENTATION NOTE: One can not assume that both start & stop
      % timestaps are sorted. One event may entirely contain another
      % event (with different QRCID) in time. Therefore enforcing only
      % sorted start values, but not sorted stop values.
      % IMPLEMENTATION NOTE: One can not assume "strictly ascending" values,
      % since events with separate QRCIDs may begin at the exact same
      % instant.
      if ~issorted(evtStartTt2000Array)
        % IMPLEMENTATION NOTE: Locating and printing the illegal entries in a
        % proper human-readable error message, since they could otherwise be
        % hard to manually locate and fix.
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

      % CASE: Start timestamps are time-sorted.

      %------------------------------------------------------
      % ASSERTION: Events with the same QRCID do not overlap
      %------------------------------------------------------
      uniqueEvtQrcidAr = unique(evtQrcidAr);
      for i = 1:numel(uniqueEvtQrcidAr)
        qrcid = uniqueEvtQrcidAr{i};
        b = strcmp(qrcid, evtQrcidAr);

        % NOTE: ASSUMPTION: Start timestamps are time-sorted.
        % NOTE: Transposing before 2D-->1D vector.
        % NOTE: 'strictascend' excludes ~adjacent events.
        tt2000Array = [...
          evtStartTt2000Array(b), ...
          evtStopTt2000Array(b)]';
        tt2000Array = tt2000Array(:);
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
      obj.evtQrcidAr          = evtQrcidAr;
    end



    % Determine the timestamps to which NSO events apply.
    %
    %
    % ARGUMENTS
    % =========
    % tt2000Ar
    %       Column array of TT2000 timestamps. Intended to be ZV Epoch.
    %
    %
    % RETURN VALUES
    % =============
    % NsoEventMatchAr
    %       Column array of bicas.NsoEventMatch, one per NSO event which
    %       overlaps with at least one of the timestamps in tt2000Array.
    %
    function NsoEventMatchAr = get_NSO_event_matches(obj, tt2000Ar)
      assert(isa(tt2000Ar, 'int64') && iscolumn(tt2000Ar), ...
        'tt2000Ar is not an int64 column vector.')

      NsoEventMatchAr = bicas.NsoEventMatch.empty(0, 1);
      for iEvent = 1:obj.nEvents
        tt2000start = obj.evtStartTt2000Array(iEvent);
        tt2000stop  = obj.evtStopTt2000Array(iEvent);

        qrbcAr = (tt2000start <= tt2000Ar) & (tt2000Ar <= tt2000stop);
        if any(qrbcAr)
          NsoEventMatchAr(end+1, 1) = bicas.NsoEventMatch(...
            qrbcAr, obj.evtQrcidAr(iEvent), iEvent);
        end
      end
    end    % get_NSO_event_matches



  end    % methods(Access = public)



  %#############################
  %#############################
  methods(Static, Access=public)
    %#############################
    %#############################



    % Read SolO non-standard operations (NSO) XML file for *BICAS* and
    % return the content as an instance of bicas.NsoTable.
    %
    % NOTE: The caller must supply list of legal QRCIDs to make it possible to
    % assert that the file is valid, without making this class less generic.
    %
    function NsoTable = read_file_validated(filePath, legalQrcidAr)

      [evtStartTt2000Array, evtStopTt2000Array, evtQrcidAr] = ...
        bicas.NsoTable.read_file_raw(filePath);

      % ASSERTION: No illegal QRCIDs (as specified in argument)
      % -------------------------------------------------------
      illegalEvtQrcidSet = setdiff(evtQrcidAr, legalQrcidAr);
      assert(isempty(illegalEvtQrcidSet), ...
        'NSO table file contains illegal QRCID(s): %s.',  ...
        ['"', strjoin(illegalEvtQrcidSet, '", "'), '"'])

      NsoTable = bicas.NsoTable(...
        evtStartTt2000Array, evtStopTt2000Array, evtQrcidAr);
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
    % evtQrcidAr
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    % First created 2020-09-21.
    %
    function [evtStartTt2000Array, evtStopTt2000Array, evtQrcidAr] = ...
        read_file_raw(filePath)

      RootXmlElem      = xmlread(filePath);
      MainXmlElem      = bicas.NsoTable.getXmlUniqChildElem(RootXmlElem, 'main');
      TablesXmlElem    = bicas.NsoTable.getXmlUniqChildElem(MainXmlElem, 'eventsTable');
      EventXmlElemList = TablesXmlElem.getElementsByTagName('event');

      nEvents = EventXmlElemList.getLength;

      evtStartTt2000Array = int64(zeros(nEvents, 1));
      evtStopTt2000Array  = int64(zeros(nEvents, 1));
      evtQrcidAr          = strings(nEvents, 1);

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
        evtQrcidAr{i, 1}          = qrcid;
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
