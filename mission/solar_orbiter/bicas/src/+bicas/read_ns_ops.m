%
% Read SolO non-standard operations (NSO) XML file.
%
%
% ARGUMENTS
% =========
% filePath
%
%
% RETURN VALUES
% =============
% NsoTable : Struct of arrays representing file content.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-09-21.
%
function NsoTable = read_ns_ops(filePath)
    % TODO-DEC: Time format? String? TT2000? Numeric?
    %   NOTE: Want to assert t1<t2.
    % PROPOSAL: Permit multiple forms of XML time input: (t1, t2), (t1,dt), (dt, t2)
    %
    % PROPOSAL: Class for the data.
    %   PROPOSAL: Methods for finding rcsNsoCodes for Epoch values.
    %       PROPOSAL: Submit array of Epoch values for entire dataset for speed.
    %           PROPOSAL: Return indices to Epoch. [iEpochCa{iRcsCode}(iEpoch), rcsNsoCodeCa]
    %               PRO: Can iterate over RCS codes instead of Epoch values.

    RootElem      = xmlread(filePath);
    TablesElem    = getUniqChildElem(RootElem, 'table');
    EventElemList = TablesElem.getElementsByTagName('event');
    
    nEvents = EventElemList.getLength;
    %NsoArray = EJ_library.utils.empty_struct([nEvents, 1], 'startTt2000', 'stopTt2000', 'rcsNsoCode');
    
    startTt2000Array = int64(zeros(nEvents, 1));
    stopTt2000Array  = int64(zeros(nEvents, 1));
    rcsNsoCodeCa     = cell(nEvents, 1);
    
    for i = 1:nEvents
        EventElem = EventElemList.item(i-1);    % NOTE: Subtract by one.

        startUtc   = getChildElemStr(EventElem, 'startTimeUtc');
        stopUtc    = getChildElemStr(EventElem, 'stopTimeUtc');
        rcsNsoCode = getChildElemStr(EventElem, 'rcsNsoCode');

        % DEBUG
        %fprintf('%i: %s -- %s  %s\n', iEvent, startUtc, stopUtc, rcsNsoCode);
    
        % Converting UTC strings to numerical format so that assertion can be
        % used on them.
        startTt2000 = spdfparsett2000(startUtc);
        stopTt2000  = spdfparsett2000(stopUtc);
        
        % ASSERTION
        % IMPLEMENTATION NOTE: This assertion requires converting the UTC
        % strings to a numerical format.
        assert(startTt2000 < stopTt2000, ...
            'BICAS:read_ns_ops:FailedToReadInterpretNsOps', ...
            'Start time does not preceed stop time for event beginning at UTC "%s".', ...
            startUtc)
        
        startTt2000Array(i, 1) = startTt2000;
        stopTt2000Array(i, 1)  = stopTt2000;
        rcsNsoCodeCa{i, 1}     = rcsNsoCode;
    end
    
    NsoTable = struct(...
        'startTt2000Array', {startTt2000Array}, ...
        'stopTt2000Array',  {stopTt2000Array}, ...
        'rcsNsoCodeCa',     {rcsNsoCodeCa});

end



% Elem : Element that has exactly one child in the form of an element with
%        specified tag name.
function ChildElem = getUniqChildElem(Elem, childTagName)
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
function s = getElemStr(Elem)
    ChildNodesList = Elem.getChildNodes();
    assert(ChildNodesList.getLength == 1, ...
        'BICAS:read_ns_ops:FailedToReadInterpretNsOps', ...
        'XML element does not have exactly one child node as expected.')
    
    s = char(ChildNodesList.item(0).getTextContent);
end



function s = getChildElemStr(Elem, childTagName)
    ChildElem = getUniqChildElem(Elem, childTagName);
    s         = getElemStr(ChildElem);
end
