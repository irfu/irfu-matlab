%
% Read a SOLO HK XML file as delivered by ROC (LESIA; not ESA) to RPW teams.
%
%
% ARGUMENTS
% =========
% xmlFilePath
% namesCa
%       1D cell array. List of XML element "Name" values, in practice SolO s/c
%       HK "Mnemonic A" on
%       https://confluence-lesia.obspm.fr/display/ROC/SOLO+HK+Parameter+data
%       Ex: {'NCFT55P0'}
%       --
%       Special value: [] ==> Read data for all XML "Name" values.
%
%
% RETURN VALUES
% =============
% D
%       Struct. Each field corresponds to the identically named XML elements
%       inside the ParamSampleListElement SML elements. Each one is a column
%       cell array. All cell arrays have the same size.
%       NOTE: Since field names are taken from XML tag names, they also violate
%             MATLAB variable naming conventions.
%       .Name
%           The content between <Name> and </Name>.
%       ...
% --
% NOTE: Return format may need to be changed, in particular rawValue format
% change, or extended, if the understanding of the XML files change.
%
%
% PURPOSE AND SCOPE
% =================
% Read ONE SOLO HK XML file as delivered by ROC (LESIA; not ESA) to RPW teams.
% Code is only intended to return the data from one file, not really modify it.
% Function only returns the data that has so far been needed. More return fields
% can be added as needed.
% IMPORTANT NOTE: The function should NOT convert or filter out data (more than
% by "Name"), e.g. convert date strings or only keep certain values of "name".
% The function should read the content of the file under a minimum of
% assumptions. The caller or a wrapper function should do conversions or
% filtering if needed.
%
%
% NOTES
% =====
% Function crashes, or crashes MATLAB itself, for large files due to
% java.lang.OutOfMemoryError when using MATLAB's XML functions. One might solve
% the problem by raising MATLAB's Java Heap Size.
% Ex: 2021-09-07: solo_HK_platform_20210901_V02.xml, 144 MiB
% Ex: 2021-09-07: solo_HK_platform_20210715_V05.xml, 379 MiB
%   Java Heap=1002 MiB ==> Crashes
%   Java Heap=2000 MiB ==> MATLAB crashes
%   Java Heap=4096 MiB ==> OK
% --
% Does not sort timestamps.
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2020-12-01.
%
function D = read_RSH_file(xmlFilePath, namesCa)
% PROPOSAL: Move to irfu-matlab +solo/+shk/.
%
% TODO-DEC: Return format?
%   NEED: Easy to use if one wants to merge the content from multiple files.
%   PROPOSAL: Struct array.
%       CON: Inefficient?
%   PROPOSAL: Struct with array fields.
%   PROPOSAL: Sort data by "Name".
%       CON: Wrapper should do that?
%           CON: Since specifying namesCa, maybe not. It is a ~universal
%                operation.
%
% PROPOSAL: Only read data for specified elements (specified tag names).
%   Ex: EngineeringValue, RawValue, Unit
%   PRO: Speeds up if future implementation includes many elements.
%
% PROPOSAL: Convert raw/engineering values to numeric when the type indicates it.
%   CON-PROPOSAL: Wrapper should do that.
%   PROPOSAL: Use Type=DOUBLE or RawValueType
%       CON: No longer included in XML files.

% Variable naming convention
% PSLE = <ParamSampleListElement> ("List" is included in XML tag name)

% XML elements to always copy, except "Name".
XML_ELEMENT_TAG_WO_NAME_CA = {'EngineeringValue', 'TimeStampAsciiA'};
% XML elements to always copy.
XML_ELEMENT_TAG_CA         = [{'Name'}, XML_ELEMENT_TAG_WO_NAME_CA];

% Validate namesCa.
% Set "allNameValues".
if isnumeric(namesCa) && isempty(namesCa)
  % CASE: Special value: []
  allNameValues = true;
else
  % CASE: List of "Name" values.
  allNameValues = false;
  assert(iscell(namesCa), 'Argument "namesCa" is not a cell array.')
end

RootXmlElem = xmlread(xmlFilePath);

ResponsePartXmlElem    = getXmlUniqChildElem(RootXmlElem,         'ns2:ResponsePart');
ResponseXmlElem        = getXmlUniqChildElem(ResponsePartXmlElem, 'Response');
ParamSampleListXmlElem = getXmlUniqChildElem(ResponseXmlElem,     'ParamResponse');

PsleXmlList = ParamSampleListXmlElem.getElementsByTagName(        'ParamSampleListElement');
n           = PsleXmlList.getLength();



%=======================
% Pre-allocate D, bKeep
%=======================
D = struct();
for tagNameCa = XML_ELEMENT_TAG_CA
  D.(tagNameCa{1}) = cell(n, 1);
end
bKeep = false(n, 1);

%==============================================
% Iterate over main array/list of XML elements
%==============================================
% X = "XML index convention" (as used by functions; i.e. 0=first)
for iX = 0 : n-1
  PsleXmlElem = PsleXmlList.item(iX);

  name = getXmlUniqChildElemStr(PsleXmlElem, 'Name');

  if allNameValues || ismember(name, namesCa)

    % M = MATLAB index convention (i.e. 1=first).
    iM = iX + 1;
    bKeep(iM) = true;
    D.Name{iM, 1} = name;
    for tagNameCa = XML_ELEMENT_TAG_WO_NAME_CA
      value = getXmlUniqChildElemStr(PsleXmlElem, tagNameCa{1});
      D.(tagNameCa{1}){iM, 1} = value;
    end

  end

end

% Compress arrays: Remove indices without data.
for i = 1:numel(XML_ELEMENT_TAG_CA)
  tagName = XML_ELEMENT_TAG_CA{i};
  D.(tagName) = D.(tagName)(bKeep);

  % Normalize to coumn vector.
  % IMPLEMENTATION NOTE: Above creates 0x0 (not 0x1) cell array when
  % all(bKeep == false) which is inconsistent behaviour.
  D.(tagName) = D.(tagName)(:);
end

end



function s = getXmlUniqChildElemStr(XmlElem, childTagName)
XmlElem = getXmlUniqChildElem(XmlElem, childTagName);

s = char(XmlElem.getTextContent());
end



% XmlElem
%   Element that has exactly one child in the form of an element
%   with the specified tag name.
%
function ChildXmlElem = getXmlUniqChildElem(XmlElem, childTagName)
% NOTE: Exact same function as in bicas.NsoTable.
% PROPOSAL: Turn into separate generic function?

ChildXmlElemList = XmlElem.getElementsByTagName(childTagName);

% ASSERTION
if ~(ChildXmlElemList.getLength() == 1)
  error( ...
    ['XML element (tag name "%s") does not have exactly, ', ...
    ' one child element with tag name "%s" as expected.'], ...
    XmlElem.getNodeName(), childTagName)
end

ChildXmlElem = ChildXmlElemList.item(0);
end
