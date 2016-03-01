function [maneuvers, fileInfo] = mms_read_timeline(xmlFile, scId)
% MMS_READ_TIMELINE get maneuver information from mms timeline xml files
%	[maneuvers, fileInfo] = MMS_READ_TIMELINE(xmlFile, scId) reads the FDOA
%	xmlFile and extract maneuvers for scId.
%
%   Note: This function uses XPath directly to extract the interested
%   maneuver information from the FDOA xml files and is thefore less
%   general purpose then things like XML2STRUCT. Should FDOA change their
%   file structure this function may require changes. The time required for
%   one xml file to be processed is greatly reduces compared with
%   xml2struct, (0.042 seconds as opposed to 24.93 seconds).
%
%   Input:
%     xmlFile = string of the xml file to be read
%     scId    = string of which spacecraft to be extracted, values '1', '2', '3' or '4'.
%   Output:
%     manuevers    = struct containing
%       .startTime = start time of each maneuver, in ttns.
%       .stopTime  = stop time of each maneuver, in ttns.
%     fileInfo     = optional outout, struct containing
%       .startTime = start time of xml file coverage, in ttns.
%       .stopTime  = stop time of xml file coverage, in ttns.
%
%   Example:
%
%	[maneuvers, fileInfo] = MMS_READ_TIMELINE('/path/to/mms_timeline_2016011_2016073_v02.xml','1');
%    returns a struct maneuvers containing
%       .startTime = Start time of maneuvers for MMS1, in ttns.
%       .stopTime  = Stop time of maneuvers for MMS1, in ttns.
%    and a struct fileInfo containing
%       .startTime = Start time of the xml file coverage, in ttns.
%       .stopTime  = Stop time of the xml file coverage, in ttns.
%
% 	See also MMS_MANEUVERS.
%

%% DO NOT USE
% THIS IS FOR TESTING ONLY, DO NOT USE ! ! ! 
%% FIXME: THIS IS FOR TESTING ONLY!!!
if(~strcmp(getComputerName,'thonilaptop')), error('DO NOT USE THIS FUNCTION'); end

% Verify input & output
narginchk(2,2);
if(~exist(xmlFile,'file')), error(['File does not exist: ', xmlFile]); end
if(isnumeric(scId)), scId = num2str(scId); end
if(~ismember(scId,{'1','2','3','4'})), error(['Unexpected scId: ', scId]); end
nargoutchk(1,2);

% Parse file and put it in RAM
import javax.xml.parsers.*;
domFactory = DocumentBuilderFactory.newInstance();
builder = domFactory.newDocumentBuilder();
doc = builder.parse(xmlFile);

% Process
import javax.xml.xpath.*;
factory = XPathFactory.newInstance();
xpath = factory.newXPath();

% Get some file information (coverage of the file)
query = sprintf('timeline');
expr = xpath.compile(query);
res = expr.evaluate(doc, XPathConstants.NODESET);
childNode = res.item(0).getChildNodes; % Java start with index 0
startTimeStr = char(childNode.getElementsByTagName('starttime').item(0).getTextContent);
fileInfo.startTime = convertTime(startTimeStr);
stopTimeStr = char(childNode.getElementsByTagName('stoptime').item(0).getTextContent);
fileInfo.stopTime = convertTime(stopTimeStr);

% Locate Manuevers for MMS Spacecraft "scId"
query = sprintf('timeline/node[@id=''Maneuver Plans (FDOA-11)'']/node[@id=''Maneuver - Spacecraft %s'']/event', scId);
expr = xpath.compile(query);
res = expr.evaluate(doc, XPathConstants.NODESET);

manLen = res.getLength();

maneuvers = struct('startTime', zeros(manLen, 1,'int64'), ...
  'stopTime', zeros(manLen, 1,'int64'));

for ii = 0:manLen-1
  % Note: Java start with index 0.
  childNode = res.item(ii).getChildNodes;
  % Get start time of maneuver and convert the time to ttns
  startTimeStr = char(childNode.getElementsByTagName('starttime').item(0).getTextContent);
  maneuvers.startTime(ii+1) = convertTime(startTimeStr);
  % Get stop time of maneuver
  stopTimeStr = char(childNode.getElementsByTagName('stoptime').item(0).getTextContent);
  maneuvers.stopTime(ii+1) = convertTime(stopTimeStr);
end

end

function ttns = convertTime(timeStr)
  % Convert time from format "YYYY-DOY/hh:mm:ss" (or similar) to ttns.
  doy5 = sscanf(timeStr,'%4d-%3d%c%2d:%2d:%2f');
  sec = floor(doy5(6));
  msec = floor((doy5(6)-sec)*10^3);
  usec = floor(((doy5(6) - sec)*10^3 - msec)*10^3);
  nsec = floor((((doy5(6) - sec)*10^3 - msec)*10^3 - usec) * 10^3);
  % YYYY, DOY, hh, mm, ss, msec, usec, nsec
  doy8 = [doy5(1), doy5(2), doy5(4), doy5(5), sec, msec, usec, nsec];
  ttns = irf_time(doy8, 'doy8>ttns');
end
