function [maneuvers, fileInterval] = mms_read_timeline(xmlFile, scIdstr)
% MMS_READ_TIMELINE get maneuver information from mms timeline xml files
%	[maneuvers, fileInfo] = MMS_READ_TIMELINE(xmlFile, scIdstr) reads the
%   FDOA xmlFile and extract maneuvers for scIdstr.
%
%   Note: This function uses XPath directly to extract the interested
%   maneuver information from the FDOA xml files and is thefore less
%   general purpose then things like XML2STRUCT. Should FDOA change their
%   file structure this function may require changes. The time required for
%   one xml file to be processed is greatly reduces compared with
%   xml2struct, (0.042 seconds as opposed to 24.93 seconds).
%
%   Input:
%     xmlFile = xml file to be read
%     scIdstr = string of which spacecraft to be extracted, values '1', 
%               '2', '3' or '4'. Or a combination such as '1234'.
%   Output:
%     manuevers    = struct containing
%       .mms1      = cell of maneuvers on MMS1, in irf.tint time interval.
%       .mms2      = corresponding on MMS2.
%       .mms3      = corresponding on MMS3.
%       .mms4      = corresponding on MMS4.
%     fileInterval = irf.tint time interval of file coverage.
%
%   Example:
%
%	[maneuvers, fileInterval] = MMS_READ_TIMELINE('/path/to/mms_timeline_2016011_2016073_v02.xml','1');
%    returns a struct maneuvers containing
%       .mms1 = Cell with time intervals of maneuvers for MMS1, in irf.tint
%    and a time interval fileInfo containing with start and stop time of
%    the xml file coverage, in irf.tint.
%
% 	See also MMS_MANEUVERS.
%

%% DO NOT USE
% THIS IS FOR TESTING ONLY, DO NOT USE ! ! ! 
%% FIXME: THIS IS FOR TESTING ONLY!!!
if(~strcmp(getComputerName,'thonilaptop')), error('DO NOT USE THIS FUNCTION'); end

% Verify input & output
narginchk(1,2);
if(~exist(xmlFile,'file')), error(['File does not exist: ', xmlFile]); end
if(exist('scIdstr','var'))
  if(isnumeric(scIdstr)), scIdstr = num2str(scIdstr); end
else
  irf.log('warning', 'scId not specified, processing all MMS 1234.');
  scIdstr='1234';
end
nargoutchk(1,2);

irf.log('debug', ['Loading timeline file: ', xmlFile,'.']);

% Parse file and put it in RAM
import javax.xml.parsers.*;
domFactory = DocumentBuilderFactory.newInstance();
builder = domFactory.newDocumentBuilder();
doc = builder.parse(xmlFile);

% Process
import javax.xml.xpath.*;
factory = XPathFactory.newInstance();
xpath = factory.newXPath();

if(nargout>=2)
  % Get some file information (coverage of the file)
  query = sprintf('timeline');
  expr = xpath.compile(query);
  res = expr.evaluate(doc, XPathConstants.NODESET);
  childNode = res.item(0).getChildNodes; % Java start with index 0
  startTimeStr = char(childNode.getElementsByTagName('starttime').item(0).getTextContent);
  stopTimeStr = char(childNode.getElementsByTagName('stoptime').item(0).getTextContent);
  fileInterval = irf.tint(convertTime(startTimeStr), convertTime(stopTimeStr));
end

for scId=1:length(scIdstr)
  irf.log('notice',['Looking for maneuvers on MMS', scIdstr(scId)]);
  % Verify it is one of 1,2,3,4
  if(~ismember(scIdstr(scId),{'1','2','3','4'}))
    errStr = ['Unexpected scIdstr ',scIdstr,'. Should be 1, 2 ,3 or 4 or a combination.'];
    irf.log('critical',errStr); error(errStr);
  end
  % Locate Manuevers for MMS Spacecraft "scId"
  query = sprintf('timeline/node[@id=''Maneuver Plans (FDOA-11)'']/node[@id=''Maneuver - Spacecraft %s'']/event', scIdstr(scId));
  expr = xpath.compile(query);
  res = expr.evaluate(doc, XPathConstants.NODESET);
  manLen = res.getLength();
  if(manLen==0)
    irf.log('warning',[xmlFile,' did not contain any maneuvers for MMS',scIdstr(scId),'.']);
  end
  maneuvers.(['mms',scIdstr(scId)]) = cell(manLen,1);
  for ii = 0:manLen-1
    % Note: Java start with index 0.
    childNode = res.item(ii).getChildNodes;
    % Get start time of maneuver
    startTimeStr = char(childNode.getElementsByTagName('starttime').item(0).getTextContent);
    % Get stop time of maneuver
    stopTimeStr = char(childNode.getElementsByTagName('stoptime').item(0).getTextContent);
    % Tint interval
    maneuvers.(['mms',scIdstr(scId)]){ii+1,1} = irf.tint(convertTime(startTimeStr), convertTime(stopTimeStr));
  end
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
