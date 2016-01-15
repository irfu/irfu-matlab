function maneuvers = mms_maneuvers( Tint, fileName )
% Read Timeline and return maneuvers in the interval Tint
%
% TESTING: DO NOT USE!

maneuvers = [];

% List timeline files
% /data/mms/ancillary/mms/timeline/mms_timeline_yyyyDOY_yyyyDOY_vXX.xml

% In mms_db, filename gives starttime and stoptime.

% Use the last file matching start of Tint in order to avoid predicted
% maneuvers as much as possible..

% Load file to struct using 
inp = xml2struct(fileName); % Note takes about 24 seconds per file.

% Shorten struct, or return if unexpected format...
if(isstruct(inp) && isfield(inp,'timeline')), inp = inp.timeline; else return, end
if(~isstruct(inp) || ~isfield(inp,'node')), return, end

% Loop ii=1:size(inp.node) to check which "node" (cell)
% corresponds to having its "node{ii}.Attributes.id==Maneuver Plans (FDOA-11)".
for ii=1:size(inp.node)
  % If it is not Maneuver, move to the next.
  if(~strcmpi(inp.node{ii}.Attributes.id, 'Maneuver Plans (FDOA-11)')), continue, end

  % In one example (mms_timeline_2016014_2016073_v02.xml) this first node is
  % the correct one.
  % Next "node" level appears to correspond to S/C.
  % Ie. input.timeline.node{1}.node{4}, MMS_Maneuvers for MMS4.
  % with input.timeline.node{1}.node{4}.event{1,1}.starttime.Text
  % Start time for Maneuver on MMS4 in format yyyy-doy/HH:MM:SS
  % And input.timeline.node{1}.node{4}.event{1,1}.stoptime.Text
  % being corresponding stoptime in same format.

  % Next maneuver is found in test1.timeline.node{1}.node{4}.event{1,2} and
  % so on...
end


end
