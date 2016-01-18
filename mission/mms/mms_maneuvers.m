function maneuvers = mms_maneuvers( Tint )
% Read Timeline and return maneuvers in the interval Tint
%
% TESTING: DO NOT USE!

narginchk(1,1);
maneuvers = struct('mms1',[],'mms2',[],'mms3',[],'mms4',[]);

if(~isa(Tint,'GenericTimeArray') || ~all(size(Tint)==[2,1])), error('Unexpected Tint'); end

% List timeline files (can cover up to over 100 days...) The following is a
% quick fix, should be integrated into mms db list_files.
% /data/mms/ancillary/mms/timeline/mms_timeline_yyyyDOY_yyyyDOY_vXX.xml
% In mms_db, filename gives starttime and stoptime.
dataPath = [filesep, 'data', filesep, 'mms', filesep, 'ancillary', filesep,...
  'mms', filesep, 'timeline', filesep];
for jj=0:100 % Begin looking the year and day of Tint start
  doy = irf_time(Tint.start.ttns - 86400*jj*10^9, 'ttns>doy');
  list = dir([dataPath, 'mms_timeline_',num2str(doy(1)),num2str(doy(2)),'_*v*.xml']);
  if(~isempty(list)), break; end
end

if(isempty(list)), return, end % No results
fileName = [dataPath, list(end).name]; % Last version to match.

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
  for jj=1:length(inp.node{ii}.node)
    scID = inp.node{ii}.node{jj}.Attributes.id(end);
    if(~ismember(str2double(scID), [1, 2, 3, 4])), error('UNEXPECTED SC'); end
    if(jj>4), warning('Something else unexpected..'); end
    maneuvers.(['mms',scID]) = inp.node{ii}.node{jj};
    % For each event (maneuver)
    for kk=1:length(inp.node{ii}.node{jj}.event)
      % Convert time to TT2000.
      maneuvers.(['mms',scID]).event{kk}.starttime = convertToTT2000(maneuvers.(['mms',scID]).event{kk}.starttime.Text);
      maneuvers.(['mms',scID]).event{kk}.stoptime =  convertToTT2000(maneuvers.(['mms',scID]).event{kk}.stoptime.Text);
    end
  end

end

end

function epochTT = convertToTT2000(timeStr)
  % Convert timeStr ("yyyy-doy/hh:mm:ss") to TT2000 (int64).
  doy = sscanf(timeStr,'%4d-%3d/%2d:%2d:%2d');
  if(~all(size(doy)==[5, 1])), warning('Unexpected time'); end;
  epochTT = irf_time([doy(1), doy(2), doy(3), doy(4), doy(5), 0, 0, 0], 'doy8>ttns');
end
