function out=caa_msp(varargin)
% CAA_MSP read master science plan from C?_JP_XXX
%
%	CAA_MSP(tint,cl_id,'PSE:BS') display predicted bow shock crossings
%	CAA_MSP(tint,cl_id,'PSE:BS',filename) write result to filename

%% check inputs
writeToFile = false;
if nargin >= 3 
	tint = varargin{1};
	cl_id = varargin{2};
	filterString = varargin{3};
else
	return;
end
if nargin == 4 % filename specified
  writeToFile = true;
  fileName = varargin{4};
end
%% check filter
if strfind(filterString,'PSE:') %#ok<STRIFCND>
	dataObj = ['C' num2str(cl_id) '_JP_PSE'];
	regionFilter = filterString(strfind(filterString,':')+1 : end);
else 
	return
end

%% Read data

eventCode = local.c_read(['event_code__' dataObj],tint,'caa');
eventSubCode = local.c_read(['event_sub_code__' dataObj],tint,'caa');
orbitNumber = local.c_read(['orbit_num__' dataObj],tint);
%timeEpoch = local.c_read(['time_tags__' dataObj],tint,'caa');
timeEpoch = orbitNumber(:,1);
tIso = irf_time(timeEpoch,'epoch>utc');

eventCode = eventCode.data;
eventCode = squeeze(eventCode);
eventCode(:,3:end) = [];

eventSubCode = eventSubCode.data;
eventSubCode = squeeze(eventSubCode);
eventSubCode(:,2:end) = [];

%% Filter region and display
indRegionLogical = ~any(eventCode - repmat(regionFilter,size(eventCode,1),1),2);
indRegion = find(indRegionLogical);
if writeToFile
  fid = fopen(fileName,'w');
else
  fid =[];
end
outPut(fid,'#################################')
outPut(fid,['  Data object: ' dataObj]);
outPut(fid,['Region filter: ' regionFilter]);
outPut(fid,'#################################')
for j=1:numel(indRegion)
	outPut(fid,[num2str(orbitNumber(indRegion(j),2)) ' ' tIso(indRegion(j),1:16) ' ' ...
		eventCode(indRegion(j),:) ' ' eventSubCode(indRegion(j),:)]);
end
outPut(fid,'#################################')
if writeToFile
  fclose(fid);
end

function out = outPut(fid,text)
if isempty(fid)
  disp(text)
else
  fprintf(fid,'%s\n',text);
end

