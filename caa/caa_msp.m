function out=caa_msp(varargin)
% CAA_MSP read master science plan from C?_JP_XXX
%
%	out = CAA_MSP(tint,cl_id,'PSE:BS') get predicted bow shock crossings

%% check inputs
if nargin == 3, 
	tint = varargin{1};
	cl_id = varargin{2};
	filterString = varargin{3};
else
	return;
end

%% check filter
if strfind(filterString,'PSE:')
	dataObj = [C num2str(cl_id) '_JP_PSE'];
	regionFilter = filterString(strfind(filterString,':')+1 : end);
else 
	return
end

%% Read data

eventCode = local.c_read(['event_code__' dataObj],tint,'caa');
eventSubCode = local.c_read(['event_sub_code__' dataObj],tint,'caa');
timeEpoch = local.c_read(['time_tags__' dataObj],tint);
tIso = irf_time(timeEpoch,'iso');

eventCode = eventCode.data;
eventCode = squeeze(eventCode);
eventCode(:,3:end) = [];

eventSubCode = eventSubCode.data;
eventSubCode = squeeze(eventSubCode);
eventSubCode(:,2:end) = [];

%% Filter region and display
indRegion = ~any(eventCode - repmat(regionFilter,size(eventCode,1),1),2);

disp('#################################')
disp(['  Data object: ' dataObj]);
disp(['Region filter: ' regionFilter]);
disp('#################################')
for j=1:numel(indRegion)
	disp([tIso(indRegion(j),:) ' ' ...
		eventCode(indRegion(j)) ' ' eventSubCode(indRegion(j))]);
end
disp('#################################')
