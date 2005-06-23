function caa_get_batch_l1(iso_t,dt,sdir)
%CAA_GET_BATCH_L1 CAA wrapper for c_get_batch
%
% caa_get_batch_l1(iso_t,dt,sdir)
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

REQ_INT = 60; % intervals (sec) for which we request FDM
DB_S = 'disco:10';
DP_S = '/data/cluster';
SPLIT_INT = 90*60;

if ~exist(sdir,'dir'), error(['directory ' sdir ' does not exist']), end

st = iso2epoch(iso_t);

% First we check if we have any EFW HX data
% and check for NM/BM1
for cl_id=1:4
	st_tmp = st;
	tm = [];
	while st_tmp<st+dt
		[t,data] = caa_is_get(DB_S,st_tmp,REQ_INT,cl_id,'efw','FDM');
		if ~isempty(data), c_eval('tm_cur = data(5,1);',cl_id);
		else
			irf_log('dsrc',['Cannot fetch FDM for C' num2str(cl_id) ...
				' at ' epoch2iso(st_tmp,1)])
			c_eval('tm_cur = [];',cl_id)
		end
		if isempty(tm) & ~isempty(tm_cur), tm = [st_tmp tm_cur]; end
		if ~isempty(tm) & ~isempty(tm_cur)
			if tm(end,2)~=tm_cur, tm(end+1,:) = [st_tmp tm_cur]; end
		end
		st_tmp = st_tmp + REQ_INT;
	end
	% Throw away all modes>1
	disp(sprintf('THROWING AWAY %d intervals',length(find(tm(:,2)>1))))
	tm(find(tm(:,2)>1),:) = [];
	c_eval('tm? = tm;')
end
clear tm tm_cur

% Make SC_LIST from SC for which TM is not empty 
sc_list = [];
c_eval('if ~isempty(tm?), sc_list = [sc_list ?]; end')

if isempty(sc_list), irf_log('dsrc','No data'), return, end

%Create the storage directory if it does not exist
maindir = [sdir '/' irf_fname(st)];
if ~exist(maindir, 'dir')
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(maindir);
	if SUCCESS, irf_log('save',['Created storage directory ' maindir])
	else, error(MESSAGE)
	end
end

for cli=sc_list
	cdir = [maindir '/C' num2str(cli)];
	if ~exist(cdir, 'dir')
		[SUCCESS,MESSAGE,MESSAGEID] = mkdir(cdir);
		if SUCCESS, irf_log('save',['Created storage directory ' cdir])
		else, error(MESSAGE)
		end
	end
	c_eval('tm=tm?;',cli);
	
	% Split long intervals into SPLIT_INT hour chunks
	j = 1;
	while 1
		st_tmp = tm(j,1);
		tm_cur = tm(j,2);
		if j==size(tm,1), dt_tmp = st + dt - st_tmp;
		else, dt_tmp = tm(j+1,1) - st_tmp;
		end
		% For BM1 we take SPLIT_INT/3 min intervals
		if dt_tmp > SPLIT_INT*4/3*(1-tm_cur*2/3)
			if j==size(tm,1)
				tm(j+1,:) = [tm(j,1)+SPLIT_INT*(1-tm_cur*2/3) tm(j,2)];
			else
				tm(j+1:end+1,:) = [tm(j,1)+SPLIT_INT*(1-tm_cur*2/3) tm(j,2); tm(j+1:end,:)];
			end
		end
		if j>=size(tm,1), break
		else, j = j + 1;
		end
	end
	int_tmp = [];
	for inter=1:size(tm,1)
		t1 = tm(inter,1);
		if inter==size(tm,1), dt1 = st +dt -t1;
		else, dt1 = tm(inter+1,1) -t1;
		end
		if inter~=1 & inter~=size(tm,1) & dt1<500
			disp(['SKIPPING: ' epoch2iso(t1) ' - ' epoch2iso(t1+dt1)])
		else, disp(['INTERVAL: ' epoch2iso(t1) ' - ' epoch2iso(t1+dt1)])
		end
		
		% Keep track of intervals we process
		if isempty(int_tmp), int_tmp = [t1 dt1];
		else, int_tmp(end+1,:) = [t1 dt1];
		end
		
		% Get data
		c_get_batch(t1,dt1,'sc_list',cli,'sdir',cdir,...
			'vars','fdm|ibias|p|e|a','noproc')
		c_get_batch(t1,dt1,'sc_list',cli,'sdir',cdir,...
			'vars','whip|sweep|bdump|badbias|probesa|p|ps|dief','nosrc')
	end
	
	if ~isempty(int_tmp)
		% Save intervals
		c_eval('INTERVALS?=int_tmp;',cli)
		if exist([cdir '/mINTER.mat'],'file')
			c_eval(['save ' cdir '/mINTER.mat INTERVALS? -append'],cli)
		else
			c_eval(['save ' cdir '/mINTER.mat INTERVALS?'],cli)
		end
		irf_log('save',irf_ssub('INTERVALS? -> mINTER',cli))
	end
end

