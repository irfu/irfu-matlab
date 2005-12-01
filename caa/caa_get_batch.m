function caa_get_batch(iso_t,dt,sdir)
%CAA_GET_BATCH  CAA wrapper for c_get_batch
%
% caa_get_batch(iso_t,dt,sdir)
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

REQ_INT = 60; % intervals (sec) for which we request FDM
DB_S = c_ctl(0,'isdat_db');
DP_S = c_ctl(0,'data_path');
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
	c_eval('tm? = tm;')
end
clear tm tm_cur

% Make SC_LIST from SC for which TM is not empty 
sc_list = [];
c_eval('if ~isempty(tm?), sc_list = [sc_list ?]; end')

if isempty(sc_list), irf_log('dsrc','No data'), return, end

c_eval('tm = tm?;',sc_list(1))

% Now we check if all TMs are the same
% TODO: handle the situation of different TMs
if length(sc_list)>1
	for j=2:length(sc_list)
		c_eval('tm_cur = tm?;',sc_list(j))
		if ~tm==tm_cur
			irf_log('dsrc',...
				'TM intervals are different for different SC. quitting')
			return
		end
	end
end

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

% Save mTMode and run c_get_batch on each interval
for j=1:size(tm,1)
	st_tmp = tm(j,1);
	sp = [sdir '/' irf_fname(st_tmp)];
	c_eval('mTMode?=tm(j,2);',sc_list)
	
	old_pwd = pwd;
	%Create the storage directory if it does not exist
	if ~exist(sp, 'dir')
		[SUCCESS,MESSAGE,MESSAGEID] = mkdir(sp);
		if SUCCESS, irf_log('save',['Created storage directory ' sp])
		else, error(MESSAGE)
		end
	end
	
	cd(sp)
	for cl_id=sc_list
		if exist('./mEFWR.mat','file')
			eval(irf_ssub('save -append mEFWR mTMode?;',cl_id))
		else, eval(irf_ssub('save mEFWR mTMode?;',cl_id))
		end
	end
	cd(old_pwd)
	
	cdb = ClusterDB(DB_S, DP_S, sp);
	if j==size(tm,1), dt_tmp = st + dt - st_tmp;
	else, dt_tmp = tm(j+1,1) - st_tmp;
	end
	c_get_batch(st_tmp,dt_tmp,'cdb',cdb)
end
