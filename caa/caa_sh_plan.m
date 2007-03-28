function caa_sh_plan(yyyy)
%CAA_SH_PLAN  identify magnetosheath/sw intervals
%
% caa_sh_plan(yyyy)
%
% See also CAA_FIND_MP
%
% $Id$

if yyyy==2001, mm=2:12;
else mm = 1:12;
end

for cl_id = 1:4
	v_s = sprintf('R%dY%d',cl_id,yyyy);
    cl_id_s = num2str(cl_id);
	if exist('./mR.mat','file'), eval(['load ./mR.mat ' v_s]), end
	if ~exist(v_s,'var')
		et = []; R = [];
		for mo=mm
		
			irf_log('proc',sprintf('Getting R%d for %d-%d',cl_id,yyyy,mo))
			if ~isempty(et), st = et;
			else st = toepoch([yyyy mo 01 00 00 00]);
			end
			
			if mo==12, et = toepoch([yyyy+1 01 01 00 00 00]);
            else et = toepoch([yyyy mo+1 01 00 00 00]);
			end
			
			data = getData(ClusterDB, st, et-st, cl_id, 'r', 'nosave');
			if isempty(data), error('cannot fetch position'), end

			R_tmp = irf_abs(data{2});
			R = [R; R_tmp];
			clear R_tmp
		end
		
		irf_log('save',[v_s ' -> mR.mat'])
		if exist('./mR.mat','file'), eval([v_s '=R; save ./mR.mat ' v_s ' -append'])
		else eval([v_s '=R; save ./mR.mat ' v_s])
		end
	else
		eval([ 'R=' v_s ';'])
	end
	
	v_s = sprintf('ORB%dY%d',cl_id,yyyy);
	if exist('./mPlan.mat','file'), eval(['load ./mPlan.mat ' v_s]), end
	if ~exist(v_s,'var')
		dR = diff(R(:,5))./diff(R(:,1));
		ii = find(dR(1:end-1,1).*dR(2:end,1)<0) +1; %shift to hit the extremum
		ORB = []; t_prev = [];
		for o=1:length(ii)
			if R(ii(o)+1,5)>R(ii(o),5)
				%Perigy
				if isempty(ORB) && ~isempty(t_prev)
					ORB = [2*t_prev-R(ii(o),1) R(ii(o),1)-t_prev];
					irf_log('proc',['C' cl_id_s ' first perigy at ' epoch2iso(2*t_prev -R(ii(o),1),1)])
					ORB = [ORB; [t_prev  R(ii(o),1)-t_prev]];
					irf_log('proc',['C' cl_id_s ' perigy at ' epoch2iso(t_prev,1)])
				end
				if ~isempty(ORB)
					ORB = [ORB; [R(ii(o),1)  R(ii(o),1)-t_prev]];
					irf_log('proc',['C' cl_id_s ' perigy at ' epoch2iso(R(ii(o),1),1)])
				end
				t_prev = R(ii(o),1);
			end
		end
		if isempty(ORB), irf_log('proc',['no data for ' v_s]), continue, end
		irf_log('save',[v_s ' -> mPlan.mat'])
		if exist('./mPlan.mat','file')
			eval([v_s '=ORB; save ./mPlan.mat ' v_s ' -append'])
		else eval([v_s '=ORB; save ./mPlan.mat ' v_s])
		end
	else eval([ 'ORB=' v_s ';'])
	end
	
	v_s = sprintf('MP%dY%d',cl_id,yyyy);
	if exist('./mPlan.mat','file'), eval(['load ./mPlan.mat ' v_s]), end
	if ~exist(v_s,'var')
		MP = [];
		for o=1:length(ORB)
			[t_out, t_in] = caa_find_mp(ORB(o,1),ORB(o,2),cl_id, R);
			if ~isempty(t_out) && ~isempty(t_in)
				if isempty(MP), MP = [t_out, t_in];
				else MP = [MP; [t_out, t_in]];
				end
			end
		end
		
		if isempty(MP), irf_log('proc',['no data for ' v_s]), continue, end
		
		irf_log('save',[v_s ' -> mPlan.mat'])
		if exist('./mPlan.mat','file')
			eval([v_s '=MP; save ./mPlan.mat ' v_s ' -append'])
		else eval([v_s '=MP; save ./mPlan.mat ' v_s])
		end
	else eval([ 'MP=' v_s ';'])
	end
	
	v_s = sprintf('MPauseY%d',yyyy);
	if exist('./mPlan.mat','file'), eval(['load ./mPlan.mat ' v_s]), end
	if ~exist(v_s,'var'), MP3h = [];
	else eval([ 'MP3h=' v_s ';'])
	end
	
	figure(cl_id), clf
	for o=1:length(ORB)
		irf_plot([ORB(o,1) ORB(o,1)+ORB(o,2); 0 ORB(o,2)/3600]'), hold on
		if ~isempty(MP3h)
			ii = find(MP3h(:,1)>=ORB(o,1) & MP3h(:,1)<=ORB(o,1)+ORB(o,2));
			if ~isempty(MP3h(ii,:))
				irf_plot([MP3h(ii,:);  (MP3h(ii,:)-ORB(o,1))/3600]','gx-')
			end
		end
		ii = find(MP(:,1)>=ORB(o,1) & MP(:,1)<=ORB(o,1)+ORB(o,2));
		if ~isempty(MP(ii,:))
			irf_plot([MP(ii,:);  (MP(ii,:)-ORB(o,1))/3600]','ro-')
		end
	end
	irf_zoom([toepoch([yyyy 1 1 0 0 0]) toepoch([yyyy+1 1 1 0 0 0])],'x')
	ylabel('time [hours] from perigy')
	title(sprintf('Cluster %d',cl_id))
	clear MP R
end
