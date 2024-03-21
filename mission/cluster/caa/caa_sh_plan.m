function caa_sh_plan(yyyy,mm)
%CAA_SH_PLAN  identify magnetosheath/sw intervals
%
% caa_sh_plan(yyyy,[month_list])
%
% Example:
%     caa_sh_plan(2012,1:5) will generate MP crossings for Jan-May 2012
%
% See also CAA_FIND_MP

if nargin <2
  if yyyy==2001, mm=2:12;
  else, mm = 1:12;
  end
end

c_ctl('set',5,'isdat_db','db:9');
force_orbit_read=[0 0 0 0];
force_orbit_splitting=[0 0 0 0];
force_MP_determination=[0 0 0 0];
solar_wind_data='omni2'; %'ace'

for cl_id = 1:4
  % Fetch position
  v_s = sprintf('R%dY%d',cl_id,yyyy);
  cl_id_s = num2str(cl_id);
  if (~force_orbit_read(cl_id)), if exist('./mR.mat','file'), eval(['load ./mR.mat ' v_s]), end, end
  if ~exist(v_s,'var')
    et = []; R = [];
    for mo=mm

      irf_log('proc',sprintf('Getting R%d for %d-%d',cl_id,yyyy,mo))
      if ~isempty(et), st = et;
      else, st = toepoch([yyyy mo 01 00 00 00]);
      end

      if mo==12, et = toepoch([yyyy+1 01 01 00 00 00]);
      else, et = toepoch([yyyy mo+1 01 00 00 00]);
      end

      data = getData(ClusterDB, st, et-st, cl_id, 'r', 'nosave');
      if isempty(data) || (data{2}(end,1)-data{2}(1,1)) < (et-st-3600*4)
        data=[];
        irf_log('proc',['Error in position data for month ' num2str(mo) '. Fetching one day at a time.'])
        for st1=st:86400:et-86399
          data2=getData(ClusterDB, st1, 86400, cl_id, 'r', 'nosave');
          if ~isempty(data2)
            if isempty(data),data=data2;
            else, data{2}=[data{2}' data2{2}']'; end
          end
        end
      end
      if isempty(data), error('cannot fetch position'), end

      R_tmp = irf_abs(data{2});
      R = [R; R_tmp];
      clear R_tmp
    end

    irf_log('save',[v_s ' -> mR.mat'])
    if exist('./mR.mat','file'), eval([v_s '=R; save ./mR.mat ' v_s ' -append'])
    else, eval([v_s '=R; save ./mR.mat ' v_s])
    end
  else
    eval([ 'R=' v_s ';'])
  end

  % Find the perigee and split into orbits
  v_s = sprintf('ORB%dY%d',cl_id,yyyy);
  if (~force_orbit_splitting(cl_id)), if exist('./mPlan.mat','file'), eval(['load ./mPlan.mat ' v_s]), end, end
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
          irf_log('proc',['C' cl_id_s ' perigee at ' epoch2iso(t_prev,1)])
        end
        if ~isempty(ORB)
          dt_expected=ORB(1,2);
          if(abs(R(ii(o),1)-t_prev-dt_expected) < 5000)
            ORB = [ORB; [R(ii(o),1)  R(ii(o),1)-t_prev]];
            irf_log('proc',['C' cl_id_s ' perigee at ' epoch2iso(R(ii(o),1),1)])
          else
            n_skipped=round((R(ii(o+1),1)-R(ii(o-1),1))/dt_expected);
            for i=1:n_skipped, ORB = [ORB; [ORB(end,1)+dt_expected  dt_expected]]; end
            R(ii(o),1)=ORB(end,1);
            irf_log('proc',['Fixing '  num2str(n_skipped) ' orbits on C' cl_id_s ' near ' epoch2iso(R(ii(o),1),1)])
          end
        end
        t_prev = R(ii(o),1);
      end
    end
    if isempty(ORB), irf_log('proc',['no data for ' v_s]), continue, end
    irf_log('save',[v_s ' -> mPlan.mat'])
    if exist('./mPlan.mat','file')
      eval([v_s '=ORB; save ./mPlan.mat ' v_s ' -append'])
    else, eval([v_s '=ORB; save ./mPlan.mat ' v_s])
    end
  else, eval([ 'ORB=' v_s ';'])
  end

  % Find magnetopause crossings
  v_s = sprintf('MP%dY%d',cl_id,yyyy);
  if (~force_MP_determination(cl_id)), if exist('./mPlan.mat','file'), eval(['load ./mPlan.mat ' v_s]), end, end
  if ~exist(v_s,'var')
    MP = [];
    for o=1:length(ORB)
      [t_out, t_in] = caa_find_mp(ORB(o,1),ORB(o,2),cl_id, R,solar_wind_data);
      if ~isempty(t_out) && ~isempty(t_in)
        if t_in>t_out, MP = [MP; [t_out, t_in]];
        else
          irf_log('proc',['ERROR: inbound (' epoch2iso(t_in,1)...
            ') is before the outboud (' epoch2iso(t_out,1) ')'])
        end
      end
    end

    if isempty(MP), irf_log('proc',['no data for ' v_s]), continue, end

    irf_log('save',[v_s ' -> mPlan.mat'])
    if exist('./mPlan.mat','file')
      eval([v_s '=MP; save ./mPlan.mat ' v_s ' -append'])
    else, eval([v_s '=MP; save ./mPlan.mat ' v_s])
    end
  else, eval([ 'MP=' v_s ';'])
  end

  v_s = sprintf('MPauseY%d',yyyy);
  if exist('./mPlan.mat','file'), eval(['load ./mPlan.mat ' v_s]), end
  if ~exist(v_s,'var'), MP3h = [];
  else, eval([ 'MP3h=' v_s ';'])
  end

  % Plot
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
  irf_zoom('x',[toepoch([yyyy 1 1 0 0 0]) toepoch([yyyy+1 1 1 0 0 0])]);
  ylabel('time [hours] from perigee')
  title(sprintf('Cluster %d',cl_id))
  clear MP R
end
