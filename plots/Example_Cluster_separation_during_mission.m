%%%%%% Cluster configuration during mission %%%%%%%
%
% can be adopted to plot also other properties of Cluster
%

%% Define flags
% mR.mat file has all positions, not needed if only apogee/perigee values needed
flag_get_all_Cluster_positions_from_isdat=0; % construct mR.mat file with Cluster position for all mission
flag_get_all_Cluster_positions_from_www=0;   % get mR.mat file from internet
% mRcluster.mat file, includes all information on perigee/apogee/separation
flag_get_Cluster_file_from_www=0; % get mRcluster.mat from www
flag_create_Cluster_file=1;       % create from mR.mat file
% figure creating/printing
flag_create_figure=1;             % create figure
flag_print_figure=1;              % print figure

%% Check if there is no mRcluster file whether it will be obtained
if ((flag_get_Cluster_file_from_www==0) && ...
    (flag_create_Cluster_file==0))
  if ~exist('mRcluster.mat','file')
    flag_get_Cluster_file_from_www=1;
  end
end

%% Defined time interval and stepping
%tint=[[2000 01 01 00 0 0 ];irf_time(now,'date>vector')]; % vector first row start time, 2nd row end time
tint=[[2000 10 01 00 0 0 ];irf_time(now,'date>vector')];
step_request=3600*100; % step with 10h, seems largest possible (why?)
step_save=3600;       % how often save points
tstart=irf_time(tint(1,:),'vector>epoch');
tend=irf_time(tint(2,:),'vector>epoch');

%% download the data using command line (works only when isdat/iscmd installed)
if flag_get_all_Cluster_positions_from_isdat
  [status,result]=system('iscmd'); %#ok<UNRCH>
  if status~=0
    disp('WARNING!!!')
    disp('You dont have isdat/iscmd installed!');
    flag_get_all_Cluster_positions_from_isdat=0;
    flag_get_all_Cluster_positions_from_www=1;
  end
end
if flag_get_all_Cluster_positions_from_www
  disp('Downloading from internet file with all Cluster orbits ...');  %#ok<UNRCH>
  if verLessThan('matlab', '8.4')
    [f,status]=urlwrite('https://space.irfu.se/~andris/cluster/mR.mat','mR.mat'); %#ok<URLWR> websave introduced in R2014b
  else
    try
      f = websave('mR.mat', 'https://space.irfu.se/~andris/cluster/mR.mat');
      status = true;
    catch
      status = false;
    end
  end
end

%% Get all Cluster position from isdat
if flag_get_all_Cluster_positions_from_isdat==1 % using command line
  R1=zeros(ceil((tend-tstart)/step_save),4);
  R2=R1;R3=R1;R4=R1;
  ii=1;icstr='1234';
  for tst=tstart:step_request:tend
    disp(irf_time(tst,'isoshort'));
    for ic=1:4
      cml=['iscmd d isdat://db:0/Cluster/' icstr(ic) '/ephemeris/position ' irf_time(tst,'epoch>utc_yyyy-mm-dd HH:MM:SS') ' ' num2str(step_request)];
      [status,result]=system(cml);
      r=textscan(result,'%f%f%f%f','CommentStyle','#','collectoutput',1);
      c_eval('r?=r{1};',ic);
    end
    c_eval('[ii1,ii2]=irf_find_comm_idx(r1,r?);r1=r1(ii1,:);r?=r?(ii2,:);',2:4);
    c_eval('[ii1,ii2]=irf_find_comm_idx(r1,r?);r1=r1(ii1,:);r?=r?(ii2,:);',2:4);
    t=r1(:,1);
    it=zeros(size(t));
    ilast=1;
    for i=1:numel(it)
      if t(i)-t(ilast)>=step_save
        ilast=i;
        it(i)=1;
      end
    end
    ind_keep=find(it==1);
    nkeep=numel(ind_keep);
    c_eval('r?(:,1)=r?(:,1)+tst; R?(ii:ii+nkeep-1,:)=r?(ind_keep,:);');
    ii=ii+nkeep;
  end
  c_eval('R?(ii:end,:)=[];');
  save mR R1 R2 R3 R4
end

%% Construct Cluster file if needed
if flag_create_Cluster_file
  %% Find apogee times for all s/c
  load mR
  izero=find(R1(:,1)==0);
  c_eval('R?(izero,:)=[];R?=irf_abs(R?);');
  c_eval('dR?=diff(R?);');
  c_eval('iap?=find(dR1(1:end-1,1) < 12*3600 & dR1(1:end-1,5) > 0 & dR1(1:end-1,5).*dR1(2:end,5) <0 );')
  c_eval('iper?=find(dR1(1:end-1,1) < 1.1*3600 & dR1(2:end,1) < 1.1*3600 & dR1(1:end-1,5) < 0 & dR1(1:end-1,5).*dR1(2:end,5) <0 );')
  c_eval('Rap?=R?(iap?+1,:);');
  c_eval('Rper?=R?(iper?+1,:);');
  save mRcluster Rap1 % initiate saving file
  c_eval('save -append mRcluster Rap? Rper?');

  %% Find separations at perigee and apogee
  c_eval('drper?!=R?(iper?,1:4)-R!(iper?,1:4);');
  c_eval('drper?!=irf_abs(drper?!);');
  c_eval('drper?!(:,1)=R?(iper?,1);');
  save -append mRcluster drper12 drper13 drper14 drper23 drper24 drper34;
  c_eval('drap?!=R?(iap?,1:4)-R!(iap?,1:4);');
  c_eval('drap?!=irf_abs(drap?!);');
  c_eval('drap?!(:,1)=R?(iap?,1);');
  save -append mRcluster drap12 drap13 drap14 drap23 drap24 drap34;

  %% Find tail season an dayside season
  c_eval('itail?=find(Rap?(2:end,2)<0 & (Rap?(1:end-1,3).*Rap?(2:end,3) <0) );')
  c_eval('ttail?=Rap?(itail?,1);');
  c_eval('iday?=find(Rap?(2:end,2)>0 & (Rap?(1:end-1,3).*Rap?(2:end,3) <0 ));')
  c_eval('tday?=Rap?(iday?,1);');
  c_eval('save -append mRcluster ttail? tday?');
end

%% Get Cluster file from internet
if flag_get_Cluster_file_from_www
  disp('Downloading from internet Cluster orbit parameter file mRcluster.mat ...');
  if verLessThan('matlab', '8.4')
    [f,status]=urlwrite('https://space.irfu.se/~andris/cluster/mRcluster.mat','mRcluster.mat'); %#ok<URLWR> websave introduced in R2014b
  else
    try
      f = websave('mRcluster.mat', 'https://space.irfu.se/~andris/cluster/mRcluster.mat');
      status = true;
    catch
      status = false;
    end
  end
end

%% Create the figure
if flag_create_figure % create the figure
  load mRcluster
  dtmark_season=365/24*2*24*3600; % seconds corresponding to 2h drift of perigee
  cluster_color=[[0 0 0];[1 0 0];[0 0.5 0];[0 0 1]];
  hca=irf_plot(1,'newfigure');
  set(gcf,'position',[10 10 700 450]);
  set(hca,'yscale','log');
  hold on;
  for ic1=1:4
    for ic2=ic1+1:4
      c_eval('irf_plot(hca,drap?!(:,[1 5]),''.'',''color'',cluster_color(!,:),''markersize'',10)',ic1,ic2);
      c_eval('irf_plot(hca,drap?!(:,[1 5]),''color'',cluster_color(?,:),''linewidth'',1)',ic1,ic2);
    end
  end
  irf_zoom(hca,'x',[drap12(1,1) drap12(end,1)]);
  irf_pl_mark(hca,[ttail1-dtmark_season ttail1+dtmark_season],[0.5 0.5 0.5]);
  irf_pl_mark(hca,[tday1-dtmark_season tday1+dtmark_season],'yellow');
  irf_legend(hca,{'C1','C2','C3','C4'},[0.05 0.95],'color','cluster')
  ylabel(hca,'Separation [km]');
  title(hca,'Cluster separation at apogee. Seasons - tail (gray), dayside (yellow).');
  irf_legend(0,'irfu-matlab Example_Cluster_separation_during_mission ',[0.02 0.02],'interpreter','none','color',[0.7 0.7 0.7],'fontsize',8)
end

%% Print figure
if flag_print_figure
  set(gcf,'paperpositionmode','auto')
  print -dpng Cluster_separation.png
end


