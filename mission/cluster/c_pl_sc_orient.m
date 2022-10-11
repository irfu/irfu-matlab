function hout=c_pl_sc_orient(spacecraft,time,phase_time_series,magnetic_field,velocity,action)
%C_PL_SC_ORIENT   Plots the orientation of the EFW probes
%
%   h = C_PL_SC_ORIENT;
%   h = C_PL_SC_ORIENT(ic);
%   h = C_PL_SC_ORIENT(ic,t);
%   h = C_PL_SC_ORIENT(ic,t,a);
%   h = C_PL_SC_ORIENT(ic,t,a,b);
%   h = C_PL_SC_ORIENT(ic,t,a,b,v);
%   ic - spacecraft number
%   t  - time in isdat epoch
%   a  - time vector of the satellite phase in degrees
%   b  - magnetic field in GSE reference frame
%   v  - velocity vector [vx vy vz] in GSE which will be marked in the
%        plots, e.g. magnetopause velocity

%% Default values
defaultTime = [2010 12 31 01 01 01];
defaultSC   = 1;
%% Check inputs
if nargin==1 && isnumeric(spacecraft)
  time = defaultTime;
  action='initialize';
elseif nargin==1 && ischar(spacecraft)
  action=spacecraft;
elseif   (nargin < 6)
  action='initialize';
end
if nargin==0
  spacecraft=defaultSC;
  isTimeOK = false;
  if evalin('caller','exist(''tint'')')
    ttt = evalin('caller','tint(1)');
    if ttt > irf_time([2001 1 1 0 0 0])
      isTimeOK = true;
      time=irf_time(ttt,'vector');
    end
  end
  if ~isTimeOK && exist('CAA','dir')
    R=irf_get_data('sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M','caa','mat');
    if numel(R)
      time=0.5*(R(1,1)+R(end,1)); % first point in center of position time series
    end
  end
  if ~isTimeOK % define default time
    time=defaultTime;
  end
end
if isempty(spacecraft)
  ic=defaultSC;
else
  ic=spacecraft;
end

%% Choose action
irf.log('debug',['action=' action]);
switch lower(action)
  case 'initialize'
    % See if spacecraft orientation figures is open
    figNumber=figure( ...
      'Name',['Cluster s/c' num2str(ic) ' orientation'], ...
      'Tag','cplscor');
    set(figNumber,'Position',[10 10 600 600]);
    delete(findall(gcf,'Type','axes'))
    h=[];
    data.figNumber=figNumber;
    data.ic=ic;
    data.getScPhase=1;
    data.getB=1;
    if length(time)==1 % define time of interest when initializing
      data.t=time;
    elseif length(time)==6
      data.t=irf_time(time);
    elseif exist('t','var')
      % use existing t
    else
      irf.log('critical','unknown time format');
      return;
    end
    if nargin<6, data.flag_v=1; end
    if nargin<5, data.flag_v=0; end
    if nargin>4 % use mangetic field that is given as input
      data.b=c_coord_trans('GSE','ISR2',magnetic_field,'cl_id',data.ic);
    end
    if nargin>2 % use phase given as input
      data.a=phase_time_series;
    end
    if data.flag_v == 1, data.v=velocity; end
    set(data.figNumber,'defaultAxesFontSize',12);
    set(data.figNumber,'defaultTextFontSize',10);
    set(data.figNumber,'defaultLineLineWidth',1);
    h(1)=subplot(2,2,1);axis equal;axis([-50 50 -50 50]);axis manual;title('Spin plane');
    h(2)=subplot(2,2,2);axis equal;axis([-50 50 -50 50]);axis manual;title('View along B');
    h(3)=subplot(2,2,3);axis equal;axis([-50 50 -50 50]);axis manual;title('View towards sun');
    h(4)=subplot(2,2,4);axis off;
    data.h=h;
    set(figNumber,'userdata',data);
    
    %====================================
    % The vector 1 entering
    labelStr='0';
    callbackStr='c_pl_sc_orient(''plot'')';
    data.vec1flag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.2 .2 .05],'string','vector 1 [GSE]','Callback',callbackStr);
    data.vec1Hndl=uicontrol( ...
      'Style','edit', ...
      'Units','normalized', ...
      'Position',[0.7 0.2 .2 .05], ...
      'String',labelStr, ...
      'Callback',callbackStr);
    %====================================
    % The B field reading
    callbackStr='c_pl_sc_orient(''read_phase_and_b'')';
    data.bflag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.3 .2 .05],'string','read B field','Callback',callbackStr);
    %====================================
    % The vector 2 entering
    labelStr='0';
    callbackStr='c_pl_sc_orient(''plot'')';
    data.vec2flag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.25 .2 .05],'string','vector 2 [GSE]','Callback',callbackStr);
    data.vec2Hndl=uicontrol( ...
      'Style','edit', ...
      'Units','normalized', ...
      'Position',[0.7 0.25 .2 .05], ...
      'String',labelStr, ...
      'Callback',callbackStr);
    %====================================
    % The phase entering
    labelStr='0';
    callbackStr='c_pl_sc_orient(''phase'')';
    uicontrol('style','text','units','normalized','Position',[0.5 0.15 .2 .05],'string','phase [deg]')
    data.phaseHndl=uicontrol( ...
      'Style','edit', ...
      'Units','normalized', ...
      'Position',[0.7 0.15 .2 .05], ...
      'String',labelStr, ...
      'Callback',callbackStr);
    %====================================
    % The time entering
    labelStr=mat2str(irf_time(data.t,'vector'));
    callbackStr='c_pl_sc_orient(''time'')';
    data.timeHndl=uicontrol( ...
      'Style','edit', ...
      'Units','normalized', ...
      'Position',[0.5 0.1 .3 .05], ...
      'String',labelStr, ...
      'Callback',callbackStr);
    %====================================
    % The CLOSE button
    labelStr='Close';
    callbackStr='close(gcf)';
    data.closeHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',[0.5 0 .1 .05], ...
      'String',labelStr, ...
      'Callback',callbackStr);
    menus
    set(gcf,'userdata',data);
    c_pl_sc_orient('read_phase_and_b');
  case 'read_phase_and_b'
    data=get(gcf,'userdata');
    tint = data.t + [-120 120]; % read in 4min data
    if get(data.bflag,'value')==1
      data.getB = true;
      if ~isfield(data,'b')
        data.b=[1 0 0 NaN]; % first col is time
      end
      if ~isempty(data.b) % use existing b if there is one
        if data.t>=data.b(1,1) && data.t<=data.b(end,1) % time within interval of B
          data.getB = false;
        elseif strcmp(get(data.bflag,'value'),1) % get B data
          data.getB = true;
        else
          data.b=[1 0 0 NaN]; % first col is time
        end
      end
      if data.getB % try to get B data from disk mat files
        irf.log('debug','Trying to to read B from disk as diB?');
        [ok,b]=c_load('diB?',data.ic);
        if any(ok) && data.t>=b(1,1) && data.t<=b(end,1) % b loaded and time within interval of B
          data.b = b;
          data.b(:,3)=-b(:,3);data.b(:,4)=-b(:,4); % go to DS reference frame instead of DSI
          data.getB = false;
        end
      end
      if data.getB % try to get B data from caa files full resolution
        irf.log('debug','Trying to get full resolution B');
        c_eval('b=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'',''tint'',tint);',data.ic);
        if ~isempty(b)
          data.b=c_coord_trans('GSE','DSC',b,'cl_id',data.ic);
          if data.t>=data.b(1,1) && data.t<=data.b(end,1) % time within interval of B
            data.getB = false;
          end
        end
      end
      if data.getB % try to get B data from caa files 5S/s resolution
        c_eval('b=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'',''mat'',''tint'',tint);',data.ic);
        if ~isempty(b)
          data.b=c_coord_trans('GSE','DSC',b,'cl_id',data.ic);
          if data.t>=data.b(1,1) && data.t<=data.b(end,1) % time within interval of B
            data.getB = false;
          end
        end
      end
      if data.getB % try to read B in ISR2 ref frame from isdat (use CSDD PP data)
        DATABASE=c_ctl(0,'isdat_db');
        isdatdata = getData(ClusterDB(DATABASE,c_ctl(0,'data_path')),data.t-5,5,data.ic,'b','nosave');
        if ~isempty(isdatdata)
          b=isdatdata{3};
          data.b=c_coord_trans('GSE','DSC',b,'cl_id',data.ic);
          if data.t>=data.b(1,1) && data.t<=data.b(end,1) % time within interval of B
            data.getB = false;
          end
        end
      end
      if data.getB % did not succeed to read B data
        irf.log('warning','Could not read B field data');
        data.b=[1 0 0 NaN]; % first col is time
      end
    end
    if ~data.getScPhase % can use old phase data if they are valid
      if data.t<=data.a(1,1) || data.t>=data.a(end,1) % time outside interval of phase
        data.getScPhase = true;
      end
    end
    if data.getScPhase % try to read phase data from disk matlab files
      irf.log('debug','Trying to to read phase from disk as Atwo?');
      [ok,phase_time_series]=c_load('Atwo?',data.ic);
      if any(ok) % check if phase info is on disk
        data.a=phase_time_series;
        if data.t>=data.a(1,1) && data.t<=data.a(end,1) % time within phase interval of B
          data.getScPhase = false;
        end
      end
    end
    if data.getScPhase % try to read caa phase if exists
      irf.log('debug','Trying to to read phase from CP_AUX_SPIN_TIME');
      spinPeriod=c_caa_var_get(['spin_period__C' num2str(data.ic) '_CP_AUX_SPIN_TIME'],'mat','tint',tint);
      if ~isempty(spinPeriod) % could not read
        spinTime=c_caa_var_get(['time_tags__C' num2str(data.ic) '_CP_AUX_SPIN_TIME'],'mat','tint',tint);
        refTime=[-.5 .5]; % part of spin
        refPhase=refTime*360+180; % spin period center corresponds to phase 0
        dt= repmat(refTime,size(spinTime,1),1).*...
          repmat(double(spinPeriod(:,2)),1,length(refTime));
        tmat=repmat(spinTime(:,1),1,length(refTime))+dt;
        amat=repmat(refPhase,size(spinTime,1),1);
        tmat=reshape(tmat',numel(tmat),1);
        difftmat=diff(tmat);
        ii=find(difftmat<0);
        tmat(ii)=tmat(ii+1);
        amat=reshape(amat',numel(amat),1);
        data.a=[tmat amat];
        if data.t>=data.a(1,1) && data.t<=data.a(end,1) % time within phase interval of B
          data.getScPhase = false;
        end
      end
    end
    if data.getScPhase % try to read phase data from isdat server
      try
        c_eval('[tt,phase_data] = irf_isdat_get([''Cluster/?/ephemeris/phase_2''], data.t-5, 10);',data.ic);
        %                [tt,phase_data] = caa_is_get('db.irfu.se:0',t-5,10,ic,'ephemeris','phase_2');
        phase_time_series=[tt phase_data];
      catch
        phase_time_series=[]; % cannot connect to internet
      end
      if ~isempty(phase_time_series)
        data.getScPhase = false;
        data.a=phase_time_series;
      end
    end
    if data.getScPhase % still no phase info, use default 0
      data.phase=0; % default using 0 phase
    else
      data.phase=c_phase(data.t,data.a);data.phase(1)=[];
    end
    set(data.phaseHndl,'string',num2str(data.phase,'%3.1f'));
    set(gcf,'userdata',data);
    c_pl_sc_orient('plot');
  case 'time'
    data=get(gcf,'userdata');
    data.t=toepoch(eval(get(data.timeHndl, 'string')));
    set(gcf,'userdata',data);
    c_pl_sc_orient('read_phase_and_b');
  case 'phase'
    data=get(gcf,'userdata');
    data.phase=str2double(get(data.phaseHndl, 'string'));
    set(gcf,'userdata',data);
    c_pl_sc_orient('plot');
  case 'plot'
    data=get(gcf,'userdata');
    h=data.h;
    data.flag_v1=get(data.vec1flag, 'value');
    data.flag_v2=get(data.vec2flag, 'value');
    data.flag_b=get(data.bflag, 'value');
    
    if data.flag_v1==1
      data.v1=eval(['[' get(data.vec1Hndl,'string') ']']);
      if length(data.v1)==1, data.flag_v1=0;end
    end
    if data.flag_v2==1
      data.v2=eval(['[' get(data.vec2Hndl,'string') ']']);
      if length(data.v2)==1, data.flag_v2=0;end
    end
    
    phase_p1=data.phase/180*pi + 3*pi/4 ;
    phase_p3=phase_p1     - pi/2   ;
    phase_p2=phase_p1     + pi     ;
    phase_p4=phase_p1     + pi/2 ;
    phase_heea=data.phase/180*pi-(30)/180*pi;
    phase_leea=phase_heea+pi;
    phase_rapid=data.phase/180*pi + 60.167/180*pi; % rapid phase
    phase_sunsensor=data.phase/180*pi + 26.367/180*pi; % the location o fsun sensor
    
    rp1=[44*cos(phase_p1) 44*sin(phase_p1) 0]; %#ok<NASGU> % in DS reference frame
    rp2=[44*cos(phase_p2) 44*sin(phase_p2) 0]; %#ok<NASGU>
    rp3=[44*cos(phase_p3) 44*sin(phase_p3) 0]; %#ok<NASGU>
    rp4=[44*cos(phase_p4) 44*sin(phase_p4) 0]; %#ok<NASGU>
    dphi=5/180*pi; % the half size of heea leea rapid azimuthal sectors that are plotted
    sec_length=15; % the length of plotted sectors (the length of efw booms is 44)
    rheea=sec_length*[cos(phase_heea) sin(phase_heea);cos(phase_heea-dphi) sin(phase_heea-dphi);cos(phase_heea+dphi) sin(phase_heea+dphi)];
    rleea=sec_length*[cos(phase_leea) sin(phase_leea);cos(phase_leea-dphi) sin(phase_leea-dphi);cos(phase_leea+dphi) sin(phase_leea+dphi)];
    rrapid=sec_length*[cos(phase_rapid) sin(phase_rapid);cos(phase_rapid-dphi) sin(phase_rapid-dphi);cos(phase_rapid+dphi) sin(phase_rapid+dphi)];
    rsunsensor=sec_length*[cos(phase_sunsensor) sin(phase_sunsensor)];
    
    for ip=1:4, c_eval('rp?_gse=c_coord_trans(''DSC'',''GSE'',[data.t rp?],''cl_id'',data.ic);rp?_gse(1)=[];',ip),end
    
    if data.flag_b==1
      bfield=irf_resamp(data.b,data.t);
      bxs=irf_norm(irf_cross(bfield,[0 0 0 1]));
      bxsxb=irf_norm(irf_cross(bxs,bfield)); %#ok<NASGU> % (bxs)xb
      bn=irf_norm(bfield);
      bn_gse=c_coord_trans('DSC','GSE',bn,'cl_id',data.ic);
      b_elevation=-asin(bn(4))*180/pi;
      angle_deg_p34_vs_b=acos(bn(2)*cos(phase_p4)+bn(3)*sin(phase_p4))*180/pi; % acos(bx*rx+by*ry)
      angle_deg_p12_vs_b=acos(bn(2)*cos(phase_p2)+bn(3)*sin(phase_p2))*180/pi;
      for ip=1:4,c_eval('rp?_b=[irf_dot(rp?,bxs,1) irf_dot(rp?,bxsxb,1) irf_dot(rp?,bn,1)];',ip),end
    end
    if data.flag_v1==1
      vn1_gse=[bn(1,1) irf_norm(data.v1)];
      vn1_ds=c_coord_trans('GSE','DSC',vn1_gse,'cl_id',data.ic);
      vn1_elevation=-asin(vn1_ds(4))*180/pi;
    end
    if data.flag_v2==1
      vn2_gse=[bn(1,1) irf_norm(data.v2)];
      vn2_ds=c_coord_trans('GSE','DSC',vn2_gse,'cl_id',data.ic);
      vn2_elevation=-asin(vn2_ds(4))*180/pi;
    end
    
    aa=0:.1:2*pi;x_circle=cos(aa);y_circle=sin(aa);
    
    axes(h(1));cla
    text(0,50,'sun','verticalalignment','top','horizontalalignment','center','fontweight','demi');
    text(50,0,'dawn','rotation',90,'verticalalignment','bottom','horizontalalignment','center','fontweight','demi');
    patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
    patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
    patch([0 rheea(2,2) rheea(3,2)], [0 rheea(2,1) rheea(3,1)],'b'); % heea
    patch([0 rleea(2,2) rleea(3,2)], [0 rleea(2,1) rleea(3,1)],'g'); % leea
    patch([0 rrapid(2,2) rrapid(3,2)], [0 rrapid(2,1) rrapid(3,1)],'k'); % rapid
    line([0 rsunsensor(2)],[0 rsunsensor(1)],[0 0],'Color',[1 0.5 .5],'LineWidth',2); % sun sensor direction
    text(50,50,'Sun sensor','verticalalignment','top','horizontalalignment','right','fontweight','demi','color',[1 0.5 .5]);
    text(50,45,'LEEA','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','g');
    text(50,40,'HEEA','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','b');
    text(50,35,'Rapid','verticalalignment','top','horizontalalignment','right','fontweight','demi','color','k');
    
    if data.flag_b==1
      bnproj=[0 bn(2)/norm(bn(2:3)) bn(3)/norm(bn(2:3))];
      hl=line([0 bnproj(3)*25],[0 bnproj(2)*25]);set(hl,'color','red','linewidth',.4);       % B direction
      hl=line([0 bn(3)*25],[0 bn(2)*25]);set(hl,'color','red','linewidth',2);       % B direction
      text(30*bnproj(3),30*bnproj(2),'B');           % label
      text(-49,-48,['B elevation=' num2str(b_elevation,2) ' deg'],'fontsize',8);           % label
    end
    if data.flag_v1==1 % plot v1 vector
      vn1proj=[0 vn1_ds(2)/norm(vn1_ds(2:3)) vn1_ds(3)/norm(vn1_ds(2:3))];
      hl=line([0 vn1proj(3)*25],[0 vn1proj(2)*25]);set(hl,'color','k','linewidth',.4);       % V direction
      hl=line([0 vn1_ds(3)*25],[0 vn1_ds(2)*25]);set(hl,'color','k','linewidth',2);       % V direction
      text(30*vn1proj(3),30*vn1proj(2),'V');           % label
      text(-49,44,['V1 elevation=' num2str(vn1_elevation,2) ' deg'],'fontsize',8);           % label
    end
    if data.flag_v2==1 % plot v1 vector
      vn2proj=[0 vn2_ds(2)/norm(vn2_ds(2:3)) vn2_ds(3)/norm(vn2_ds(2:3))];
      hl=line([0 vn2proj(3)*25],[0 vn2proj(2)*25]);set(hl,'color','b','linewidth',.4);       % V direction
      hl=line([0 vn2_ds(3)*25],[0 vn2_ds(2)*25]);set(hl,'color','b','linewidth',2);       % V direction
      text(30*vn2proj(3),30*vn2proj(2),'V');           % label
      text(-49,38,['V2 elevation=' num2str(vn2_elevation,2) ' deg'],'fontsize',8);           % label
    end
    
    for aa=0:pi/12:2*pi % plot grid
      hl=line([0 100*cos(aa)],[0 100*sin(aa)]);
      set(hl,'linestyle',':','color','green','linewidth',.2);
    end
    for ip=1:4 % plot probes
      c_eval('line([0 rp?(2)],[0 rp?(1)]);',ip);
      c_eval('patch(rp?(2)+x_circle*0.4,rp?(1)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip);
      c_eval('text(rp?(2)*.9,rp?(1)*.9,num2str(?),''fontweight'',''bold'');',ip);
    end
    
    axes(h(3));cla
    text(0,50,'Z_{GSE}','verticalalignment','top','horizontalalignment','center','fontweight','demi');
    text(50,0,'-Y_{GSE}','rotation',90,'verticalalignment','bottom','horizontalalignment','center','fontweight','demi');
    patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
    patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
    if data.flag_b==1
      bnproj=[0 0 bn_gse(3)/norm(bn_gse(3:4)) bn_gse(4)/norm(bn_gse(3:4))];
      hl=line([0 -bnproj(3)*25],[0 bnproj(4)*25]);set(hl,'color','red','linewidth',.4);       % B direction
      hl=line([0 -bn_gse(3)*25],[0 bn_gse(4)*25]);set(hl,'color','red','linewidth',2);       % B direction
      text(-30*bnproj(3),30*bnproj(4),'B');           % label
    end
    if data.flag_v1==1 % plot v vector
      vn1proj=[0 0 vn1_gse(3)/norm(vn1_gse(3:4)) vn1_gse(4)/norm(vn1_gse(3:4))];
      hl=line([0 -vn1proj(3)*25],[0 vn1proj(4)*25]);set(hl,'color','k','linewidth',.4);       % V direction
      hl=line([0 -vn1_gse(3)*25],[0 vn1_gse(4)*25]);set(hl,'color','k','linewidth',2);       % V direction
      text(-30*vn1proj(3),30*vn1proj(4),'V');           % label
    end
    if data.flag_v2==1 % plot v vector
      vn2proj=[0 0 vn2_gse(3)/norm(vn2_gse(3:4)) vn2_gse(4)/norm(vn2_gse(3:4))];
      hl=line([0 -vn2proj(3)*25],[0 vn2proj(4)*25]);set(hl,'color','b','linewidth',.4);       % V direction
      hl=line([0 -vn2_gse(3)*25],[0 vn2_gse(4)*25]);set(hl,'color','b','linewidth',2);       % V direction
      text(-30*vn2proj(3),30*vn2proj(4),'V');           % label
    end
    
    for aa=0:pi/12:pi/2
      hl=line(x_circle*44*sin(aa),y_circle*44*sin(aa));
      set(hl,'linestyle',':','color','green','linewidth',.2);
    end
    for ip=1:4
      c_eval('line([0 -rp?_gse(2)],[0 rp?_gse(3)]);',ip);
      c_eval('patch(-rp?_gse(2)+x_circle*0.4,rp?_gse(3)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip);
      c_eval('text(-rp?_gse(2)*.8,rp?_gse(3)*.8,num2str(?));',ip);
    end
    
    axes(h(2));cla
    text(0,50,'(BxS)xS','verticalalignment','top','horizontalalignment','center','fontweight','demi');
    text(50,0,'BxS','rotation',90,'verticalalignment','bottom','horizontalalignment','center','fontweight','demi');
    patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
    patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
    for aa=0:pi/12:pi/2
      hl=line(x_circle*44*sin(aa),y_circle*44*sin(aa));
      set(hl,'linestyle',':','color','green','linewidth',.2);
    end
    if data.flag_b==1
      for ip=1:4
        c_eval('line([0 rp?_b(1)],[0 rp?_b(2)]);',ip);
        c_eval('patch(rp?_b(1)+x_circle*0.4,rp?_b(2)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip);
        c_eval('text(rp?_b(1)*.8,rp?_b(2)*.8,num2str(?));',ip);
      end
    end
    
    % add text
    cla(h(4))
    irf_legend(h(4),['c_pl_sc_orient() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],[0,1],'fontsize',8,'interpreter','none','color',[0.5 0.5 0.5]);
    irf_legend(h(4),['Cluster spacecraft C' num2str(data.ic)],[0,.9],'fontsize',10);
    if data.flag_b==1
      irf_legend(h(4),['angle beteen p34 (Ey WBD) and B ' num2str(angle_deg_p34_vs_b,'%3.1f') 'deg'],[0,.82],'interpreter','none','fontsize',10);
      irf_legend(h(4),['angle beteen p12 (Ez WBD) and B ' num2str(angle_deg_p12_vs_b,'%3.1f') 'deg'],[0,.74],'interpreter','none','fontsize',10);
    end
    set(gcf,'userdata',data);
  case 'c1'
    data=get(gcf,'userdata');
    if ic~=1
      data.getScPhase=1;data.getB=1;data.b=[];data.ic=1;
      set(gcf,'Name',['Cluster s/c' num2str(data.ic) ' orientation']);
      set(gcf,'userdata',data);
      c_pl_sc_orient('read_phase_and_b');
    end
  case 'c2'
    data=get(gcf,'userdata');
    if ic~=2
      data.getScPhase=1;data.getB=1;data.b=[];data.ic=2;
      set(gcf,'Name',['Cluster s/c' num2str(data.ic) ' orientation']);
      set(gcf,'userdata',data);
      c_pl_sc_orient('read_phase_and_b');
    end
  case 'c3'
    data=get(gcf,'userdata');
    if ic~=3
      data.getScPhase=1;data.getB=1;data.b=[];data.ic=3;
      set(gcf,'Name',['Cluster s/c' num2str(data.ic) ' orientation']);
      set(gcf,'userdata',data);
      c_pl_sc_orient('read_phase_and_b');
    end
  case 'c4'
    data=get(gcf,'userdata');
    if ic~=4
      data.getScPhase=1;data.getB=1;data.b=[];data.ic=4;
      set(gcf,'Name',['Cluster s/c' num2str(ic) ' orientation']);
      set(gcf,'userdata',data);
      c_pl_sc_orient('read_phase_and_b');
    end
end
if nargout,hout=h;end

function menus
% generate menus
if isempty(findobj(gcf,'type','uimenu','label','&Options'))
  hcoordfigmenu=uimenu('Label','&Options');
  uimenu(hcoordfigmenu,'Label','C1','Callback','c_pl_sc_orient(''C1'')')
  uimenu(hcoordfigmenu,'Label','C2','Callback','c_pl_sc_orient(''C2'')')
  uimenu(hcoordfigmenu,'Label','C3','Callback','c_pl_sc_orient(''C3'')')
  uimenu(hcoordfigmenu,'Label','C4','Callback','c_pl_sc_orient(''C4'')')
  user_data = get(gcf,'userdata');
  user_data.coordfigmenu=1;
  set(gcf,'userdata',user_data);
end

