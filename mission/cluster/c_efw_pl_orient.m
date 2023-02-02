function h=c_efw_pl_orient(spacecraft,time,phase_time_series,magnetic_field,velocity,action)
%C_EFW_PL_ORIENT   Plot the orientation of EFW probes
%
%   h = C_EFW_PL_ORIENT;
%   h = C_EFW_PL_ORIENT(ic);
%   h = C_EFW_PL_ORIENT(ic,t);
%   h = C_EFW_PL_ORIENT(ic,t,a);
%   h = C_EFW_PL_ORIENT(ic,t,a,b);
%   h = C_EFW_PL_ORIENT(ic,t,a,b,v);
%   ic - spacecraft number
%   t  - time in isdat epoch
%   a  - time vector of the satellite phase in degrees
%   b  - magnetic field in despinned reference frame
%   v  - velocity vector [vx vy vz] in GSE which will be marked in the plots,
%        e.g. magnetopause velocity
%

%   figuserdata=[h];
eval_figuserdata='figuserdata={h};';

persistent t a b phase v ic phaseHndl timeHndl figNumber ...
  vec1Hndl vec2Hndl vec1flag vec2flag ...
  flag_v1 flag_v2 v1 v2;
if       (nargin==1 && ischar(spacecraft)), action=spacecraft;irf_log('proc',['action=' action]);
elseif   (nargin < 6)                   , action='initialize';
end

if strcmp(action,'initialize')
  if nargin<1, help c_efw_pl_orient;return;                                            end
  ic=spacecraft;
  if nargin<6, flag_v=1;                                                                   end
  if nargin<5, flag_v=0;                                                                   end
  if nargin<4
    if  exist('mB.mat','file')
      ok=c_load('B?',ic);
      if ok, c_eval('magnetic_field=c_gse2dsc(B?,?);clear B?',ic);end
    elseif exist('mBPP.mat','file'), eval(irf_ssub('load mBPP dBPP?;magnetic_field=dBPP?;clear dBPP?',ic));
    else,  irf_log('load','Could not read B field, using B=[0 0 1] nT in DS ref frame');magnetic_field=[1 0 0 1]; % first col is time
    end
  end
  if nargin<3, eval(irf_ssub('load mA.mat Atwo?;phase_time_series=Atwo?;clear Atwo?',ic));          end
  if nargin<2, time=phase_time_series(1,1);                                                end
  t=time;a=phase_time_series;b=magnetic_field;
  if flag_v == 1, v=velocity; end
  % See if spacecraft orientation figures is open
  ch = get(0,'ch');indx=[];
  if ~isempty(ch)
    chTags = get(ch,'Tag');
    indx = find(strcmp(chTags,'cplscor'));
  end
  if isempty(indx)
    figNumber=figure( ...
      'Name',['Cluster s/c' num2str(ic) ' orientation'], ...
      'Tag','cplscor');
  else
    figure(ch(indx));clf;figNumber=gcf;
  end
  
  h(1)=subplot(2,2,1);axis equal;axis([-50 50 -50 50]);axis manual;title(['s/c' num2str(spacecraft) ' spin plane']);
  h(2)=subplot(2,2,2);axis equal;axis([-50 50 -50 50]);axis manual;title('View along B');
  h(3)=subplot(2,2,3);axis equal;axis([-50 50 -50 50]);axis manual;title('View towards sun');
  h(4)=subplot(2,2,4);axis off;
  %====================================
  % The vector 1 entering
  labelStr='0';
  callbackStr='c_efw_pl_orient(''plot'')';
  vec1flag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.2 .2 .05],'string','show vector 1 [GSE]','Callback',callbackStr);
  vec1Hndl=uicontrol( ...
    'Style','edit', ...
    'Units','normalized', ...
    'Position',[0.7 0.2 .2 .05], ...
    'String',labelStr, ...
    'Callback',callbackStr);
  %====================================
  % The vector 2 entering
  labelStr='0';
  callbackStr='c_efw_pl_orient(''plot'')';
  vec2flag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.25 .2 .05],'string','show vector 2 [GSE]','Callback',callbackStr);
  vec2Hndl=uicontrol( ...
    'Style','edit', ...
    'Units','normalized', ...
    'Position',[0.7 0.25 .2 .05], ...
    'String',labelStr, ...
    'Callback',callbackStr);
  %====================================
  % The phase entering
  labelStr='0';
  callbackStr='c_efw_pl_orient(''phase'')';
  uicontrol('style','text','units','normalized','Position',[0.5 0.15 .2 .05],'string','phase [deg]')
  phaseHndl=uicontrol( ...
    'Style','edit', ...
    'Units','normalized', ...
    'Position',[0.7 0.15 .2 .05], ...
    'String',labelStr, ...
    'Callback',callbackStr);
  %====================================
  % The time entering
  labelStr=[epoch2iso(t)];
  callbackStr='c_efw_pl_orient(''time'')';
  timeHndl=uicontrol( ...
    'Style','edit', ...
    'Units','normalized', ...
    'Position',[0.5 0.1 .3 .05], ...
    'String',labelStr, ...
    'Callback',callbackStr);
  %====================================
  % The CLOSE button
  labelStr='Close';
  callbackStr='close(gcf)';
  closeHndl=uicontrol( ...
    'Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[0.5 0 .1 .05], ...
    'String',labelStr, ...
    'Callback',callbackStr);
  eval(eval_figuserdata);
  set(figNumber,'UserData',figuserdata);
  c_efw_pl_orient('time');
  
elseif strcmp(action,'time')
  t=iso2epoch(get(timeHndl, 'string'));
  phase = c_phase(t,a);
  phase = phase(2);
  set(phaseHndl,'string',num2str(phase,'%.1f'));
  c_efw_pl_orient('plot');
  
elseif strcmp(action,'phase')
  phase=str2num(get(phaseHndl, 'string'));
  c_efw_pl_orient('plot');
  
  
elseif strcmp(action,'plot')
  figuserdata=get(figNumber,'userdata');
  h=figuserdata{1};
  
  flag_v1=get(vec1flag, 'value');
  if flag_v1==1, v1=eval(get(vec1Hndl,'string'));if length(v1)==1, flag_v1=0;end;end
  flag_v2=get(vec2flag, 'value');
  if flag_v2==1, v2=eval(get(vec2Hndl,'string'));if length(v2)==1, flag_v2=0;end;end
  
  phase_p1=phase/180*pi + 3*pi/4 ;
  phase_p3=phase_p1     - pi/2   ;
  phase_p2=phase_p1     + pi     ;
  phase_p4=phase_p1     + pi/2 ;
  rp1=[44*cos(phase_p1) 44*sin(phase_p1) 0]; %#ok<NASGU> % in DS reference frame
  rp2=[44*cos(phase_p2) 44*sin(phase_p2) 0]; %#ok<NASGU>
  rp3=[44*cos(phase_p3) 44*sin(phase_p3) 0]; %#ok<NASGU>
  rp4=[44*cos(phase_p4) 44*sin(phase_p4) 0]; %#ok<NASGU>
  
  for ip=1:4,eval(irf_ssub('rp?_gse=c_gse2dsc([t rp?],ic,-1);rp?_gse(1)=[];',ip)),end
  bfield=irf_resamp(b,t);
  bxs=irf_norm(irf_cross(bfield,[0 0 0 1]));
  bxsxb=irf_norm(irf_cross(bxs,bfield)); %#ok<NASGU>  % (bxs)xb
  bn=irf_norm(bfield);
  bn_gse=c_gse2dsc(bn,ic,-1);
  b_elevation=-asin(bn(4))*180/pi;
  
  if flag_v1==1
    vn1_gse=[bn(1,1) irf_norm(v1)];
    vn1_ds=c_gse2dsc(vn1_gse,ic);
    vn1_elevation=-asin(vn1_ds(4))*180/pi;
  end
  if flag_v2==1
    vn2_gse=[bn(1,1) irf_norm(v2)];
    vn2_ds=c_gse2dsc(vn2_gse,ic);
    vn2_elevation=-asin(vn2_ds(4))*180/pi;
  end
  
  for ip=1:4,eval(irf_ssub('rp?_b=[irf_dot(rp?,bxs,1) irf_dot(rp?,bxsxb,1) irf_dot(rp?,bn,1)];',ip)),end
  
  aa=0:.1:2*pi;x_circle=cos(aa);y_circle=sin(aa);
  
  axes(h(1));cla
  text(0,50,'sun','verticalalignment','middle','horizontalalignment','center','fontweight','demi');
  text(50,0,'dawn','rotation',90,'verticalalignment','middle','horizontalalignment','center','fontweight','demi');
  patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
  patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
  bnproj=[0 bn(2)/norm(bn(2:3)) bn(3)/norm(bn(2:3))];
  hl=line([0 bnproj(3)*25],[0 bnproj(2)*25]);set(hl,'color','red','linewidth',.4);       % B direction
  hl=line([0 bn(3)*25],[0 bn(2)*25]);set(hl,'color','red','linewidth',2);       % B direction
  text(30*bnproj(3),30*bnproj(2),'B');           % label
  text(-49,48,['B elevation=' num2str(b_elevation,2) ' deg'],'fontsize',8);           % label
  if flag_v1==1 % plot v1 vector
    vn1proj=[0 vn1_ds(2)/norm(vn1_ds(2:3)) vn1_ds(3)/norm(vn1_ds(2:3))];
    hl=line([0 vn1proj(3)*25],[0 vn1proj(2)*25]);set(hl,'color','k','linewidth',.4);       % V direction
    hl=line([0 vn1_ds(3)*25],[0 vn1_ds(2)*25]);set(hl,'color','k','linewidth',2);       % V direction
    text(30*vn1proj(3),30*vn1proj(2),'V');           % label
    text(-49,44,['V1 elevation=' num2str(vn1_elevation,2) ' deg'],'fontsize',8);           % label
  end
  if flag_v2==1 % plot v1 vector
    vn2proj=[0 vn2_ds(2)/norm(vn2_ds(2:3)) vn2_ds(3)/norm(vn2_ds(2:3))];
    hl=line([0 vn2proj(3)*25],[0 vn2proj(2)*25]);set(hl,'color','b','linewidth',.4);       % V direction
    hl=line([0 vn2_ds(3)*25],[0 vn2_ds(2)*25]);set(hl,'color','b','linewidth',2);       % V direction
    text(30*vn2proj(3),30*vn2proj(2),'V');           % label
    text(-49,38,['V2 elevation=' num2str(vn2_elevation,2) ' deg'],'fontsize',8);           % label
  end
  
  for aa=0:pi/12:2*pi
    hl=line([0 100*cos(aa)],[0 100*sin(aa)]);
    set(hl,'linestyle',':','color','green','linewidth',.2);
  end
  for ip=1:4
    eval(irf_ssub('line([0 rp?(2)],[0 rp?(1)]);',ip));
    eval(irf_ssub('patch(rp?(2)+x_circle*0.4,rp?(1)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip));
    eval(irf_ssub('text(rp?(2)*.8,rp?(1)*.8,num2str(?));',ip));
  end
  
  axes(h(3));cla
  text(0,50,'Z_{GSE}','verticalalignment','middle','horizontalalignment','center','fontweight','demi');
  text(50,0,'-Y_{GSE}','rotation',90,'verticalalignment','middle','horizontalalignment','center','fontweight','demi');
  patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
  patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
  bnproj=[0 0 bn_gse(3)/norm(bn_gse(3:4)) bn_gse(4)/norm(bn_gse(3:4))];
  hl=line([0 -bnproj(3)*25],[0 bnproj(4)*25]);set(hl,'color','red','linewidth',.4);       % B direction
  hl=line([0 -bn_gse(3)*25],[0 bn_gse(4)*25]);set(hl,'color','red','linewidth',2);       % B direction
  text(-30*bnproj(3),30*bnproj(4),'B');           % label
  if flag_v1==1 % plot v vector
    vn1proj=[0 0 vn1_gse(3)/norm(vn1_gse(3:4)) vn1_gse(4)/norm(vn1_gse(3:4))];
    hl=line([0 -vn1proj(3)*25],[0 vn1proj(4)*25]);set(hl,'color','k','linewidth',.4);       % V direction
    hl=line([0 -vn1_gse(3)*25],[0 vn1_gse(4)*25]);set(hl,'color','k','linewidth',2);       % V direction
    text(-30*vn1proj(3),30*vn1proj(4),'V');           % label
  end
  if flag_v2==1 % plot v vector
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
    eval(irf_ssub('line([0 -rp?_gse(2)],[0 rp?_gse(3)]);',ip));
    eval(irf_ssub('patch(-rp?_gse(2)+x_circle*0.4,rp?_gse(3)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip));
    eval(irf_ssub('text(-rp?_gse(2)*.8,rp?_gse(3)*.8,num2str(?));',ip));
  end
  
  axes(h(2));cla
  text(0,50,'(BxS)xS','verticalalignment','middle','horizontalalignment','center','fontweight','demi');
  text(50,0,'BxS','rotation',90,'verticalalignment','middle','horizontalalignment','center','fontweight','demi');
  patch(x_circle*1.5,y_circle*1.5,x_circle*0+1);hold on; % plot spacecraft
  patch(x_circle*1.5,y_circle*1.5,x_circle*0-1);         % plot spacecraft
  for aa=0:pi/12:pi/2
    hl=line(x_circle*44*sin(aa),y_circle*44*sin(aa));
    set(hl,'linestyle',':','color','green','linewidth',.2);
  end
  for ip=1:4
    eval(irf_ssub('line([0 rp?_b(1)],[0 rp?_b(2)]);',ip));
    eval(irf_ssub('patch(rp?_b(1)+x_circle*0.4,rp?_b(2)+y_circle*0.4,x_circle*0+1,''facecolor'',''black'',''edgecolor'',''none'');',ip));
    eval(irf_ssub('text(rp?_b(1)*.8,rp?_b(2)*.8,num2str(?));',ip));
  end
  
  axes(h(4))
  ht=irf_pl_info(['c_efw_pl_orient() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],gca,[0,1 ]); set(ht,'interpreter','none');
  xp=0;yp=.9;dyp=-0.1;
  yp=yp+dyp;
  %if flag_v==1, % add v velocity
  %  vstr=[' V=' num2str(irf_abs(v,1),3) ' [' num2str(irf_norm(v),' %5.2f') '] km/s. '];
  %  text(xp,yp,vstr);
  %  yp=yp+dyp;
  %end
  
else
  disp(sprintf( ...
    'c_efw_pl_orient: action string ''%s'' not recognized, no action taken.',action))
end

