function h=c_pl_sc_configuration_LMN(time,R1,R2,R3,R4,L,M,N);
% C_PL_SC_CONFIGURATION_LMN plots the configuration of CLuster in LMN coordinates
%   h = C_PL_SC_CONFIGURATION_LMN;
%   h = C_PL_SC_CONFIGURATION_LMN(t);
%   h = C_PL_SC_CONFIGURATION_LMN(t,n);
%   h = C_PL_SC_CONFIGURATION_LMN(t,n,b);
%   h = C_PL_SC_CONFIGURATION_LMN(t,l,m,n);
%   h = C_PL_SC_CONFIGURATION_LMN(t,R1,R2,R3,R4,n);
%   ic - spacecraft number
%   t  - time in isdat epoch
%   R1 - position vectors, columns are [t X Y Z]
%   n  - the direction of normal in GSE [nx ny nz]
%   b  - magnetic field in GSE [bx by bz]
%   n  - N direction in LMN reference system (N along normal, outward;L closest to B;M=LxN)
%   l  - L
%   m  - M

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'c_pl_sc_conf_lmn')

%   figuserdata=[h];
eval_figuserdata='figuserdata={h};';
LMN_Ltitle='L along B, N closest to normal, M=LxN';
LMN_Ntitle='N along normal, L closest to B, M=LxN';

persistent t b l m n B r1 r2 r3 r4 phaseHndl timeHndl figNumber ...
            resHndl NHndl LHndl Lflag LMN_Lflag LMN_Nflag ...
            flag_v1 flag_v2 v1 v2;
if       (nargin==1 & isstr(time)), action=time;disp(['action=' action]);
elseif   (nargin < 9)                   , action='initialize';
end

if strcmp(action,'initialize'),
  initLflag=0;
  if   exist('mB.mat'), load mB B3; B=B3;
  else                  disp('No magnetic field data available, reading from DDS');B=c_get_bfgm(time+[0 1],3);
  end

  if nargin<1, help c_pl_sc_configuration_LMN;return;                                      end
  if nargin==1, n=[0 0 1]; b=av_interp(B,time); l=b(1,[2 3 4]);                                                                    end
  if nargin==2, n=R1;      b=av_interp(B,time); l=b(1,[2 3 4]);                                   end
  if nargin==3,
    n=R1;
    if size(R2,1)>1, B=R2;
    else,
      initLflag=1;
      if size(R2,2)>3 l=R2(2:4);
      else, l=R2(1:3);
      end;
    end
  end
  if nargin==4,
    l=R1;m=R2;n=R3;b=av_interp(B,time); 
    initLflag=1;
  end
  if nargin<5,  % load position vectors if not given
    if   exist('mR.mat'), load mR R1 R2 R3 R4;for ic=1:4,eval(av_ssub('r?=R?;clear R?;',ic)),end
    else                  disp('No position data available');return;
    end
  end
  if nargin==5,
    r1=R1;r2=R2;r3=R3;r4=R4;
    l=[1 0 0];n=[0 0 1];
  end
  if nargin==6,
    r1=R1;r2=R2;r3=R3;r4=R4;
    n=L;l=B(1,[2 3 4])
  end
  t=time;

  % See if spacecraft configuration LMN figure is open
  ch = get(0,'ch');indx=[];
  if ~isempty(ch),
        chTags = get(ch,'Tag');
        indx = find(strcmp(chTags,'cplscconfLMN'));
  end
  if isempty(indx),
  figNumber=figure( ...
        'Name',['Cluster s/c configuration in LMN'], ...
        'Tag','cplscconfLMN');
  else
        figure(ch(indx));clf;figNumber=gcf;
  end

  h(1)=subplot(2,2,1);axis([-50 50 -50 50]);
  h(2)=subplot(2,2,2);axis([-50 50 -50 50]);
  h(3)=subplot(2,2,3);axis([-50 50 -50 50]);
  h(4)=subplot(2,2,4);axis off;
  	ht=av_pl_info(['c_pl_sc_configuration_LMN() ' datestr(now)],gca,[0,1 ]); set(ht,'interpreter','none');
    resHndl=text(0,0.8,'result');
  %====================================
  % The normal entering
  labelStr=num2str(n);
  callbackStr='c_pl_sc_configuration_LMN(''plot'')';
  uicontrol('style','text','units','normalized','Position',[0.5 0.2 .2 .05],'string','N, normal vector [GSE]','Callback',callbackStr)
  NHndl=uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'Position',[0.7 0.2 .2 .05], ...
        'String',labelStr, ...
        'Callback',callbackStr);
  %====================================
  % The magnetic field or L entering
  labelStr=num2str(l);
  callbackStr='c_pl_sc_configuration_LMN(''plot'')';
  Lflag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.25 .2 .05],'string','L (or B) [GSE]','Callback',callbackStr);
  if initLflag==1, set(Lflag,'value',1);end
  LHndl=uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'Position',[0.7 0.25 .2 .05], ...
        'String',labelStr, ...
        'Callback',callbackStr);
  %====================================
  % The time entering
  labelStr=[datestr(datenum(fromepoch(t))) '.' num2str(floor(mod(t,1)*100),'%.2d')];
  callbackStr='c_pl_sc_configuration_LMN(''plot'')';
  uicontrol('style','text','units','normalized','Position',[0.5 0.15 .2 .05],'string','Time')
  timeHndl=uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'Position',[0.7 0.15 .2 .05], ...
        'String',labelStr, ...
        'Callback',callbackStr);
  %====================================
  % The definition of LMN coordinates
  callbackStr='c_pl_sc_configuration_LMN(''fixLMN_N'')';
  LMN_Lflag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.10 .3 .05],'string',LMN_Ltitle,'Callback',callbackStr);
  callbackStr='c_pl_sc_configuration_LMN(''fixLMN_L'')';
  LMN_Nflag=uicontrol('style','checkbox','units','normalized','Position',[0.5 0.05 .3 .05],'string',LMN_Ntitle,'Callback',callbackStr);
  set(LMN_Lflag,'value',1);

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
  c_pl_sc_configuration_LMN('plot');


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%% action fix LMN %%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 elseif strcmp(action,'fixLMN_N'),
    if get(LMN_Lflag,'value'),
      set(LMN_Nflag,'value',0);
    else
      set(LMN_Nflag,'value',1);
    end
    c_pl_sc_configuration_LMN('plot')
 elseif strcmp(action,'fixLMN_L'),
    if get(LMN_Nflag,'value'),
      set(LMN_Lflag,'value',0);
    else
      set(LMN_Lflag,'value',1);
    end
    c_pl_sc_configuration_LMN('plot')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(action,'plot'),
  figuserdata=get(figNumber,'userdata');
  h=figuserdata{1};

  t=toepoch(datevec(get(timeHndl, 'string')));

  flag_L=get(Lflag, 'value');
  if flag_L==0,
    b=av_interp(B,t); l=b(1,[2 3 4]);set(LHndl,'string',num2str(l));
  else
    l=eval(['[ ' get(LHndl,'string') ' ]']);
  end
  n=eval(['[ ' get(NHndl,'string') ' ]']);

  if get(LMN_Nflag,'value'),
    n=n./norm(n);
    l=cross(n,cross(l,n));l=l./norm(l);
    m=cross(l,n);
    titlestr=LMN_Ntitle;
  else
    l=l./norm(l);
    n=cross(l,cross(n,l));n=n./norm(n);
    m=cross(l,n);
    titlestr=LMN_Ltitle;
  end

  for ic=1:4,eval(av_ssub('rr?=av_interp(r?,t);',ic)),end
  R=(rr1+rr2+rr3+rr4)/4;
  for ic=1:4,eval(av_ssub('dr?=rr?-R;dr?(1)=t;dr?=av_abs(dr?);drlnm?=av_newxyz(dr?,l,n,m);x?=drlnm?;',ic)),end
  drref=max([dr1(5) dr2(5) dr3(5) dr4(5)]);

	%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
  axes(h(1));
    plot(x1(3),x1(2),'ks', x2(3),x2(2),'rd', x3(3),x3(2),'go', x4(3),x4(2),'bv')
    xlabel('N [km]');ylabel('L [km]');
    grid on;axis([-drref drref -drref drref]);
  axes(h(1));      title(titlestr);
  axes(h(2));
    plot(x1(4),x1(2),'ks', x2(4),x2(2),'rd', x3(4),x3(2),'go', x4(4),x4(2),'bv')
    xlabel('M [km]');ylabel('L [km]');
    grid on;axis([-drref drref -drref drref]);
  axes(h(3));
    plot(x1(3),x1(4),'ks', x2(3),x2(4),'rd', x3(3),x3(4),'go', x4(3),x4(4),'bv')
    xlabel('N [km]');ylabel('M [km]');
    grid on;axis([-drref drref -drref drref]);
	axes(h(4));
    L_str=[' L=[' num2str(l,'%6.2f') ']'];
    if exist('b'),
      L_str=[L_str ', angle between L and B is ' num2str(av_angle(l,b,1),2) ' deg'];
    end
    N_str=['\newline N=[' num2str(n,'%6.2f') ']'];
    M_str=['\newline M=[' num2str(m,'%6.2f') ']'];
    set(resHndl,'string',[L_str N_str M_str],'verticalalignment','top');


else
  disp(sprintf( ...
     'c_pl_sc_configuration_LMN: action string ''%s'' not recognized, no action taken.',action))
end

