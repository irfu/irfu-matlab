function h=c_pl_sc_conf_xyz(time,coord_sys);
%C_PL_SC_CONF_XYZ   Plot the configuration of CLuster in XYZ coordinates
%
%   h = C_PL_SC_CONF_XYZ;
%   h = C_PL_SC_CONF_XYZ(t);
%   h = C_PL_SC_CONF_XYZ(t,coord_sys);
%   t  - time in isdat epoch
% coord_sys - 'GSE' or 'GSM', default is 'GSE'
%
% $Id$

%   figuserdata=[h];
eval_figuserdata='figuserdata={h};';

persistent t b l m n B r1 r2 r3 r4 phaseHndl timeHndl figNumber ...
            resHndl NHndl LHndl Lflag LMN_Lflag LMN_Nflag ...
            coord_label;
if       (nargin==1 & isstr(time)), action=time;irf_log('fcal',['action=' action]);
elseif   (nargin < 9)                   , action='initialize';
end

if strcmp(action,'initialize'),
  initLflag=0;
  if nargin<1, help c_pl_sc_conf_xyz;return;end
  if nargin==1, coord_label='GSE';else, coord_label=coord_sys;end
  ok=c_load('R?');
  if  min(ok) == 1,
      c_eval('r?=R?;clear R?;');
  else
      irf_log('fcal','No position data available');return;
  end
  t=time;

  % See if spacecraft configuration XYZ figure is open
  ch = get(0,'ch');indx=[];
  if ~isempty(ch),
        chTags = get(ch,'Tag');
        indx = find(strcmp(chTags,'cplscconfXYZ'));
  end
  if isempty(indx),
  figNumber=figure( ...
        'Name',['Cluster s/c configuration in XYZ'], ...
        'Tag','cplscconfXYZ');
  else
        figure(ch(indx));clf;figNumber=gcf;
  end
  set(figNumber,'Position',[10 10 600 1000])

  h(1)=subplot(4,2,1);axis([-19.99 9.99 -14.99 14.99]);hold on;
  h(2)=subplot(4,2,2);axis([-19.99 19.99 -19.99 19.99]);hold on;
  h(3)=subplot(4,2,3);axis([-19.99 9.99 -19.99 19.99]);hold on;
  h(4)=subplot(4,2,4);axis off;
  h(5)=subplot(4,2,5);axis([-50 50 -50 50]);
  h(6)=subplot(4,2,6);axis([-50 50 -50 50]);
  h(7)=subplot(4,2,7);axis([-50 50 -50 50]);
  h(8)=subplot(4,2,8);axis off;

  axes(h(4));
  ht=irf_pl_info(['c_pl_sc_conf_xyz() ' datestr(now)],gca,[0,1 ]); set(ht,'interpreter','none');
  htime=irf_pl_info(['Cluster configuration\newline ' epoch2iso(time,1)],gca,[0,.5 ]);set(htime,'fontsize',12); 


 eval(eval_figuserdata);
 set(figNumber,'UserData',figuserdata);
 c_pl_sc_conf_xyz('plot');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(action,'plot'),
  figuserdata=get(figNumber,'userdata');
  h=figuserdata{1};

  c_eval('rr?=irf_resamp(r?,t);');
  R=(rr1+rr2+rr3+rr4)/4;
  c_eval('XRe?=irf_tappl(rr?,''/6372'');dr?=rr?-R;dr?(1)=t;dr?=irf_abs(dr?);x?=dr?;');
  drref=max([dr1(5) dr2(5) dr3(5) dr4(5)]);

	%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
  axes(h(1));
    plot(XRe1(2),XRe1(4),'ks', XRe2(2),XRe2(4),'rd', XRe3(2),XRe3(4),'go', XRe4(2),XRe4(4),'bv');
    xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] '  coord_label]);
    grid on;
    set(gca,'xdir','reverse')
%  axes(h(1));      title(titlestr);
  axes(h(2));
    plot(XRe1(3),XRe1(4),'ks', XRe2(3),XRe2(4),'rd', XRe3(3),XRe3(4),'go', XRe4(3),XRe4(4),'bv');
    xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
    grid on;
  axes(h(3));
    plot(XRe1(2),XRe1(3),'ks', XRe2(2),XRe2(3),'rd', XRe3(2),XRe3(3),'go', XRe4(2),XRe4(3),'bv');
    xlabel(['X [R_E] ' coord_label]);ylabel(['Y [R_E] ' coord_label]);
    grid on;
    set(gca,'xdir','reverse')
    
  axes(h(5));
    plot(x1(3),x1(2),'ks', x2(3),x2(2),'rd', x3(3),x3(2),'go', x4(3),x4(2),'bv');
    xlabel(['X [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
    grid on;axis([-drref drref -drref drref]);
    set(gca,'xdir','reverse')
%  axes(h(1));      title(titlestr);
  axes(h(6));
    plot(x1(4),x1(2),'ks', x2(4),x2(2),'rd', x3(4),x3(2),'go', x4(4),x4(2),'bv')
    xlabel(['Y [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
    grid on;axis([-drref drref -drref drref]);
  axes(h(7));
    plot(x1(3),x1(4),'ks', x2(3),x2(4),'rd', x3(3),x3(4),'go', x4(3),x4(4),'bv')
    xlabel(['X [km] ' coord_label]);ylabel(['Y [km] ' coord_label]);
    grid on;axis([-drref drref -drref drref]);
    set(gca,'xdir','reverse')    

else
  disp(sprintf( ...
     'c_pl_sc_conf_lmn: action string ''%s'' not recognized, no action taken.',action))
end

