function h=c_pl_sc_conf_xyz(time,coord_sys,flag)
%C_PL_SC_CONF_XYZ   Plot the configuration of CLuster in XYZ coordinates
%
%   h = C_PL_SC_CONF_XYZ;
%   h = C_PL_SC_CONF_XYZ(t);
%   h = C_PL_SC_CONF_XYZ(t,coord_sys);
%   t  - time in isdat epoch
% coord_sys - 'GSE' or 'GSM', default is 'GSE'
% flag - 'default','compact'(to make a compact plot with double axis)
%
% $Id$

%   figuserdata=[h];
eval_figuserdata='figuserdata={h};';

persistent t r1 r2 r3 r4 figNumber coord_label plot_type;
if       (nargin==1 && ischar(time)), action=time;irf_log('fcal',['action=' action]);
elseif   (nargin==3), plot_type=flag;action='initialize';
elseif   (nargin < 9)                   , action='initialize';
end

if strcmp(action,'initialize'),
  if nargin<1, help c_pl_sc_conf_xyz;return;end
  if nargin==1, coord_label='GSE';else coord_label=coord_sys;end
  ok=c_load('R?');
  if  min(ok) == 1,
      switch coord_label
          case {'GSE','gse'}
              c_eval('r?=R?;clear R?;');
          case {'GSM','gsm'}
              c_eval('r?=irf_gse2gsm(R?);clear R?;');
      end
  else
      irf_log('fcal','No position data available');return;
  end
  t=time;

  % See if spacecraft configuration XYZ figure is open
  ch = get(0,'ch');indx=[];
  if ~isempty(ch),
=======
    if nargin<1, help c_pl_sc_conf_xyz;return;end
    if nargin==1, coord_label='GSE';else coord_label=coord_sys;end
    ok=c_load('R?');
    if  min(ok) == 1,
        switch coord_label
            case 'GSE'
                c_eval('r?=R?;clear R?;');
            case 'GSM'
                c_eval('r?=irf_gse2gsm(R?);clear R?;');
            case 'gsm'
                c_eval('coord_label=''GSM'';r?=irf_gse2gsm(R?);clear R?;');
            otherwise
                c_eval('coord_label=''GSE'';r?=R?;clear R?;');
        end
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

    switch plot_type
        case 'default'
            set(figNumber,'Position',[10 10 600 1000])

            h(1)=subplot(4,2,1);axis([-19.99 14.99 -14.99 14.99]);hold on;
            h(2)=subplot(4,2,2);axis([-19.99 19.99 -19.99 19.99]);hold on;
            h(3)=subplot(4,2,3);axis([-19.99 14.99 -19.99 19.99]);hold on;
            h(4)=subplot(4,2,4);axis off;
            h(5)=subplot(4,2,5);axis([-50 50 -50 50]);
            h(6)=subplot(4,2,6);axis([-50 50 -50 50]);
            h(7)=subplot(4,2,7);axis([-50 50 -50 50]);
            h(8)=subplot(4,2,8);axis off;

        case 'compact'
            set(figNumber,'Position',[10 10 700 700]);clear h;

            h(1)=subplot(2,2,1);axis([-50 50 -50 50]);
            h(2)=subplot(2,2,2);axis([-50 50 -50 50]);
            h(3)=subplot(2,2,3);axis([-50 50 -50 50]);
            h(4)=subplot(2,2,4);axis off;

    end

    axes(h(4));
    plot(0,.3,'ks',.2,.3,'rd',.4,.3,'go',.6,.3,'bv','LineWidth',1.5);
    text(0.03,.3,'C1');text(.23,.3,'C2');text(.43,.3,'C3');text(.63,.3,'C4');
    axis off;
    ht=irf_pl_info(['c_pl_sc_conf_xyz() ' datestr(now)],gca,[0,1 ]); set(ht,'interpreter','none');
    htime=irf_pl_info(['Cluster configuration\newline ' epoch2iso(time,1)],gca,[0,.7 ]);set(htime,'fontsize',12);

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
switch plot_type
    case 'default'
        axes(h(1));
        plot(XRe1(2),XRe1(4),'ks', XRe2(2),XRe2(4),'rd', XRe3(2),XRe3(4),'go', XRe4(2),XRe4(4),'bv','LineWidth',1.5);
        xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] '  coord_label]);
        grid on;
        set(gca,'xdir','reverse')
        %  axes(h(1));      title(titlestr);
        axes(h(2));
        plot(XRe1(3),XRe1(4),'ks', XRe2(3),XRe2(4),'rd', XRe3(3),XRe3(4),'go', XRe4(3),XRe4(4),'bv','LineWidth',1.5);
        xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
        grid on;
        axes(h(3));
        plot(XRe1(2),XRe1(3),'ks', XRe2(2),XRe2(3),'rd', XRe3(2),XRe3(3),'go', XRe4(2),XRe4(3),'bv','LineWidth',1.5);
        xlabel(['X [R_E] ' coord_label]);ylabel(['Y [R_E] ' coord_label]);
        grid on;
        set(gca,'xdir','reverse')
        set(gca,'ydir','reverse')


        axes(h(5));
        plot(x1(2),x1(4),'ks', x2(2),x2(4),'rd', x3(2),x3(4),'go', x4(2),x4(4),'bv','LineWidth',1.5);
        xlabel(['X [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
        grid on;axis([-drref drref -drref drref]);
        set(gca,'xdir','reverse')
        %  axes(h(1));      title(titlestr);
        axes(h(6));
        plot(x1(3),x1(4),'ks', x2(3),x2(4),'rd', x3(3),x3(4),'go', x4(3),x4(4),'bv','LineWidth',1.5)
        xlabel(['Y [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
        grid on;axis([-drref drref -drref drref]);
        axes(h(7));
        plot(x1(2),x1(3),'ks', x2(2),x2(3),'rd', x3(2),x3(3),'go', x4(2),x4(3),'bv','LineWidth',1.5)
        xlabel(['X [km] ' coord_label]);ylabel(['Y [km] ' coord_label]);
        grid on;axis([-drref drref -drref drref]);
        set(gca,'xdir','reverse')
        set(gca,'ydir','reverse')
    case 'compact'
        axes(h(1));ax1=gca;
        plot(x1(2),x1(4),'ks', x2(2),x2(4),'rd', x3(2),x3(4),'go', x4(2),x4(4),'bv','LineWidth',1.5);
        xlabel(['{\Delta}X [km] ' coord_label]);ylabel(['{\Delta}Z [km] ' coord_label]);
        grid on;axis([-drref drref -drref drref]);
        set(gca,'xdir','reverse')
        grid on;axis([-drref drref -drref drref]);

        xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
        xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
        xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
        axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
        set(ax1,'position',axpos); % narrow axis
        ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        xlabel(['X [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
        xtlax2=num2str((str2num(xticklabel_ax1)+R(2))/6372,'%5.2f');
        ytlax2=num2str((str2num(yticklabel_ax1)+R(4))/6372,'%5.2f');
        set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);


        %  axes(h(1));      title(titlestr);
        axes(h(2));ax1=gca;
        plot(x1(3),x1(4),'ks', x2(3),x2(4),'rd', x3(3),x3(4),'go', x4(3),x4(4),'bv','LineWidth',1.5)
        xlabel(['{\Delta}Y [km] ' coord_label]);ylabel(['{\Delta}Z [km] ' coord_label]);
        grid on;axis([-drref drref -drref drref]);
        
        xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
        xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
        xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
        axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
        set(ax1,'position',axpos); % narrow axis
        ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        xlabel(['Y [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
        set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
        xtlax2=num2str((str2num(xticklabel_ax1)+R(3))/6372,'%5.2f');
        ytlax2=num2str((str2num(yticklabel_ax1)+R(4))/6372,'%5.2f');
        set(ax1_2,'xticklabel',xtlax2,'yticklabel',ytlax2);

        axes(h(3));ax1=gca;
        plot(x1(2),x1(3),'ks', x2(2),x2(3),'rd', x3(2),x3(3),'go', x4(2),x4(3),'bv','LineWidth',1.5)
        xlabel(['{\Delta}X [km] ' coord_label]);ylabel(['{\Delta}Y [km] ' coord_label]);
        grid on;axis([-drref drref -drref drref]);
        set(gca,'xdir','reverse')
        set(gca,'ydir','reverse')
             
        xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
        xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
        xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
        axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
        set(ax1,'position',axpos); % narrow axis
        ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        xlabel(['X [km] ' coord_label]);ylabel(['Y [km] ' coord_label]);
        set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);set(ax1_2,'xdir','reverse')
        xtlax2=num2str((str2num(xticklabel_ax1)+R(2))/6372,'%5.2f');
        ytlax2=num2str((str2num(yticklabel_ax1)+R(3))/6372,'%5.2f');
        set(ax1_2,'ydir','reverse','xticklabel',xtlax2,'yticklabel',ytlax2);

end
else
    disp(sprintf( ...
        'c_pl_sc_conf_lmn: action string ''%s'' not recognized, no action taken.',action))
end

