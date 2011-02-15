function h=c_pl_sc_conf_xyz(time,coord_sys,flag,spacecraft)
%C_PL_SC_CONF_XYZ   Plot the configuration of CLuster in XYZ coordinates
%
%   h = C_PL_SC_CONF_XYZ;
%   h = C_PL_SC_CONF_XYZ(t);
%   h = C_PL_SC_CONF_XYZ(t,coord_sys);
%   h = C_PL_SC_CONF_XYZ(t,coord_sys,flag);
%   h = C_PL_SC_CONF_XYZ(t,coord_sys,flag,spacecraft);
%   t  - time in isdat epoch or vector [year  month day hour min seconds]
% coord_sys - 'GSE' or 'GSM', default is 'GSE'
% flag - 'default','compact','supercompact'
% spacecraft - default 1:4, in other case specify
% $Id$

%   figuserdata=[h];
eval_figuserdata='figuserdata={h};';
cluster_marker={'ks','rd','go','bv'};

persistent t r1 r2 r3 r4 figNumber coord_label plot_type sc_list;
if       (nargin==1 && ischar(time)), action=time;irf_log('fcal',['action=' action]);
elseif   (nargin==3), plot_type=flag;action='initialize';
elseif   (nargin==4), plot_type=flag;action='initialize';
elseif   (nargin < 9),plot_type='default';action='initialize';
end
if nargin==4, sc_list=spacecraft;
elseif ~exist('sc_list') sc_list=1:4;
else 
    if isempty(sc_list), sc_list=1:4;end
end


if strcmp(action,'initialize'),
    if length(time)==1,
        start_time=fromepoch(time);
    elseif length(time)==6,
        start_time=time;
        time=toepoch(time);
    else
        irf_log('fcal','Check time format');return;
    end
    if nargin<1, help c_pl_sc_conf_xyz;return;end
    if nargin==1, coord_label='GSE';else coord_label=coord_sys;end
    ok=c_load('R?',sc_list);
    if ~any(ok),
        DATABASE = c_ctl(0,'isdat_db');
        db = Mat_DbOpen(DATABASE);
        for ic=1:4,
            [tr,r] = isGetDataLite( db, start_time, 60,'Cluster', num2str(ic), 'ephemeris', 'position', ' ', ' ', ' ');
            c_eval('R?=[double(tr) double(r)''];',ic);clear tr r;
        end
    end
    switch coord_label
        case 'GSE'
            c_eval('r?=R?;clear R?;',sc_list);
        case 'GSM'
            c_eval('r?=irf_gse2gsm(R?);clear R?;',sc_list);
        case 'gsm'
            c_eval('coord_label=''GSM'';r?=irf_gse2gsm(R?);clear R?;',sc_list);
        otherwise
            c_eval('coord_label=''GSE'';r?=R?;clear R?;',sc_list);
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
    
    flag_show_cluster_description=1; % per default make Cluster description labels
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
        case 'supercompact'
            set(figNumber,'Position',[10 10 700 350]);clear h;
            set(gcf,'PaperUnits','centimeters')
            xSize = 16; ySize = 8;
            xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
            set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
            
            h(1)=axes('position',[0.09 0.18 0.33 0.66]); % [x y dx dy]
            h(2)=axes('position',[0.57 0.18 0.33 0.66]); % [x y dx dy]
            h(3)=axes('position',[0 0 1 1]);axis off;
            %h(2)=subplot(1,2,2);axis([-50 50 -50 50]);
            flag_show_cluster_description=0;
    end
    
    if flag_show_cluster_description==1,
    axes(h(4));
    plot(0,.3,'ks',.2,.3,'rd',.4,.3,'go',.6,.3,'bv','LineWidth',1.5);
    text(0.03,.3,'C1');text(.23,.3,'C2');text(.43,.3,'C3');text(.63,.3,'C4');
    axis off;
    ht=irf_pl_info(['c_pl_sc_conf_xyz() ' datestr(now)],gca,[0,1 ]); set(ht,'interpreter','none');
    htime=irf_pl_info(['Cluster configuration\newline ' epoch2iso(time,1)],gca,[0,.7 ]);set(htime,'fontsize',12);
    end
    
    eval(eval_figuserdata);
    set(figNumber,'UserData',figuserdata);
    c_pl_sc_conf_xyz('plot');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(action,'plot'),
    
    figuserdata=get(figNumber,'userdata');
    h=figuserdata{1};
    
    c_eval('rr?=irf_resamp(r?,t);',sc_list);
    R=0; c_eval('R=R+rr?/length(sc_list);',sc_list);
    c_eval('XRe?=irf_tappl(rr?,''/6372'');dr?=rr?-R;dr?(1)=t;dr?=irf_abs(dr?);x?=dr?;',sc_list);
    drref=0; c_eval('drref=max([drref dr?(5)]);',sc_list);
    
    %%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
    switch plot_type
        case 'default'
            axes(h(1));
            c_eval('plot(XRe?(2),XRe?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] '  coord_label]);
            grid on;
            set(gca,'xdir','reverse')
            %  axes(h(1));      title(titlestr);
            axes(h(2));
            c_eval('plot(XRe?(3),XRe?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
            grid on;
            axes(h(3));
            c_eval('plot(XRe?(2),XRe?(3),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['X [R_E] ' coord_label]);ylabel(['Y [R_E] ' coord_label]);
            grid on;
            set(gca,'xdir','reverse')
            set(gca,'ydir','reverse')
            
            
            axes(h(5));
            c_eval('plot(x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['X [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
            grid on;axis([-drref drref -drref drref]);
            set(gca,'xdir','reverse')
            %  axes(h(1));      title(titlestr);
            axes(h(6));
            c_eval('plot(x?(3),x?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['Y [km] ' coord_label]);ylabel(['Z [km] ' coord_label]);
            grid on;axis([-drref drref -drref drref]);
            axes(h(7));
            c_eval('plot(x?(2),x?(3),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['X [km] ' coord_label]);ylabel(['Y [km] ' coord_label]);
            grid on;axis([-drref drref -drref drref]);
            set(gca,'xdir','reverse')
            set(gca,'ydir','reverse')
        case 'compact'
            axes(h(1));ax1=gca;
            c_eval('plot(x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['{\Delta}X [km] ' coord_label]);ylabel(['{\Delta}Z [km] ' coord_label]);
            set(gca,'xdir','reverse')
            grid on;axis([-drref drref -drref drref]);
            if drref>1000, REform='%6.1f';
            elseif drref<100, REform='%6.3f';
            else REform='%6.2f';
            end
            
            axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
            set(ax1,'position',axpos); % narrow axis
            ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
            xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
            xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
            xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
            xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
            xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
            ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
            set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                        
            %  axes(h(1));      title(titlestr);
            axes(h(2));ax1=gca;
            c_eval('plot(x?(3),x?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['{\Delta}Y [km] ' coord_label]);ylabel(['{\Delta}Z [km] ' coord_label]);
            grid on;axis([-drref drref -drref drref]);
            
            axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
            set(ax1,'position',axpos); % narrow axis
            ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
            xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
            xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
            xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
            xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
            set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
            xtlax2=num2str((xtick_ax1'+R(3))/6372,REform);
            ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
            set(ax1_2,'xticklabel',xtlax2,'yticklabel',ytlax2);
            
            axes(h(3));ax1=gca;
            c_eval('plot(x?(2),x?(3),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['{\Delta}X [km] ' coord_label]);ylabel(['{\Delta}Y [km] ' coord_label]);
            grid on;axis([-drref drref -drref drref]);
            set(gca,'xdir','reverse')
            set(gca,'ydir','reverse')
            
            axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
            set(ax1,'position',axpos); % narrow axis
            ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
            xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
            xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
            xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
            xlabel(['X [R_E] ' coord_label]);ylabel(['Y [R_E] ' coord_label]);
            set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);set(ax1_2,'xdir','reverse')
            xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
            ytlax2=num2str((ytick_ax1'+R(3))/6372,REform);
            set(ax1_2,'ydir','reverse','xticklabel',xtlax2,'yticklabel',ytlax2);
            
        case 'supercompact'
            axes(h(1));ax1=gca;
            c_eval('plot(x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['{\Delta}X [km] ' coord_label]);ylabel(['{\Delta}Z [km] ' coord_label]);
            set(gca,'xdir','reverse')
            grid on;axis([-drref drref -drref drref]);
            if drref>10000, REform='%6.1f';
            elseif drref<100, REform='%6.3f';
            else REform='%6.2f';
            end
            
            axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
            set(ax1,'position',axpos); % narrow axis
            ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
            xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
            xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
            xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
            xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
            xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
            ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
            set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                        
            %  axes(h(1));      title(titlestr);
            axes(h(2));ax1=gca;
            c_eval('plot(x?(3),x?(4),cluster_marker{?},''LineWidth'',1.5);hold on;',sc_list);
            xlabel(['{\Delta}Y [km] ' coord_label]);ylabel(['{\Delta}Z [km] ' coord_label]);
            grid on;axis([-drref drref -drref drref]);
            c_eval('text(x?(3),x?(4),''  C?'',''HorizontalAlignment'',''Left'');',sc_list);
            text(0.02,0.94,epoch2iso(t,1),'units','normalized','horizontalalignment','left','parent',h(1),'fontsize',9);
            
            axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
            set(ax1,'position',axpos); % narrow axis
            ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
            xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
            xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
            xticklabel_ax1=get(ax1,'XTickLabel');yticklabel_ax1=get(ax1,'YTickLabel');
            xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
            set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
            xtlax2=num2str((xtick_ax1'+R(3))/6372,REform);
            ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
            set(ax1_2,'xticklabel',xtlax2,'yticklabel',ytlax2);


    end
else
    disp(sprintf( ...
        'c_pl_sc_conf_lmn: action string ''%s'' not recognized, no action taken.',action))
end

