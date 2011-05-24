function hout=c_pl_sc_conf_xyz(time,coord_sys,flag,spacecraft)
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
% eval_figuserdata='figuserdata={h};';
cluster_marker={'ks','rd','go','bv'};

persistent t r1 r2 r3 r4 R1 R2 R3 R4 figNumber coord_label plot_type sc_list h flag_show_cluster_description;

if       (nargin==1 && ischar(time)),
    action=time;
    %irf_log('fcal',['action=' action]);
elseif   (nargin==3), plot_type=flag;action='initialize';
elseif   (nargin==4), plot_type=flag;action='initialize';
elseif   (nargin < 9),plot_type='default';action='initialize';
end
if nargin==0, time=[2010 12 31 01 01 01];end
if nargin==4, sc_list=spacecraft;
elseif ~exist('sc_list','var'), sc_list=1:4;
elseif isempty(sc_list), sc_list=1:4;
end
if nargin>=2, % t,coord_sys
    coord_label=upper(coord_sys);
    if ~(strcmp(coord_label,'GSE') || strcmp(coord_label,'GSM')),
        coord_label='GSE'; % default reference frame GSE if does not recognize coord system
    end
end
if exist('coord_label','var'), % define coord label if not defined so far
    if isempty(coord_label),
        coord_label='GSE';
    end
else % in case coord_sys not specified
    coord_label='GSE';
end

switch lower(action)
    case 'initialize' % read in all data and open figure
        if length(time)==1,
            start_time=fromepoch(time);
            t=time;
        elseif length(time)==6,
            start_time=time;
            t=toepoch(time);
        elseif exist('t','var'),
            start_time=fromepoch(t);
        else
            irf_log('fcal','Check time format');return;
        end
        if ~is_R_ok,     % try reading from disk mat files
            ok=c_load('R?',sc_list);
            if ~is_R_ok,  % try reading from isdat server
                c_eval('[tr,r] = irf_isdat_get([''Cluster/?/ephemeris/position''], toepoch(start_time), 60);R?=[tr r];clear tr r;',sc_list);
                if ~is_R_ok,% no idea
                    disp('NO POSITION DATA!');
                end
            end
        end
        %        if ~any(ok),
        %            for ic=sc_list,
        %                [tr,r] = caa_is_get('db.irfu.se:0', toepoch(start_time), 60, ic, 'ephemeris', 'position');
        %                c_eval('R?=[double(tr) double(r)''];',ic);clear tr r;
        %            end
        %        end
        % See if spacecraft configuration XYZ figure is open
        ch = get(0,'ch');indx=[];
        if ~isempty(ch),
            chTags = get(ch,'Tag');
            indx = find(strcmp(chTags,'cplscconfXYZ'));
        end
        if isempty(indx),
            figNumber=figure( ...
                'Name','Cluster s/c configuration in XYZ', ...
                'Tag','cplscconfXYZ');
        else
            figure(ch(indx));clf;figNumber=gcf;
        end
        menus;
        c_pl_sc_conf_xyz(plot_type);
    case 'gse'
        coord_label='GSE';
        c_eval('coord_label=''GSE'';r?=R?;',sc_list);
        c_pl_sc_conf_xyz('plot');
    case 'gsm'
        coord_label='GSM';
        c_eval('r?=irf_gse2gsm(R?);',sc_list);
        c_pl_sc_conf_xyz('plot');
    case 'default'
        set(figNumber,'Position',[10 10 600 1000]);
        delete(findall(gcf,'Type','axes'))
        h=[];
        h(1)=subplot(4,2,1);axis([-19.99 14.99 -14.99 14.99]);hold on;
        h(2)=subplot(4,2,2);axis([-19.99 19.99 -19.99 19.99]);hold on;
        h(3)=subplot(4,2,3);axis([-19.99 14.99 -19.99 19.99]);hold on;
        h(4)=subplot(4,2,4);axis off;
        h(5)=subplot(4,2,5);axis([-50 50 -50 50]);
        h(6)=subplot(4,2,6);axis([-50 50 -50 50]);
        h(7)=subplot(4,2,7);axis([-50 50 -50 50]);
        h(8)=subplot(4,2,8);axis off;
        flag_show_cluster_description=1; % show cluster description
        plot_type='default';
        c_pl_sc_conf_xyz(coord_label);
        
    case 'compact'
        set(figNumber,'Position',[10 10 700 700]);
        flag_show_cluster_description=1; % show cluster description
        plot_type='compact';
        c_pl_sc_conf_xyz(coord_label);
        
    case 'supercompact'
        set(figNumber,'Position',[10 10 700 350]);
        flag_show_cluster_description=0;
        plot_type='supercompact';
        c_pl_sc_conf_xyz(coord_label);
        
    case 'supercompact2'
        set(figNumber,'Position',[10 10 350 650]);
        flag_show_cluster_description=0;
        plot_type='supercompact2';
        c_pl_sc_conf_xyz(coord_label);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'plot'
        c_eval('rr?=irf_resamp(r?,t);',sc_list);
        R=0; c_eval('R=R+rr?/length(sc_list);',sc_list);
        c_eval('XRe?=irf_tappl(rr?,''/6372'');dr?=rr?-R;dr?(1)=t;dr?=irf_abs(dr?);x?=dr?;',sc_list);
        drref=0; c_eval('drref=max([drref dr?(5)]);',sc_list);
        if drref==0, drref=1; end % in case 1 satellite or satellites in the same location:)
        %%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
        switch plot_type
            case 'default'
                cla(h(1));
                c_eval('plot(h(1),XRe?(2),XRe?(4),cluster_marker{?},''LineWidth'',1.5);hold(h(1),''on'');',sc_list);
                xlabel(h(1),['X [R_E] ' coord_label]);
                ylabel(h(1),['Z [R_E] '  coord_label]);
                grid(h(1),'on')
                set(h(1),'xdir','reverse')
                
                %  axes(h(1));      title(titlestr);
                cla(h(2));
                c_eval('plot(h(2),XRe?(3),XRe?(4),cluster_marker{?},''LineWidth'',1.5);hold(h(2),''on'');',sc_list);
                xlabel(h(2),['Y [R_E] ' coord_label]);
                ylabel(h(2),['Z [R_E] ' coord_label]);
                grid(h(2),'on');
                
                cla(h(3));
                c_eval('plot(h(3),XRe?(2),XRe?(3),cluster_marker{?},''LineWidth'',1.5);hold(h(3),''on'');',sc_list);
                xlabel(h(3),['X [R_E] ' coord_label]);
                ylabel(h(3),['Y [R_E] ' coord_label]);
                grid(h(3),'on');
                set(h(3),'xdir','reverse')
                set(h(3),'ydir','reverse')
                
                cla(h(5));
                c_eval('plot(h(5),x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(h(5),''on'');',sc_list);
                xlabel(h(5),['X [km] ' coord_label]);
                ylabel(h(5),['Z [km] ' coord_label]);
                grid(h(5),'on');
                axis(h(5),[-drref drref -drref drref]);
                set(h(5),'xdir','reverse')
                
                cla(h(6));
                c_eval('plot(h(6),x?(3),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(h(6),''on'');',sc_list);
                xlabel(h(6),['Y [km] ' coord_label]);
                ylabel(h(6),['Z [km] ' coord_label]);
                grid(h(6),'on');
                axis(h(6),[-drref drref -drref drref]);
                
                cla(h(7));
                c_eval('plot(h(7),x?(2),x?(3),cluster_marker{?},''LineWidth'',1.5);hold(h(7),''on'');',sc_list);
                xlabel(h(7),['X [km] ' coord_label]);
                ylabel(h(7),['Y [km] ' coord_label]);
                grid(h(7),'on');
                axis(h(7),[-drref drref -drref drref]);
                set(h(7),'xdir','reverse')
                set(h(7),'ydir','reverse')
            case 'compact'
                delete(findall(gcf,'Type','axes'))
                h=[];
                h(1)=subplot(2,2,1);axis([-50 50 -50 50]);
                h(2)=subplot(2,2,2);axis([-50 50 -50 50]);
                h(3)=subplot(2,2,3);axis([-50 50 -50 50]);
                h(4)=subplot(2,2,4);axis off;
                hold(h(1),'off');
                c_eval('plot(h(1),x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(h(1),''on'');',sc_list);
                xlabel(h(1),['{\Delta}X [km] ' coord_label]);
                ylabel(h(1),['{\Delta}Z [km] ' coord_label]);
                set(h(1),'xdir','reverse')
                grid(h(1),'on');
                axis(h(1),[-drref drref -drref drref]);
                if drref>1000, REform='%6.1f';
                elseif drref<100, REform='%6.3f';
                else REform='%6.2f';
                end
                
                axpos=get(h(1),'position');
                axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(h(1),'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(h(1),'XLim');ylim_ax1=get(h(1),'YLim');
                xtick_ax1=get(h(1),'XTick');ytick_ax1=get(h(1),'YTick');
                xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                %  axes(h(1));      title(titlestr);
                cla(h(2));ax1=h(2);
                c_eval('plot(h(2),x?(3),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(h(2),''on'');',sc_list);
                xlabel(h(2),['{\Delta}Y [km] ' coord_label]);
                ylabel(h(2),['{\Delta}Z [km] ' coord_label]);
                grid(h(2),'on');
                axis(h(2),[-drref drref -drref drref]);
                
                axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(ax1,'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
                set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
                xtlax2=num2str((xtick_ax1'+R(3))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                cla(h(3));ax1=h(3);
                c_eval('plot(ax1,x?(2),x?(3),cluster_marker{?},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);ylabel(ax1,['{\Delta}Y [km] ' coord_label]);
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                set(ax1,'xdir','reverse')
                set(ax1,'ydir','reverse')
                
                axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(ax1,'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(['X [R_E] ' coord_label]);ylabel(['Y [R_E] ' coord_label]);
                set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);set(ax1_2,'xdir','reverse')
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(3))/6372,REform);
                set(ax1_2,'ydir','reverse','xticklabel',xtlax2,'yticklabel',ytlax2);
                
            case 'supercompact'
                delete(findall(gcf,'Type','axes'))
                h=[];
                set(gcf,'PaperUnits','centimeters')
                xSize = 16; ySize = 8;
                xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
                set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
                
                h(1)=axes('position',[0.09 0.18 0.33 0.66]); % [x y dx dy]
                h(2)=axes('position',[0.57 0.18 0.33 0.66]); % [x y dx dy]
                h(3)=axes('position',[0 0 1 1]);axis off;
                %h(2)=subplot(1,2,2);axis([-50 50 -50 50]);
                
                ax1=h(1);
                c_eval('plot(ax1,x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Z [km] ' coord_label]);
                set(ax1,'xdir','reverse')
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                if drref>10000, REform='%6.1f';
                elseif drref<100, REform='%6.3f';
                else REform='%6.2f';
                end
                
                axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(ax1,'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                %  axes(h(1));      title(titlestr);
                ax1=h(2);
                c_eval('plot(ax1,x?(3),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}Y [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Z [km] ' coord_label]);
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                c_eval('text(x?(3),x?(4),''  C?'',''parent'',ax1,''HorizontalAlignment'',''Left'');',sc_list);
                text(0.02,0.94,epoch2iso(t,1),'parent',ax1,'units','normalized','horizontalalignment','left','parent',h(1),'fontsize',9);
                
                axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(ax1,'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(['Y [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
                set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
                xtlax2=num2str((xtick_ax1'+R(3))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xticklabel',xtlax2,'yticklabel',ytlax2);
            case 'supercompact2'
                delete(findall(gcf,'Type','axes'))
                h=[];
                set(gcf,'PaperUnits','centimeters')
                xSize = 16; ySize = 8;
                xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
                set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
                
                h(2)=axes('position',[0.18 0.08 0.66 0.4]); % [x y dx dy]
                h(1)=axes('position',[0.18 0.58 0.66 0.4]); % [x y dx dy]
                h(3)=axes('position',[0 0 1 1]);axis off;
                %h(2)=subplot(1,2,2);axis([-50 50 -50 50]);
                
                ax1=h(1);
                c_eval('plot(ax1,x?(2),x?(4),cluster_marker{?},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Z [km] ' coord_label]);
                set(ax1,'xdir','reverse')
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                if drref>10000, REform='%6.1f';
                elseif drref<100, REform='%6.3f';
                else REform='%6.2f';
                end
                
                axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(ax1,'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(['X [R_E] ' coord_label]);ylabel(['Z [R_E] ' coord_label]);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                ax1=h(2);
                c_eval('plot(ax1,x?(2),x?(3),cluster_marker{?},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Y [km] ' coord_label]);
                set(ax1,'xdir','reverse','ydir','reverse')
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                c_eval('text(x?(2),x?(3),''  C?'',''parent'',ax1,''HorizontalAlignment'',''Left'');',sc_list);
                text(0.02,0.94,epoch2iso(t,1),'parent',ax1,'units','normalized','horizontalalignment','left','parent',h(1),'fontsize',9);
                
                axpos=get(ax1,'position');axpos(1)=axpos(1)+0.02;axpos(3)=axpos(3)-0.05;axpos(4)=axpos(4)-0.05;
                set(ax1,'position',axpos); % narrow axis
                ax1_2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(['X [R_E] ' coord_label]);ylabel(['Y [R_E] ' coord_label]);
                set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(3))/6372,REform);
                set(ax1_2,'xdir','reverse','ydir','reverse','xticklabel',xtlax2,'yticklabel',ytlax2);
        end
        if flag_show_cluster_description==1,
            plot(h(4),0,.3,'ks',.2,.3,'rd',.4,.3,'go',.6,.3,'bv','LineWidth',1.5);
            text(0.03,.3,'C1','parent',h(4));
            text(.23,.3,'C2','parent',h(4));
            text(.43,.3,'C3','parent',h(4));
            text(.63,.3,'C4','parent',h(4));
            axis(h(4),'off');
            ht=irf_pl_info(['c_pl_sc_conf_xyz() ' datestr(now)],h(4),[0,1 ]);
            set(ht,'interpreter','none');
            htime=irf_pl_info(['Cluster configuration\newline ' epoch2iso(t,1)],h(4),[0,.7 ]);
            set(htime,'fontsize',12);
        end
        
    case 'new_time'
        xx=inputdlg('Enter new time. [yyyy mm dd hh mm ss]','**',1,{mat2str(fromepoch(t))});
        variable_str=xx{1};
        t=toepoch(eval(variable_str));
        c_pl_sc_conf_xyz('initialize');
    case 'new_sc_list'
        xx=inputdlg('Enter new sc_list. ex. [1 3 4]','**',1,{mat2str(sc_list)});
        variable_str=xx{1};
        sc_list=eval(variable_str);
        c_pl_sc_conf_xyz('initialize');
end
if nargout,
    hout=h;
else
    clear hout;
end

function menus
% generate menus
if isempty(findobj(gcf,'type','uimenu','label','&Options'))
    hcoordfigmenu=uimenu('Label','&Options');
    uimenu(hcoordfigmenu,'Label','New time','Callback','c_pl_sc_conf_xyz(''new_time'')')
    uimenu(hcoordfigmenu,'Label','&GSE','Callback','c_pl_sc_conf_xyz(''GSE'')','Accelerator','G')
    uimenu(hcoordfigmenu,'Label','GS&M','Callback','c_pl_sc_conf_xyz(''GSM'')','Accelerator','M')
    uimenu(hcoordfigmenu,'Label','normal','Callback','c_pl_sc_conf_xyz(''default'')')
    uimenu(hcoordfigmenu,'Label','compact','Callback','c_pl_sc_conf_xyz(''compact'')')
    uimenu(hcoordfigmenu,'Label','supercompact','Callback','c_pl_sc_conf_xyz(''supercompact'')')
    uimenu(hcoordfigmenu,'Label','supercompact2','Callback','c_pl_sc_conf_xyz(''supercompact2'')')
    uimenu(hcoordfigmenu,'Label','New sc_list','Callback','c_pl_sc_conf_xyz(''new_sc_list'')')
    user_data = get(gcf,'userdata');
    user_data.coordfigmenu=1;
    set(gcf,'userdata',user_data);
end

function answer=is_R_ok
sc_list=evalin('caller','sc_list');
t=evalin('caller','t');
for ic=1:numel(sc_list)
    stric=num2str(sc_list(ic));
    if evalin('caller',['numel(R' stric ')']) < 8 % less than 2 time points
        answer=0;
        return;
    else
        tint=evalin('caller',['[ R' stric '(1,1) R' stric '(end,1)]']);
        if (tint(1)>t) || (tint(2)<t),
            answer=0;
            return;
        end
    end
end
answer=1;