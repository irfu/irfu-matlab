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
cluster_marker={{'ks','markersize',12},{'rd','markersize',12},...
  {'go','markersize',12,'color',[0 0.6 0]},{'bv','markersize',12}};
cluster_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
  {'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
cluster_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
  {'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};

if       (nargin==1 && ischar(time)),
    action=time;
    %irf_log('fcal',['action=' action]);
elseif   (nargin==3), plot_type=flag;action='initialize';
elseif   (nargin==4), plot_type=flag;action='initialize';
elseif   (nargin < 9),plot_type='default';action='initialize';
end
if nargin==0, % default time (with time can make smarter solution)
    if evalin('caller','exist(''tint'')'),
        time=irf_time(evalin('caller','tint(1)'),'vector');
    elseif exist('CAA','dir')
        R=irf_get_data('sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M','caa','mat');
        if numel(R)==0,
            time=[2010 12 31 01 01 01];
        else
            time=0.5*(R(1,1)+R(end,1)); % first point in center of position time series
        end
    else
        time=[2010 12 31 01 01 01];
    end
end
if nargin==4, sc_list=spacecraft;
else sc_list=1:4;
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
        if length(time)==1, % time given in epoch
            t=time;
        elseif length(time)==6, % time given as vector
            t=irf_time(time);
        else
            irf_log('fcal','Wrong input format of time.');
            return;
        end
        % Open new figure
        figNumber=figure( ...
            'Name','Cluster s/c configuration in XYZ', ...
            'Tag','cplscconfXYZ');
        set(figNumber,'defaultLineLineWidth', 1.5);
        set(figNumber,'defaultAxesFontSize', 12);
        set(figNumber,'defaultTextFontSize', 12);
        set(figNumber,'defaultAxesFontUnits', 'pixels');
        menus;
        data.t=t;
        data.figNumber=figNumber;
        data.coord_label=coord_label;
        data.plot_type=plot_type;
        data.sc_list=sc_list;
        c_eval('data.r?=[];data.R?=[];');
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz('read_position');
        c_pl_sc_conf_xyz(plot_type);

    case 'read_position'
        data=get(gcf,'userdata');
        c_eval('R?=data.R?;',data.sc_list);
        if ~is_R_ok,     % try reading from disk mat files
            c_load('R?',sc_list);
        end
        if ~is_R_ok,     % try reading from CAA files
            irf_log('dsrc','Trying to read CAA files...')
            c_eval('R?=irf_get_data(''sc_r_xyz_gse__C?_CP_AUX_POSGSE_1M'',''caa'',''mat'');',sc_list);
        end
        if ~is_R_ok,  % try reading from isdat server
            irf_log('dsrc','Trying to obtain satellite position from isdat server...')
            try
                c_eval('[tr,r] = irf_isdat_get([''Cluster/?/ephemeris/position''], data.t, 60);R?=[tr r];clear tr r;',data.sc_list);
                if ~is_R_ok,% no idea
                    disp('NO POSITION DATA!');
                end
            catch ME
                irf_log('dsrc',['Did not succeed! (' ME.identifier ')'] );
            end
        end
        if ~is_R_ok,     % could not obtain
            irf_log('dsrc','Could not obtain position data!')
            c_eval('R?=[];',data.sc_list);
        end
        c_eval('data.R?=R?;',data.sc_list);
        if isfield(data,'omni')
            if data.t < data.omni(1,1) || data.t > data.omni(end,1)
                omni=irf_get_data(data.t+[-2 2]*3600,'p,bx,bygsm,bzgsm','omni');
            end
        else
            omni=irf_get_data(data.t+[-2 2]*3600,'p,bx,bygsm,bzgsm','omni');
            data.omni=omni;
        end
        set(gcf,'userdata',data);
        return;
        
    case 'gse'
        data=get(gcf,'userdata');
        data.coord_label='GSE';
        c_eval('data.r?=data.R?;',data.sc_list);
        set(gcf,'userdata',data);
        if strcmp(data.plot_type,'lmn'), % need to redraw lmn text
            c_pl_sc_conf_xyz('lmn');
        else
            c_pl_sc_conf_xyz('plot');
        end
        
    case 'gsm'
        data=get(gcf,'userdata');
        data.coord_label='GSM';
        c_eval('data.r?=irf_gse2gsm(data.R?);',data.sc_list);
        set(gcf,'userdata',data);
        if strcmp(data.plot_type,'lmn'), % need to redraw lmn text
            c_pl_sc_conf_xyz('lmn');
        else
            c_pl_sc_conf_xyz('plot');
        end
        
    case 'default'
        data=get(gcf,'userdata');
        ss=get(0,'screensize');
        sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
        set(gcf,'Position',[10 10 600/sfactor 1000/sfactor]);
        delete(findall(gcf,'Type','axes'))
        data.h=[];h=zeros(1,8);
        xsize=.35;ysize=.195;dx=.13;dy=.05;
        for ix=1:2,
            for iy=1:4,
                h(iy*2-2+ix)=axes('position',[dx*ix+(ix-1)*xsize dy*(5-iy)+(4-iy)*ysize xsize ysize]);                
            end
        end
        axis(h(8),'off');
        axis(h(1:3),[-20 20 -20 20]);
        axis(h(4),[-20 20 0 20]);
        for ii=1:4, hold(h(ii),'on');daspect(h(ii),[1 1 1]);end
        data.h=h;
        data.flag_show_cluster_description=1; % show cluster description
        data.plot_type='default';
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz(data.coord_label);
        
    case 'compact'
        data=get(gcf,'userdata');
        clf;menus;
        set(gcf,'userdata',data);
        ss=get(0,'screensize');
        sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
        set(gcf,'Position',[10 ss(4)-80-700/sfactor 700/sfactor 700/sfactor]);
        initialize_figure;
        h=[];
        h(1)=axes('position',[0.1  0.56 0.3 0.36]); % [x y dx dy]
        h(2)=axes('position',[0.59 0.56 0.3 0.36]); % [x y dx dy]
        h(3)=axes('position',[0.1  0.06 0.3 0.36]); % [x y dx dy]
        h(4)=axes('position',[0.59 0.06 0.3 0.36]); % [x y dx dy]
        h(21) = axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        h(22) = axes('Position',get(h(2),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        h(23) = axes('Position',get(h(3),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        axis(h(4),'off');hold(h(4),'on');
        data.h=h;
        data.flag_show_cluster_description=1; % show cluster description
        data.plot_type='compact';
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz(data.coord_label);
    case 'config3d'
        data=get(gcf,'userdata');
        clf;menus;
        set(gcf,'userdata',data);
        ss=get(0,'screensize');
        sfactor=max([1 700/(ss(3)-80) 1000/(ss(4)-80)]);
        set(gcf,'Position',[10 ss(4)-80-500/sfactor 500/sfactor 500/sfactor]);
        initialize_figure;
        h=[];
        h(1)=axes('position',[0.15  0.16 0.7 0.7]); % [x y dx dy]
        data.h=h;
        data.flag_show_cluster_description=1; % show cluster description
        data.plot_type='config3d';
        data.flag_show_cluster_description=0;
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz(data.coord_label);
    case 'lmn'
        c_pl_sc_conf_xyz('compact');
        data=get(gcf,'userdata');
        delete(data.h(21)); % remove secondary axis
        delete(data.h(22));
        delete(data.h(23));
        data.h(5:end)=[];
        data.plot_type='lmn';
        hca=data.h(4);
        cla(hca);hold(hca,'on');
        % L vector
        callbackStr='c_pl_sc_conf_xyz(''plot'')';
        data.LMN_text_hndl=uicontrol('string',['LMN vectors in ' data.coord_label '. One of L/M/N can be zero.'],'style','text','units','normalized','Position',[0.5 0.25 .35 .05]);
        uicontrol('string','L','style','text','units','normalized','Position',[0.5 0.2 .05 .04])
        if isfield(data,'Lstr'), Lstr=data.Lstr;else Lstr='[1 0 0]';end
        data.L_hndl=uicontrol('Style','edit','Units','normalized', ...
            'Position',[0.55 0.2 .3 .05],'String',Lstr,'Callback',callbackStr);
        uicontrol('string','M','style','text','units','normalized','Position',[0.5 0.15 .05 .04])
        if isfield(data,'Lstr'), Mstr=data.Mstr;else Mstr='[0 1 0]';end
        data.M_hndl=uicontrol('Style','edit','Units','normalized', ...
            'Position',[0.55 0.15 .3 .05],'String',Mstr,'Callback',callbackStr);
        uicontrol('string','N','style','text','units','normalized','Position',[0.5 0.1 .05 .04])
        if isfield(data,'Nstr'), Nstr=data.Nstr;else Nstr='0';end
        data.N_hndl=uicontrol('Style','edit','Units','normalized', ...
            'Position',[0.55 0.1 .3 .05],'String',Nstr,'Callback',callbackStr);
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz('plot');
        
    case 'supercompact'
        data=get(gcf,'userdata');
        ss=get(0,'screensize');
        sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
        set(gcf,'Position',[10 ss(4)-80-350/sfactor 750/sfactor 350/sfactor]);
        initialize_figure;
        h=[];
        h(1)=axes('position',[0.09 0.13 0.32 0.74]); % [x y dx dy]
        h(2)=axes('position',[0.59 0.13 0.32 0.74]); % [x y dx dy]
        h(3)=axes('position',[0 0 1 1]);
        h(21)=axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        h(22)=axes('Position',get(h(2),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        axis(h(3),'off');
        data.h=h;
        data.flag_show_cluster_description=0;
        data.plot_type='supercompact';
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz(data.coord_label);
        
    case 'supercompact2'
        data=get(gcf,'userdata');
        ss=get(0,'screensize');
        sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
        set(gcf,'Position',[10 ss(4)-80-650/sfactor 350/sfactor 650/sfactor]);
        initialize_figure;
        h=[];
        h(2)=axes('position',[0.18 0.06 0.63 0.36]); % [x y dx dy]
        h(1)=axes('position',[0.18 0.57 0.63 0.36]); % [x y dx dy]
        h(3)=axes('position',[0 0 1 1]);axis off;
        h(21)=axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        h(22)=axes('Position',get(h(2),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
        axis(h(3),'off');
        data.h=h;
        data.flag_show_cluster_description=0;
        data.plot_type='supercompact2';
        set(gcf,'userdata',data);
        c_pl_sc_conf_xyz(data.coord_label);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'plot'
        data=get(gcf,'userdata');
        c_eval('rr?=irf_resamp(data.r?,data.t);',data.sc_list);
        R=0; c_eval('R=R+rr?/length(data.sc_list);',data.sc_list);
        c_eval('XRe?=irf_tappl(rr?,''/6372'');dr?=rr?-R;dr?(1)=data.t;dr?=irf_abs(dr?);x?=dr?;',data.sc_list);
        x=[];
        for ic=1:numel(data.sc_list),
            eval(['x{ic}=x' num2str(data.sc_list(ic)) ';']);
        end
        drref=0; c_eval('drref=max([drref dr?(5)]);',data.sc_list);
        if drref==0, drref=1; end % in case 1 satellite or satellites in the same location:)
        set(gcf,'userdata',data);
        %%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
        h=data.h;
        coord_label=data.coord_label;
        sc_list=data.sc_list;
        switch data.plot_type
            case 'default'
                cla(h(1));
                c_eval('plot(h(1),XRe?(2),XRe?(4),cluster_marker_small{?}{:});hold(h(1),''on'');',sc_list);
                xlabel(h(1),['X [R_E] ' coord_label]);
                ylabel(h(1),['Z [R_E] '  coord_label]);
                grid(h(1),'on')
                set(h(1),'xdir','reverse')
                add_magnetopause(h(1));
                add_bowshock(h(1));
                add_Earth(h(1));
                
                %  axes(h(1));      title(titlestr);
                cla(h(2));
                c_eval('plot(h(2),XRe?(3),XRe?(4),cluster_marker_small{?}{:});hold(h(2),''on'');',sc_list);
                xlabel(h(2),['Y [R_E] ' coord_label]);
                ylabel(h(2),['Z [R_E] ' coord_label]);
                grid(h(2),'on');
                 
                cla(h(3));
                c_eval('plot(h(3),XRe?(2),XRe?(3),cluster_marker_small{?}{:});hold(h(3),''on'');',sc_list);
                xlabel(h(3),['X [R_E] ' coord_label]);
                ylabel(h(3),['Y [R_E] ' coord_label]);
                grid(h(3),'on');
                set(h(3),'xdir','reverse')
                set(h(3),'ydir','reverse')
                add_magnetopause(h(3));
                add_bowshock(h(3));
                add_Earth(h(3));
                
                cla(h(4));
                c_eval('plot(h(4),XRe?(2),sqrt(XRe?(3)^2+XRe?(4)^2),cluster_marker_small{?}{:});hold(h(4),''on'');',sc_list);
                xlabel(h(4),['X [R_E] ' coord_label]);
                ylabel(h(4),['sqrt (Y^2+Z^2) [R_E] ' coord_label]);
                grid(h(4),'on');
                set(h(4),'xdir','reverse')
                add_magnetopause(h(4));
                add_bowshock(h(4));
                add_Earth(h(4));
                
                cla(h(5));
                c_eval('plot(h(5),x?(2),x?(4),cluster_marker{?}{:});hold(h(5),''on'');',sc_list);
                xlabel(h(5),['X [km] ' coord_label]);
                ylabel(h(5),['Z [km] ' coord_label]);
                grid(h(5),'on');
                axis(h(5),[-drref drref -drref drref]);
                set(h(5),'xdir','reverse')
                
                cla(h(6));
                c_eval('plot(h(6),x?(3),x?(4),cluster_marker{?}{:});hold(h(6),''on'');',sc_list);
                xlabel(h(6),['Y [km] ' coord_label]);
                ylabel(h(6),['Z [km] ' coord_label]);
                grid(h(6),'on');
                axis(h(6),[-drref drref -drref drref]);
                
                cla(h(7));
                c_eval('plot(h(7),x?(2),x?(3),cluster_marker{?}{:});hold(h(7),''on'');',sc_list);
                xlabel(h(7),['X [km] ' coord_label]);
                ylabel(h(7),['Y [km] ' coord_label]);
                grid(h(7),'on');
                axis(h(7),[-drref drref -drref drref]);
                set(h(7),'xdir','reverse')
                set(h(7),'ydir','reverse')

            case 'compact'
                
                hold(h(1),'off');
                c_eval('plot(h(1),x?(2),x?(4),cluster_marker{?}{:});hold(h(1),''on'');',sc_list);
                xlabel(h(1),['{\Delta}X [km] ' coord_label]);
                ylabel(h(1),['{\Delta}Z [km] ' coord_label]);
                set(h(1),'xdir','reverse');
                grid(h(1),'on');
                axis(h(1),[-drref drref -drref drref]);
                if drref>1000, REform='%6.1f';
                elseif drref<100, REform='%6.3f';
                else REform='%6.2f';
                end
                
                xlim_ax1=get(h(1),'XLim');ylim_ax1=get(h(1),'YLim');
                xtick_ax1=get(h(1),'XTick');ytick_ax1=get(h(1),'YTick');
                xlabel(h(21),['X [R_E] ' coord_label]);
                ylabel(h(21),['Z [R_E] ' coord_label]);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(h(21),'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                cla(h(2));
                c_eval('plot(h(2),x?(3),x?(4),cluster_marker{?}{:});hold(h(2),''on'');',sc_list);
                xlabel(h(2),['{\Delta}Y [km] ' coord_label]);
                ylabel(h(2),['{\Delta}Z [km] ' coord_label]);
                grid(h(2),'on');
                axis(h(2),[-drref drref -drref drref]);
                
                xlim_ax1=get(h(2),'XLim');ylim_ax1=get(h(2),'YLim');
                xtick_ax1=get(h(2),'XTick');ytick_ax1=get(h(2),'YTick');
                xlabel(h(22),['Y [R_E] ' coord_label]);
                ylabel(h(22),['Z [R_E] ' coord_label]);
                set(h(22),'xlim',xlim_ax1,'ylim',ylim_ax1);
                xtlax2=num2str((xtick_ax1'+R(3))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(h(22),'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                cla(h(3));
                c_eval('plot(h(3),x?(2),x?(3),cluster_marker{?}{:});hold(h(3),''on'');',sc_list);
                xlabel(h(3),['{\Delta}X [km] ' coord_label]);
                ylabel(h(3),['{\Delta}Y [km] ' coord_label]);
                grid(h(3),'on');
                axis(h(3),[-drref drref -drref drref]);
                set(h(3),'xdir','reverse')
                set(h(3),'ydir','reverse')
                
                xlim_ax1=get(h(3),'XLim');ylim_ax1=get(h(3),'YLim');
                xtick_ax1=get(h(3),'XTick');ytick_ax1=get(h(3),'YTick');
                xlabel(h(23),['X [R_E] ' coord_label]);
                ylabel(h(23),['Y [R_E] ' coord_label]);
                set(h(23),'xlim',xlim_ax1,'ylim',ylim_ax1,'xdir','reverse')
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(3))/6372,REform);
                set(h(23),'ydir','reverse','xticklabel',xtlax2,'yticklabel',ytlax2);
            case 'config3d'
                
                hold(h(1),'off');
                c_eval('plot3(h(1),-drref,x?(3),x?(4),cluster_marker_shaded{?}{:});hold(h(1),''on'');',sc_list);
                c_eval('plot3(h(1),x?(2),-drref,x?(4),cluster_marker_shaded{?}{:});hold(h(1),''on'');',sc_list);
                c_eval('plot3(h(1),x?(2),x?(3),-drref,cluster_marker_shaded{?}{:});hold(h(1),''on'');',sc_list);
                axis(h(1),[-drref drref -drref drref -drref drref ]);
                for ii=1:4, 
                  for jj=ii+1:4,
                    if any(find(sc_list==ii)) && any(find(sc_list==jj)),
                      c_eval('line([x?(2) x!(2)],[x?(3) x!(3)],[x?(4) x!(4)],''parent'',h(1),''linewidth'',2,''linestyle'',''-'',''color'',[0.6 0.6 0.6])',ii,jj);
                    end
                  end
                end
                c_eval('line([x?(2) -drref],[x?(3) x?(3)],[x?(4) x?(4)],''parent'',h(1),''linestyle'','':'',''linewidth'',0.6);',sc_list);
                c_eval('line([x?(2) x?(2)],[x?(3) -drref],[x?(4) x?(4)],''parent'',h(1),''linestyle'','':'',''linewidth'',0.6);',sc_list);
                c_eval('line([x?(2) x?(2)],[x?(3) x?(3)],[x?(4) -drref],''parent'',h(1),''linestyle'','':'',''linewidth'',0.6);',sc_list);
                line([-drref drref],[-drref -drref],[-drref -drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
                line([-drref -drref],[-drref drref],[-drref -drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
                line([-drref -drref],[-drref -drref],[-drref drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
                c_eval('plot3(h(1),x?(2),x?(3),x?(4),cluster_marker{?}{:});hold(h(1),''on'');',sc_list);
                text(0.1,1,0,irf_time(data.t,'isoshort'),'parent',h(1),'units','normalized','horizontalalignment','center','fontsize',9);
                xlabel(h(1),['{\Delta}X [km] ' coord_label]);
                ylabel(h(1),['{\Delta}Y [km] ' coord_label]);
                zlabel(h(1),['{\Delta}Z [km] ' coord_label]);
                set(h(1),'xdir','reverse');
                set(h(1),'ydir','reverse');
                grid(h(1),'on');
                axis(h(1),[-drref drref -drref drref]);
                
            case 'lmn'
                cla(h(1));
                x=get_in_lmn(x);
                for ic=1:numel(sc_list);
                    plot(h(1),x{ic}(2),x{ic}(4),cluster_marker{sc_list(ic)}{:});
                    hold(h(1),'on');
                end
                xlabel(h(1),'L [km]');
                ylabel(h(1),'N [km] ');
                set(h(1),'xdir','reverse');
                grid(h(1),'on');
                axis(h(1),[-drref drref -drref drref]);
                              
                cla(h(2));
                for ic=1:numel(sc_list);
                    plot(h(2),x{ic}(3),x{ic}(4),cluster_marker{sc_list(ic)}{:});
                    hold(h(2),'on');
                end
                xlabel(h(2),'M [km] ');
                ylabel(h(2),'N [km] ');
                grid(h(2),'on');
                axis(h(2),[-drref drref -drref drref]);
                
                cla(h(3));
                for ic=1:numel(sc_list);
                    plot(h(3),x{ic}(2),x{ic}(3),cluster_marker{sc_list(ic)}{:});
                    hold(h(3),'on');
                end
                xlabel(h(3),'L [km] ');
                ylabel(h(3),'M [km] ');
                grid(h(3),'on');
                axis(h(3),[-drref drref -drref drref]);
                set(h(3),'xdir','reverse')
                set(h(3),'ydir','reverse')
                                
            case 'supercompact'
                ax1=h(1);ax1_2=h(21);
                hold(ax1,'off');
                c_eval('plot(ax1,x?(2),x?(4),cluster_marker{?}{:},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Z [km] ' coord_label]);
                set(ax1,'xdir','reverse')
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                if drref>10000, REform='%6.1f';
                elseif drref<100, REform='%6.3f';
                else REform='%6.2f';
                end
                
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(ax1_2,['X [R_E] ' coord_label]);ylabel(ax1_2,['Z [R_E] ' coord_label]);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                %  axes(h(1));      title(titlestr);
                ax1=h(2);ax1_2=h(22);
                hold(ax1,'off');
                c_eval('plot(ax1,x?(3),x?(4),cluster_marker{?}{:});hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}Y [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Z [km] ' coord_label]);
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                c_eval('text(x?(3),x?(4),''   C?'',''parent'',ax1,''HorizontalAlignment'',''Left'');',sc_list);
                text(0.02,0.94,irf_time(data.t,'isoshort'),'parent',ax1,'units','normalized','horizontalalignment','left','parent',h(1),'fontsize',9);
                
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(ax1_2,['Y [R_E] ' coord_label]);ylabel(ax1_2,['Z [R_E] ' coord_label]);
                set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
                xtlax2=num2str((xtick_ax1'+R(3))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xticklabel',xtlax2,'yticklabel',ytlax2);

            case 'supercompact2'
                
                ax1=h(1);ax1_2=h(21);
                hold(ax1,'off');
                c_eval('plot(ax1,x?(2),x?(4),cluster_marker{?}{:},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Z [km] ' coord_label]);
                set(ax1,'xdir','reverse')
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                if drref>10000, REform='%6.1f';
                elseif drref<100, REform='%6.3f';
                else REform='%6.2f';
                end
                
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(ax1_2,['X [R_E] ' coord_label]);ylabel(ax1_2,['Z [R_E] ' coord_label]);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(4))/6372,REform);
                set(ax1_2,'xdir','reverse','xlim',xlim_ax1,'ylim',ylim_ax1,'xticklabel',xtlax2,'yticklabel',ytlax2);
                
                ax1=h(2);ax1_2=h(22);
                hold(ax1,'off');
                c_eval('plot(ax1,x?(2),x?(3),cluster_marker{?}{:},''LineWidth'',1.5);hold(ax1,''on'');',sc_list);
                xlabel(ax1,['{\Delta}X [km] ' coord_label]);
                ylabel(ax1,['{\Delta}Y [km] ' coord_label]);
                set(ax1,'xdir','reverse','ydir','reverse')
                grid(ax1,'on');
                axis(ax1,[-drref drref -drref drref]);
                c_eval('text(x?(2),x?(3),''   C?'',''parent'',ax1,''HorizontalAlignment'',''Left'');',sc_list);
                text(0.02,0.94,irf_time(data.t,'isoshort'),'parent',ax1,'units','normalized','horizontalalignment','left','parent',h(1),'fontsize',9);
                
                xlim_ax1=get(ax1,'XLim');ylim_ax1=get(ax1,'YLim');
                xtick_ax1=get(ax1,'XTick');ytick_ax1=get(ax1,'YTick');
                xlabel(ax1_2,['X [R_E] ' coord_label]);ylabel(ax1_2,['Y [R_E] ' coord_label]);
                set(ax1_2,'xlim',xlim_ax1,'ylim',ylim_ax1);
                xtlax2=num2str((xtick_ax1'+R(2))/6372,REform);
                ytlax2=num2str((ytick_ax1'+R(3))/6372,REform);
                set(ax1_2,'xdir','reverse','ydir','reverse','xticklabel',xtlax2,'yticklabel',ytlax2);
        end
        if data.flag_show_cluster_description==1,
            if strcmpi(data.plot_type,'compact') || ...
                    strcmpi(data.plot_type,'lmn')  % decide in which axes write labels
                hca=h(4);
            else
                hca=h(8);
            end
            cla(hca);
            axis(hca,[0 1 0 1]);
            hold(hca,'on');
            yy=1;
            plot(hca,0,yy,'ks',.2,yy,'rd',.4,yy,'go',.6,yy,'bv','LineWidth',1.5);
            text(0.03,yy,'C1','parent',hca);
            text(.23,yy,'C2','parent',hca);
            text(.43,yy,'C3','parent',hca);
            text(.63,yy,'C4','parent',hca);
            axis(hca,'off');
            ht=irf_legend(hca,['c_pl_sc_conf_xyz() ' datestr(now)],[0,0],'fontsize',8);
            set(ht,'interpreter','none');
            htime=irf_legend(hca,['Cluster configuration\newline ' irf_time(data.t,'isoshort')],[0,.95]);
            set(htime,'fontsize',12);
            if isfield(data,'omni')
                if ~isempty(data.omni) && ~strcmpi(data.plot_type,'lmn'),
                    fft=irf_resamp(data.omni,data.t);
                    irf_legend(hca,['IMF from OMNI 1h database:\newline P=' num2str(fft(2),'%6.1f') '[nPa],\newline Bx=' num2str(fft(3),'%6.1f') ',By=' num2str(fft(4),'%6.1f') ',Bz=' num2str(fft(5),'%6.1f') '[nT] GSM' ],[0,0.7]);
                end
            end
        end

    case 'new_time'
        data=get(gcf,'userdata');
        xx=inputdlg('Enter new time. [yyyy mm dd hh mm ss]','**',1,{mat2str(irf_time(data.t,'vector'))});
        if ~isempty(xx),
            variable_str=xx{1};
            data.t=irf_time(eval(variable_str));
            set(gcf,'userdata',data);
            c_pl_sc_conf_xyz('read_position');
            c_pl_sc_conf_xyz(data.coord_label);
        end
    case 'new_sc_list'
        xx=inputdlg('Enter new sc_list. ex. [1 3 4]','**',1,{mat2str(sc_list)});
        if ~isempty(xx),
            variable_str=xx{1};
            sc_list=eval(variable_str);
            data=get(gcf,'userdata');
            data.sc_list=sc_list;
            set(gcf,'userdata',data);
            c_pl_sc_conf_xyz('read_position');
            c_pl_sc_conf_xyz(data.coord_label);
        end
end
if nargout,
    hout=h;
else
    clear hout;
end

function initialize_figure
delete(findall(gcf,'Type','axes'))
set(gcf,'PaperUnits','centimeters')
pos=get(gcf,'position');
factor=min(16./pos(3:4));
xSize = pos(3)*factor; ySize = pos(4)*factor; % should not exceed 16cm in x or y direction
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'defaultLineLineWidth', 1.5);
set(gcf,'defaultAxesFontUnits', 'pixels');
set(gcf,'defaultAxesFontSize', 12);


function menus
% generate menus
if isempty(findobj(gcf,'type','uimenu','label','&Options'))
    hcoordfigmenu=uimenu('Label','&Options');
    uimenu(hcoordfigmenu,'Label','New time','Callback','c_pl_sc_conf_xyz(''new_time'')')
    uimenu(hcoordfigmenu,'Label','&GSE','Callback','c_pl_sc_conf_xyz(''GSE'')','Accelerator','G')
    uimenu(hcoordfigmenu,'Label','GS&M','Callback','c_pl_sc_conf_xyz(''GSM'')','Accelerator','M')
    uimenu(hcoordfigmenu,'Label','normal','Callback','c_pl_sc_conf_xyz(''default'')')
    uimenu(hcoordfigmenu,'Label','config3D','Callback','c_pl_sc_conf_xyz(''config3D'')')
    uimenu(hcoordfigmenu,'Label','compact','Callback','c_pl_sc_conf_xyz(''compact'')')
    uimenu(hcoordfigmenu,'Label','supercompact','Callback','c_pl_sc_conf_xyz(''supercompact'')')
    uimenu(hcoordfigmenu,'Label','supercompact2','Callback','c_pl_sc_conf_xyz(''supercompact2'')')
    uimenu(hcoordfigmenu,'Label','LMN','Callback','c_pl_sc_conf_xyz(''lmn'')')
    uimenu(hcoordfigmenu,'Label','New sc_list','Callback','c_pl_sc_conf_xyz(''new_sc_list'')')
    user_data = get(gcf,'userdata');
    user_data.coordfigmenu=1;
    set(gcf,'userdata',user_data);
end

function answer=is_R_ok
data=get(gcf,'userdata');
sc_list=data.sc_list;
t=data.t;
for ic=1:numel(sc_list)
    stric=num2str(sc_list(ic));
    if evalin('caller',['numel(R' stric ')']) < 8 % less than 2 time points
        answer=0;
        return;
    else
        tint=evalin('caller',['[' 'R' stric '(1,1) R' stric '(end,1)]']);
        if (tint(1)>t) || (tint(2)<t),
            answer=0;
            return;
        end
    end
end
answer=1;

function add_magnetopause(h)
t=getfield(get(gcf,'userdata'),'t');
[x,y]=irf_magnetosphere('mp_shue1998',t);
x=[fliplr(x) x];
y=[fliplr(y) -y];
line(x,y,'parent',h,'linewidth',0.5,'linestyle','-','color','k');

function add_bowshock(h)
t=getfield(get(gcf,'userdata'),'t');
[x,y]=irf_magnetosphere('bs',t);
x=[fliplr(x) x];
y=[fliplr(y) -y];
line(x,y,'parent',h,'linewidth',0.5,'linestyle','-','color','k');

function add_Earth(h)
theta=0:pi/20:pi;
x=sin(theta);y=cos(theta);
patch(-x,y,'k','edgecolor','none','parent',h)
patch(x,y,'w','edgecolor','k','parent',h)

function y=get_in_lmn(x)
% get x in LMN reference frame
data=get(gcf,'userdata');
y=x; % intialize
Lstr=get(data.L_hndl, 'string');
Mstr=get(data.M_hndl, 'string');
Nstr=get(data.N_hndl, 'string');
L=eval(Lstr);
M=eval(Mstr);
N=eval(Nstr);
if numel(L)~=3, L=0; end
if numel(M)~=3, M=0; end
if numel(N)~=3, N=0; end
if numel(L)+numel(M)+numel(N) <= 5 % not enough information
    irf_log('fcal','LMN not correctly defined');
    return
end
data.Lstr=Lstr;
data.Mstr=Mstr;
data.Nstr=Nstr;
for ic=1:numel(data.sc_list),
    y{ic}=irf_newxyz(x{ic},L,M,N);
end
if L==0, data.Lstr=['[' num2str(irf_norm(cross(M,N)),'%6.2f') ']']; set(data.L_hndl,'string',data.Lstr);end
if M==0, data.Mstr=['[' num2str(irf_norm(cross(N,L)),'%6.2f') ']']; set(data.M_hndl,'string',data.Mstr);end
if N==0, data.Nstr=['[' num2str(irf_norm(cross(L,M)),'%6.2f') ']']; set(data.N_hndl,'string',data.Nstr);end
set(gcf,'userdata',data);
