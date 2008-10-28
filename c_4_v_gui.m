function c_4_v_gui(x1,x2,x3,x4,column)
%C_4_V_GUI interactive discontinuity analyzer for Cluster
%
%  status = c_4_v_gui(x1,x2,x3,x4,column)
%  status = c_4_v_gui('B?',column)   use variables B1,B2,B3,B4 from base workspace
%  status = c_4_v_gui('B?')          use column 2
%
% $Id$

persistent ud ic var1 var2 var3 var4 var_col variable_str variable_label
tlim = []; 

if       (nargin<=2 && ischar(x1)), % either action as parameter or string variable
    if findstr(x1,'?'),
        if nargin == 1, var_col=2;else var_col=x2;end;
        variable_str=x1;
        action='new_var';
    else
        action=x1;
    end
elseif   (nargin ==4),
    irf_log('fcal','Using second column');
    var1=x1;var2=x2;var3=x3;var4=x4;var_col=2;
    action='initialize';
elseif   (nargin == 5),
    var1=x1;var2=x2;var3=x3;var4=x4;var_col=column;
    action='initialize';
else      help c_4_v_gui
end

irf_log('fcal',['action=' action]);
switch action,
    case {'c1','c2','c3','c4','c5','c6'}
        var_col=str2num(action(2:end));
        c_4_v_gui('update_var_col');
    case 'update_var_col'
        axes(ud.h(1));
        xl=get(ud.h(1),'XLim');
        c_pl_tx(var1,var2,var3,var4,var_col);zoom on;
        ylabel(var_label(variable_str,var_col));
        axis([xl(1) xl(2) 0 1]);
        set(ud.h(1),'YLimMode','auto');add_timeaxis
        axes(ud.h(2));
        c_pl_tx(var1,var2,var3,var4,var_col);
        c_4_v_gui('dt');
    case 'new_var_enter'
        xx=inputdlg('Enter new variable mask. Examples: B? or R? or P?p1','**',1,{'B?'});
        variable_str=xx{1};
        c_4_v_gui('new_var');
    case 'new_var'
        evalin('base',['if ~exist(''' irf_ssub(variable_str,1) '''), c_load(''' variable_str ''');end' ]);
        var1=evalin('base',irf_ssub(variable_str,1));
        var2=evalin('base',irf_ssub(variable_str,2));
        var3=evalin('base',irf_ssub(variable_str,3));
        var4=evalin('base',irf_ssub(variable_str,4));
        if var_col > size(var1,2), var_col=2;end % in case new variable has less columns
        if isempty(ud),
            c_4_v_gui('initialize');
        else
            if ishandle(ud.h(1)),
                for j_col=2:size(var1,2)
                    if j_col<=length(ud.hcol),
                        if ishandle(ud.hcol(j_col)),
                            set(ud.hcol(j_col),'enable','on')
                        else
                            eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''c_4_v_gui(''''c' num2str(j_col) ''''')'');'];
                            eval(eval_str);
                        end
                    else
                        eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''c_4_v_gui(''''c' num2str(j_col) ''''')'');'];
                        eval(eval_str);
                    end
                    for j_col=(size(var1,2)+1):length(ud.hcol)
                        set(ud.hcol(j_col),'enable','off')
                    end
                end
                c_4_v_gui('update_var_col');
            else
                c_4_v_gui('initialize');
            end
        end
    case 'initialize'
        dgh=figure;clf;irf_figmenu;
        h(1)=subplot(3,1,1);
        c_pl_tx(var1,var2,var3,var4,var_col);zoom on;
        irf_pl_info(['c\_4\_v\_int() ' datestr(now)]); % add information to the plot
        ylabel(var_label(variable_str,var_col));

        h(2)=subplot(3,1,2);
        c_pl_tx(var1,var2,var3,var4,var_col);
        ylabel(var_label(variable_str,var_col));

        hh=h(1,1);  % use the first subplot to estimate available time interval
        xl=get(hh,'XLim');
        hc=get(hh,'Children');
        xd=get(hc(end),'XData');
        avail=[min([xl xd]) max([xl xd])];
        presel=xl;

        dt = 0.02*diff(avail);
        xlim = [avail(1)-dt avail(2)+dt];
        ttics = timeaxis(xlim);

        ud.tlim = avail;
        ud.h=h;

        xp=0.05;yp=0.25;
        uch1 = uicontrol('style', 'text', 'string', '[dt1 dt2 dt3 dt4] =','units','normalized','position', [xp yp 0.15 0.03]);
        ud.dt = uicontrol('style', 'edit', ...
            'string', '[0 0 0 0]', ...
            'callback', 'c_4_v_gui(''dt'')', ...
            'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.29 0.05]);

        xp=0.05;yp=0.2;
        uch1 = uicontrol('style', 'text', 'string', '[vx vy vz] km/s =','units','normalized','position', [xp yp 0.15 0.03]);
        ud.v = uicontrol('style', 'edit', ...
            'string', '0*[0 0 0]', ...
            'callback', 'c_4_v_gui(''v'')', ...
            'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.29 0.05]);

        xp=0.05;yp=0.15;
        uch1 = uicontrol('style', 'text', 'string', 'Low pass filter f/Fs = ','units','normalized','position', [xp yp 0.15 0.03]);
        ud.filter = uicontrol('style', 'edit', ...
            'string', '1', ...
            'callback', 'c_4_v_gui(''dt'')', ...
            'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.1 0.05]);

        xp=0.05;yp=0.10;
        ud.coord_sys = uicontrol('style', 'checkbox', ...
            'string', 'velocity in GSM', ...
            'callback', 'c_4_v_gui(''dt'')', ...
            'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.3 0.05]);

        xp=0.05;yp=0.05;
        uch1 = uicontrol('style', 'text', 'string', 'Reference satellite ','units','normalized','position', [xp yp 0.15 0.03]);
        ud.ref_satellite = uicontrol('style', 'edit', ...
            'string', '1', ...
            'callback', 'c_4_v_gui(''v'')', ...
            'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.1 0.05]);

        uimenu('label','Auto &YLim','accelerator','y','callback','c_4_v_gui(''autoY'')');
        uimenu('label','&Distance','accelerator','d','callback','c_4_v_gui(''distance'')');
        uimenu('label','Click&Times','accelerator','t','callback','c_4_v_gui(''click_times'')');
        uimenu('label','New&Variable','accelerator','v','callback','c_4_v_gui(''new_var_enter'')');
        ud.columns=uimenu('label','&Columns','accelerator','c');
        uimenu(ud.columns,'label','1 (time)');
        for j_col=2:size(var1,2)
            %    hcol(j_col)=uimenu(ud.columns,'label',num2str(j_col));
            eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''c_4_v_gui(''''c' num2str(j_col) ''''')'');'];
            eval(eval_str);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'dt'
        hh=ud.h(1);  % use the first subplot to estimate available time interval
        udfig=get(gcf,'userdata');
        xl=get(hh,'XLim')+udfig.t_start_epoch;
        yl=get(hh,'YLim');
        hc=get(hh,'Children');
        dt=eval(['[' get(ud.dt,'string') ']']);
        tstr=['[' num2str(dt,'%7.2f') '] s'];
        t=0.5*(xl(1)+xl(2))+dt;
        if max(abs(dt))==0
            vstr='0 * [0 0 0]';
        else
            v=c_v(t,coord_sys(ud));
            vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%6.2f') ']'];
        end
        set(ud.v,'string',vstr);
        axes(ud.h(2));
        if eval(get(ud.filter,'string'))<1
            x1=var1;Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
            c_eval('x?=irf_tlim(var?,xl+[-20/Fs 20/Fs]);x?=irf_filt(x?,0,flim,Fs,5);');
            c_pl_tx('x?',var_col,dt);
        else c_pl_tx('var?',var_col,dt);
        end
        irf_zoom(xl,'x');irf_zoom(yl,'y');
        add_timeaxis;
        text(.5,-.6,['t_{2nd panedel} = t_{1st panel} - dt\newline  dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s ' coord_sys(ud) ],'units','normalized','verticalalignment','top');

    case 'v'
        hh = ud.h(1);  % use the first subplot to estimate available time interval
        udf = get(gcf,'userdata');
        xl = get(hh,'XLim'); yl = get(hh,'YLim');
        hc = get(hh,'Children');
        v = eval(['[' get(ud.v,'string') ']']);
        t = udf.t_start_epoch +0.5*(xl(1) +xl(2));
        if max(abs(v))==0
            dt=[0 0 0 0];
        else
            dt=c_v([t v],coord_sys(ud));
            ref_satellite_string=get(ud.ref_satellite,'string');
            ref_satellite=str2num(ref_satellite_string);
            if ref_satellite<1 || ref_satellite>4, ref_satellite=1;end
            dt=dt-dt(ref_satellite);
        end
        tstr=['[' num2str(dt,'%7.2f') ']'];
        if norm(v) > 0,
            vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%6.2f') ']'];
        else
            vstr='0*[0 0 0]';
        end
        set(ud.dt,'string',tstr);
        axes(ud.h(2));
        if eval(get(ud.filter,'string'))<1,
            x1=var1;Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
            for ic=1:4,eval(irf_ssub('x?=irf_tlim(var?,xl+[-20/Fs 20/Fs]);x?=irf_filt(x?,0,flim,Fs,5);',ic)),end
            c_pl_tx(x1,x2,x3,x4,var_col,dt);
        else
            c_pl_tx(var1,var2,var3,var4,var_col,dt);
        end
        axis([xl yl]);add_timeaxis;
        text(.5,-.6,['t_{2nd panedel} = t_{1st panel} - dt\newline dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s ' coord_sys(ud)],'units','normalized','verticalalignment','top');
        
    case 'click_times'
        axes(ud.h(1));zoom off;
        if isempty(ic), ic=0;ud.dtv=[];end
        if ic==0,
            set(gcf,'windowbuttondownfcn', 'c_4_v_gui(''click_times'')');
            ic=1;
        else
            p = get(gca, 'currentpoint');
            ud.dtv(ic)=p(1);
            ic=ic+1;
        end
        title(['click on s/c ' num2str(ic)]);
        if ic==5,
            set(gcf,'windowbuttondownfcn', '');
            title('');
            ic=0;
            dt=ud.dtv-ud.dtv(1);
            tstr=['[' num2str(dt,'%8.2f') ']'];
            set(ud.dt,'string',tstr);
            zoom on;
            c_4_v_gui('dt');
        end
    case 'distance'
        axes(ud.h(2));
        hh=ud.h(1);  % use the first subplot to estimate available time interval
        xl=get(hh,'XLim');
        v=eval(['[' get(ud.v,'string') ']']);
        tcenter = mean(xl);distance=norm(v)*diff(xl)/2;logd=log10(distance);
        if logd>round(logd), dx=10^(round(logd))/2;
        else dx=10^(round(logd))/5;
        end
        xticks=[-30:30]*dx/norm(v)/5+tcenter;
        xticklabels=cell(size(-30:30));
        for j=-30:30, xticklabels{j+31}=' ';end
        for j=-6:6, xticklabels{j*5+31}=num2str(j*dx);end
        set(ud.h(2),'xtick',xticks,'xticklabel',xticklabels);
        xlabel('km');
    case 'autoY'
        for h=ud.h(1:2),
            set(h,'YLimMode','auto');
        end
    case 'ylabel'
        h_select=get(ud.ylabpanel,'Value')-1;
        i_select=get(ud.ylab,'Value');
        strings=get(ud.ylab,'String');
        label=strings{i_select};
        if strcmp(label,'enter'); label=input('Input Y label>','s');end
        axes(ud.h(h_select));
        ylabel(label,'verticalalignment','bottom');
    case 'toggle'
        if ud.from
            set(gcbf, 'pointer', 'right');
            set(ud.helph, 'string', 'Click on axis selects ''To'' time');
            ext = get(ud.helph, 'extent');
            set(ud.helph, 'position', [395-ext(3) 55 ext(3:4)]);
            ud.from = 0;
        else
            set(gcbf, 'pointer', 'left');
            set(ud.helph, 'string', 'Click on axis selects ''From'' time');
            ext = get(ud.helph, 'extent');
            set(ud.helph, 'position', [5 55 ext(3:4)]);
            ud.from = 1;
        end
        set(gcbf, 'userdata', ud);
end
end

function label=var_label(var_str,var_col)
try
    dd=c_desc(irf_ssub(var_str,1));
    %      variable_label=[dd.labels{var_col} '[' dd.units{var_col} ']'];
    label=[var_str '[' num2str(var_col) ']'];
catch
    label=var_str;
end
end

function coord_sys_label=coord_sys(ud)
coord_sys_flag=get(ud.coord_sys,'value');
        if coord_sys_flag == 1
            coord_sys_label='GSM';
        else 
            coord_sys_label='GSE';
        end
end
