function c_4_v_gui(x1,x2,x3,x4,column)
%C_4_V_GUI interactive discontinuity analyzer for Cluster
%
%  status = c_4_v_gui(x1,x2,x3,x4,column)
%  status = c_4_v_gui('B?',column)   use variables B1,B2,B3,B4 from base workspace
%  status = c_4_v_gui('B?')          use column 2
%
% $Id$

flag_first_call=0;

if       (nargin<=2 && ischar(x1)), % either action as parameter or string variable
	if strfind(x1,'?'),
		figure;ud=[]; % intialize
		ud.variable_str=x1;
		if nargin == 1, ud.var_col=2;else ud.var_col=x2;end;
		action='new_var';
		flag_first_call=1;
	else
		action=x1;
		ud=get(gcf,'userdata');
	end
elseif   (nargin ==4) || (nargin == 5),
	if nargin ==4, irf_log('fcal','Using second column');column=2;end
	figure;ud=[]; % intialize
	c_eval('ud.var?=x?;');
	ud.var_col=column;
	action='initialize';
else      help c_4_v_gui
end

irf_log('fcal',['action=' action]);
switch action,
	case {'c1','c2','c3','c4','c5','c6'}
		ud.var_col=str2num(action(2:end));
		set(gcf,'userdata',ud);
		c_4_v_gui('update_var_col');
	case 'update_var_col'
		hca=ud.h(1);
		xl=get(ud.h(1),'XLim');
		c_pl_tx(hca,ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col);zoom(hca,'on');
		ylabel(hca,var_label(ud.variable_str,ud.var_col));
		axis(hca,[xl(1) xl(2) 0 1]);
		irf_zoom(hca,'y');irf_timeaxis(hca)
		c_pl_tx(ud.h(2),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col,ud.dt);
%		c_4_v_gui('dt');
	case 'new_var_enter'
		xx=inputdlg('Enter new variable mask. Examples: B? or R? or P?p1','**',1,{'B?'});
		ud.variable_str=xx{1};
		set(gcf,'userdata',ud);
		c_4_v_gui('new_var');
	case 'new_var'
		evalin('base',['if ~exist(''' irf_ssub(ud.variable_str,1) '''), c_load(''' ud.variable_str ''');end' ]);
		c_eval('ud.var?=evalin(''base'',irf_ssub(ud.variable_str,?));');
		if ud.var_col > size(ud.var1,2), ud.var_col=2;end % in case new variable has less columns
		if flag_first_call,
			set(gcf,'userdata',ud);
			c_4_v_gui('initialize');
		else
			if ishandle(ud.h(1)),
				for j_col=2:size(ud.var1,2)
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
					for j_col=(size(ud.var1,2)+1):length(ud.hcol)
						set(ud.hcol(j_col),'enable','off')
					end
				end
				set(gcf,'userdata',ud);
				c_4_v_gui('update_var_col');
			else
				set(gcf,'userdata',ud);
				c_4_v_gui('initialize');
			end
		end
	case 'initialize'
		irf_figmenu;
		set(gcf,'color','white'); % white background for figures (default is grey)
		set(gcf,'userdata',ud); % because c_pl_tx can also changed userdata)

		h(1)=subplot(3,1,1);
		c_pl_tx(h(1),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col);zoom(h(1),'on');
		ylabel(h(1),var_label(ud.variable_str,ud.var_col));
		
		h(2)=subplot(3,1,2);
		c_pl_tx(h(2),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col);
		ylabel(h(2),var_label(ud.variable_str,ud.var_col));
		
		ud=get(gcf,'userdata');
		
		irf_legend(0,['c\_4\_v\_int() ' datestr(now)],[0.01 0.99],'fontsize',7); % add information to the plot
		
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
		ud.dt_input = uicontrol('style', 'edit', ...
			'string', '[0 0 0 0]', ...
			'callback', 'c_4_v_gui(''dt'')', ...
			'backgroundcolor','white','units','normalized','position', [xp+0.15 yp 0.29 0.05]);
		ud.dt=[0 0 0 0]; % default values
		
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
		for j_col=2:size(ud.var1,2)
			%    hcol(j_col)=uimenu(ud.columns,'label',num2str(j_col));
			eval_str=['ud.hcol(j_col)=uimenu(ud.columns,''label'',''' num2str(j_col) ''',''callback'',''c_4_v_gui(''''c' num2str(j_col) ''''')'');'];
			eval(eval_str);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		set(gcf,'userdata',ud);
	case 'dt'
		hh=ud.h(1);  % use the first subplot to estimate available time interval
		xl=get(hh,'XLim')+ud.t_start_epoch;
		yl=get(hh,'YLim');
		hc=get(hh,'Children');
		ud.dt=eval(['[' get(ud.dt_input,'string') ']']);
		tstr=['[' num2str(ud.dt,'%9.2f') '] s'];
		t=0.5*(xl(1)+xl(2))+ud.dt;
		if max(abs(ud.dt))==0
			vstr='0 * [0 0 0]';
		else
			v=c_v(t,coord_sys(ud));
			vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%6.2f') ']'];
		end
		set(ud.v,'string',vstr);
		if eval(get(ud.filter,'string'))<1
			x1=ud.var1;Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
			c_eval('x?=irf_tlim(ud.var?,xl+[-20/Fs 20/Fs]);x?=irf_filt(x?,0,flim,Fs,5);');
			c_pl_tx(ud.h(2),'x?',ud.var_col,ud.dt);
		else
			c_pl_tx(ud.h(2),'ud.var?',ud.var_col,ud.dt);
		end
		irf_zoom(ud.h(2),'x',xl);
		irf_zoom(ud.h(2),'y',yl);
		irf_timeaxis(ud.h(2));
		text(.5,-.6,['t_{2nd panedel} = t_{1st panel} - dt\newline  dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s ' coord_sys(ud) ],'units','normalized','verticalalignment','top','paren',ud.h(2));
		set(gcf,'userdata',ud);
	case 'v'
		hh = ud.h(1);  % use the first subplot to estimate available time interval
		xl = get(hh,'XLim'); yl = get(hh,'YLim');
		hc = get(hh,'Children');
		v = eval(['[' get(ud.v,'string') ']']);
		t = ud.t_start_epoch +0.5*(xl(1) +xl(2));
		if max(abs(v))==0
			dt=[0 0 0 0];
		else
			dt=c_v([t v],coord_sys(ud));
			ref_satellite_string=get(ud.ref_satellite,'string');
			ref_satellite=str2num(ref_satellite_string);
			if ref_satellite<1 || ref_satellite>4, ref_satellite=1;end
			dt=dt-dt(ref_satellite);
		end
		tstr=['[' num2str(dt,'%9.2f') ']'];
		if norm(v) > 0,
			vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%6.2f') ']'];
		else
			vstr='0*[0 0 0]';
		end
		set(ud.dt_input,'string',tstr);
		if eval(get(ud.filter,'string'))<1,
			x1=ud.var1;Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
			for ic=1:4,eval(irf_ssub('x?=irf_tlim(var?,xl+[-20/Fs 20/Fs]);x?=irf_filt(x?,0,flim,Fs,5);',ic)),end
			c_pl_tx(ud.h(2),x1,x2,x3,x4,ud.var_col,dt);
		else
			c_pl_tx(ud.h(2),ud.var1,ud.var2,ud.var3,ud.var4,ud.var_col,dt);
		end
		axis(ud.h(2),[xl yl]);
		irf_timeaxis(ud.h(2));
		text(.5,-.6,['t_{2nd panedel} = t_{1st panel} - dt\newline dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s ' coord_sys(ud)],'units','normalized','verticalalignment','top','paren',ud.h(2));
		set(gcf,'userdata',ud);
	case 'click_times'
		zoom(ud.h(1),'off');
		if (~isfield(ud,'ic') || isempty(ud.ic)), ud.ic=0;ud.dtv=[];end
		if ud.ic==0,
			set(gcf,'windowbuttondownfcn', 'c_4_v_gui(''click_times'')');
			ud.ic=1;
		else
			p = get(ud.h(1), 'currentpoint');
			ud.dtv(ud.ic)=p(1);
			ud.ic=ud.ic+1;
		end
		title(['click on s/c ' num2str(ud.ic)]);
		if ud.ic==5,
			set(gcf,'windowbuttondownfcn', '');
			title('');
			ud.ic=0;
			ud.dt=ud.dtv-ud.dtv(1);
			tstr=['[' num2str(ud.dt,'%10.2f') ']'];
			set(ud.dt_input,'string',tstr);
			zoom(ud.h(1),'on');
			set(gcf,'userdata',ud);
			c_4_v_gui('dt');
		end
		set(gcf,'userdata',ud);
	case 'distance'
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
		ylabel(ud.h(h_select),label,'verticalalignment','bottom');
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
		set(gcf, 'userdata', ud);
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
