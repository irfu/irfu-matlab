function varargout = c_cal_gui(varargin)
% EFW calibration GUI
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%
persistent h0;
%disp(h0)

sp = '.';

inactive_color = [0.831373 0.815686 0.784314]; %gray
active_color = [1 1 1]; %white
pxa = .07+.12; pya = .1; wa = .47; ha = .215; dya = .01; % positions

if nargin, action = varargin{1};
else, action = 'init';
end

if regexp(action,'^update_DATA.+radiobutton$')
	vs = action(12:end-11);
	action = 'update_DATAradiobutton';
elseif regexp(action,'^update_DATA.+checkbox$')
	vs = action(12:end-8);
	action = 'update_DATAcheckbox';
elseif regexp(action,'^update_C[1-4]checkbox$')
	cl_id = eval(action(9));
	action = 'update_Ccheckbox';
end

switch action
case 'init'
	% Create figure
	h0 = figure(23);
	clf
	set(h0,'Position', [25 40 990 640])
	hnd = guihandles(h0);
	hnd.DATApanel = uipanel('Position',[.85 .01 .14 .98]);
	set(hnd.DATApanel,'Title','Data')
	
	% Load data
	
	c_Edata = {'diE?p1234','diEs?p12','diEs?p34'};
	c_Ddata = {'diEDI?'};
	c_Vdata = {'diVCp?','diVCh?'};
	c_AUXdata = {'P?','NCp?','NCh?'};
	dd={c_Edata c_Ddata c_Vdata c_AUXdata};
	%d_s = {'Edata', 'Ddata', 'Vdata', 'AUXdata'};
	%for d=1:4, eval(['hnd.' d_s{d} '= {};']),end
	hnd.Data = {};
	hnd.DataList = {};
	hnd.AUXList = {};
	hnd.BPPData = {};
	hnd.BData = {};
	hnd.EFWoffset = {};
	hnd.CISHoffset = {};
	hnd.CISCoffset = {};
	hnd.mode = 1; % 0 is for V, 1 is for E.
	% hnd.last
	% empty means replotting all, 
	% otherwise contains a cell array of data in a plot queue 
	hnd.last = [];
	hnd.ang_limit = 15; % 15 degrees.
	hnd.tlim = [0 0];
	hnd.c_visible = [1 1 1 1]; % all sc are visible by default
	no_active = 1;
	ncdata = 0; % number of non-AUX variables
	
	for cl_id=1:4
		s = av_ssub('C?',cl_id);
		pxd = .1; pyd = .95-(cl_id-1)*.24;
		hhh = uicontrol(hnd.DATApanel,'Style','checkbox',...
			'Units','normalized','Position',[pxd pyd 0.4 .05],...
			'String',s,'Value',hnd.c_visible(cl_id),...
			'Callback',['c_cal_gui(''update_' s 'checkbox'')'],...
			'Tag',[s 'checkbox']);
		ndata = 0;
		
		% Load FGM data
		clear BFGM BPP
		BFGM = []; BPP = [];
		if c_load(av_ssub('diB?',cl_id))
			c_log('load',av_ssub('loaded  diB?',cl_id))
			c_eval('BFGM=diB?;clear diB?',cl_id)
		else
			if c_load(av_ssub('diBPP?',cl_id))
				c_log('load',av_ssub('loaded  diB?',cl_id))
				c_eval('BFGM=diBPP?;clear diBPP?',cl_id)
			end
		end
		hnd.BData = [hnd.BData {BFGM}];
		hnd.BPPData = [hnd.BPPData {BPP}];
		clear BFGM BPP
		
		% Load offsets
		old_pwd = pwd; cd(sp);
		clear Delta_off;
		c_eval('load mEDSI D?p12p34',cl_id)
		if exist(av_ssub('D?p12p34',cl_id),'var')
			c_eval('Delta_off=D?p12p34;clear D?p12p34',cl_id)
		else
			c_log('load',...
			'Cannot load D#p12p34. Probably we have only one probe pair')
			Delta_off = 0;
		end
		cd(old_pwd)
		
		% Load data
		for d=1:4
			for j=1:length(dd{d})
				vs = av_ssub(dd{d}{j},cl_id);
				if c_load(vs) 
					ndata = ndata + 1;
					c_log('load',['loaded ' vs])
					if eval(['any(any(isnan(' vs '(:,2:end))))']),
						c_log('load',[vs ' contains no data'])
						continue
					end
					dsc = c_desc(vs);
					clear data
					data.name = vs;
					data.cl_id = str2num(dsc.cl_id);
					data.inst = dsc.inst;
					data.label = [dsc.inst ' ' dsc.sig];
					if dsc.sen, data.label = [data.label ' (' dsc.sen ')'];, end
					data.sen = dsc.sen;
					[data.plot_style,data.plot_color] =...
						p_style(data.cl_id, data.inst,data.sen);
					
					% tlim
					c_eval(['tlxxx(1) = ' vs '(1,1); tlxxx(2) = ' vs '(end,1);'],cl_id)
					if ~hnd.tlim(1), hnd.tlim(1) = tlxxx(1);
					elseif hnd.tlim(1) > tlxxx(1), hnd.tlim(1) = tlxxx(1); 
					end
					if hnd.tlim(2) < tlxxx(2), hnd.tlim(2) = tlxxx(2); end
					clear tlxxx
					
					% Correct delta offset for Es?p12 and Es?p34
					if Delta_off
						if isreal(Delta_off) & strcmp(vs,av_ssub('diEs?p12',cl_id))
							c_log('calb','correcting delta offset in p12')
							c_eval(...
							'diEs?p12=diEs?p12-ones(length(diEs?p12),1)*[0 Delta_off 0];',...
							cl_id)
						elseif ~isreal(Delta_off) & strcmp(vs,av_ssub('diEs?p34',cl_id))
							c_log('calb','correcting delta offset in p34')
							c_eval(...
							'diEs?p34=diEs?p34-ones(length(diEs?p34),1)*[0 imag(Delta_off) 0];',...
							cl_id)
						end
					end
					eval(['data.data =' vs ';'])
					data.p_data = [];
					if ((strcmp(dsc.sen,'all') | strcmp(dsc.sen,'p12'))& ...
						strcmp(dsc.sig,'E')) | strcmp(dsc.sen,'COD')
						data.visible = 0;
					else
						data.visible = 1;
					end
					if (strcmp(dsc.inst,'EFW') & strcmp(dsc.sen,'E'))|...
						(strcmp(dsc.inst,'CIS') & (strcmp(dsc.sen,'V')|...
						strcmp(dsc.sen,'Vp'))), data.editable = 1;
					else, data.editable = 0;
					end
					if strcmp(dsc.sig,'E'), data.type = 'E';
					elseif strcmp(dsc.sig,'V') | strcmp(dsc.sig,'Vp')
						data.type = 'V';
					else, data.type = 'A';
					end
					if d==1|d==3
						hhd = uicontrol(hnd.DATApanel,'Style','radiobutton',...
							'Units','normalized',...
							'Position',[pxd pyd-.02*ndata 0.1 .02],...
							'String','',...
							'Callback',['c_cal_gui(''update_DATA' vs 'radiobutton'')'],...
							'Tag',['DATA' vs 'radiobutton']);
						if strcmp(dsc.sig,'V')|strcmp(dsc.sig,'Vp')
							set(hhd,'Enable','off')
						end
						eval(['hnd.DATA' vs 'radiobutton=hhd;clear hhd'])
					end
					hhd = uicontrol(hnd.DATApanel,'Style','checkbox',...
						'Units','normalized',...
						'Position',[pxd+.1 pyd-.02*ndata 0.8 .02],...
						'String',data.label,...
						'Callback',['c_cal_gui(''update_DATA' vs 'checkbox'')'],...
						'Tag',['DATA' vs 'checkbox']);
					if data.visible, set(hhd,'Value',1), end
					if no_active & data.visible & ...
						strcmp(dsc.inst,'EFW') & strcmp(dsc.sig,'E')
						eval(['set(hnd.DATA' vs 'radiobutton,''Value'',1)'])
						no_active = 0;
						hnd.ActiveVar = vs;
					end
					eval(['hnd.DATA' vs 'checkbox=hhd;clear hhd'])
					if d < 4
						% E and V data
						% Resample B
						if ~isempty(hnd.BData{cl_id})
							c_log('proc','resampling B GFM')
							c_eval(...
								['data.B = c_resamp(hnd.BData{cl_id},' vs ');'],...
								cl_id)
						elseif ~isempty(hnd.BPPData{cl_id})
							c_log('proc','resampling B PP')
							c_eval(...
								['data.B = c_resamp(hnd.BPPData{cl_id},' vs ');'],...
								cl_id)
						else, data.B = [];
						end
						data.aux = 0;
						ncdata = ncdata + 1;
					else
						% AUX data
						data.aux = 1;
					end
					hnd.Data = [hnd.Data, {data}];
					eval(['clear ' vs])
					
				else
					c_log('load',['cannot load ' vs])
				end
			end
		end
		if ndata==0
			set(hhh,'Value',0,'Enable','off')
			hnd.c_visible(cl_id) = 0;
		end
		eval(['hnd.' s 'checkbox=hhh;clear hhh'])
		
		% Load EFW offsets
		old_pwd = pwd; cd(sp);
		offset = [1+0i 1];
		c_eval('load mEDSI Ddsi? Damp?',cl_id)
		if exist(av_ssub('Ddsi?',cl_id),'var')
			c_eval('offset(1)=Ddsi?;clear Ddsi?',cl_id)
			c_log('load',...
			sprintf('loading DSI offset from file Ddsi=%.2f %.2f*i',...
			real(offset(1)),imag(offset(1))))
		end
		if exist(av_ssub('Damp?',cl_id),'var')
			c_eval('offset(2)=Damp?;clear Damp?',cl_id)
			c_log('load',...
			sprintf('loading amplitude correction from file Damp=%.2f',offset(2)))
		end
		cd(old_pwd)
		hnd.EFWoffset = [hnd.EFWoffset {offset}];
		
		% Load CIS offsets, must be only in Z.
		offset = [0 0 0];
		hnd.CISHoffset = [hnd.CISHoffset {offset}];
		hnd.CISHoffset = [hnd.CISCoffset {offset}];
		
		hnd.cal_updated = 1;
	end
	
	% Check if we have any data apart from AUX
	if ~ncdata, error('No usefull data loaded'), end
	
	% Create Data Axes
	hnd.Xaxes = axes('Position',[pxa pya+(ha+dya)*3 wa ha],'Tag','Xaxes');
	hnd.Yaxes = axes('Position',[pxa pya+(ha+dya)*2 wa ha],'Tag','Yaxes');
	hnd.Zaxes = axes('Position',[pxa pya+(ha+dya)*1 wa ha],'Tag','Zaxes');
	hnd.AUXaxes = axes('Position',[pxa pya wa ha],'Tag','AUXaxes');
	
	% Create Legend Axes
	hnd.DLaxes = axes('Position',[.01 pya+(ha+dya) .12 ha+(ha+dya)*2 ],...
		'Visible','off','XTick',[],'YTick',[]);
	hnd.ALaxes = axes('Position',[.01 pya .12 ha ],...
		'Visible','off','XTick',[],'YTick',[]);
		
	% Create Calibrators
	wp = .15;
	hnd.DXpanel = uipanel('Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*3 wp ha]);
	hnd.DYpanel = uipanel('Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*2 wp ha]);
	hnd.DZpanel = uipanel('Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*1 wp ha]);
	set(hnd.DXpanel,'Title','dX')
	set(hnd.DYpanel,'Title','dY')
	set(hnd.DZpanel,'Title','dZ')
	hnd.DXslider = uicontrol(hnd.DXpanel,'Style','slider',...
		'Units','normalized','Position',[0.1 0.7 0.8 0.2],...
		'Max',5,'Min',-5,'Value',0,'Callback','c_cal_gui(''update_DXslider'')');
	hnd.DXedit = uicontrol(hnd.DXpanel,'Style','edit',...
		'Units','normalized','Position',[0.1 0.1 0.4 0.2],...
		'String','0.0','BackgroundColor',active_color,...
		'Callback','c_cal_gui(''update_DXedit'')','Tag','DXedit');
	hnd.DXcheckbox = uicontrol(hnd.DXpanel,'Style','checkbox',...
		'Units','normalized','Position',[0.6 0.1 0.4 0.2],...
		'String','Lock',...
		'Callback','c_cal_gui(''update_DXcheckbox'')','Tag','DXcheckbox');

	guidata(h0,hnd);
	c_cal_gui('replot_all')
	av_figmenu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fix_plot_pos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'fix_plot_pos'
	hnd = guidata(h0);
	h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
	for ax=1:4
		set(h(ax),'Position',[pxa pya+(ha+dya)*(4-ax) wa ha]);
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'replot'
	% Plot data
	hnd = guidata(h0);
	h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
	if ~length(hnd.Data), error('no data to calibrate'),end
	if isempty(hnd.last), disp('List empty, will do nothing'), return, end

	% Update selected records in DataLegList
	% we set hnd.last to name of the variable we need to add
	d_ii = D_findByNameList(hnd.Data,hnd.last);
	
	% Sanity check
	if isempty(d_ii)
		error(['cannot find variables in the data list'])
	end
	
	% Plotting 
	for j=1:length(d_ii)
		%disp(['replot: plotting ' hnd.Data{d_ii(j)}.name])
		
		% Process only visible data
		if ~hnd.Data{d_ii(j)}.visible, continue, end
		
		if hnd.Data{d_ii(j)}.aux
			hnd.AUXList = L_add(hnd.AUXList,hnd.Data{d_ii(j)}.name);
			
			% Plotting
			hold(h(4),'on')
			hnd.Data{d_ii(j)}.ploth = plot(hnd.AUXaxes,...
				hnd.Data{d_ii(j)}.data(:,1),hnd.Data{d_ii(j)}.data(:,2),...
				hnd.Data{d_ii(j)}.plot_style,...
				'Color',hnd.Data{d_ii(j)}.plot_color);
			hold(h(4),'off')
		else
			% We update p_data only if calibrations were changed
			if isempty(hnd.Data{d_ii(j)}.p_data) | hnd.cal_updated
				%disp('replot: recalculating plot data')
				hnd.Data{d_ii(j)}.p_data =...
					get_plot_data(hnd.Data{d_ii(j)}, hnd);
			end
			if ~isempty(hnd.Data{d_ii(j)}.p_data)
				hnd.DataList = L_add(hnd.DataList,hnd.Data{d_ii(j)}.name);
				
				% Plotting
				for ax=1:3
					hold(h(ax),'on')
					hnd.Data{d_ii(j)}.ploth(ax) = plot(h(ax),...
						hnd.Data{d_ii(j)}.p_data(:,1),...
						hnd.Data{d_ii(j)}.p_data(:,1+ax),...
						hnd.Data{d_ii(j)}.plot_style,...
						'Color',hnd.Data{d_ii(j)}.plot_color);
					hold(h(ax),'off')
				end
			end
		end
	end
	hnd.last = [];
	hnd.cal_updated = 0;

	av_zoom(hnd.tlim,'x',h);
	guidata(h0,hnd);
	c_cal_gui('update_legend')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replot_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'replot_all'
	% Replot all data
	hnd = guidata(h0);
	h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
	if ~length(hnd.Data), error('no data to calibrate'),end
	
	if isempty(hnd.last)
		% Clear everything
		for ax=1:3, cla(h(ax)), end
		hnd.DataList = {};
		hnd.AUXList = {};
		
		% Plot
		for j=1:length(hnd.Data), hnd.last = L_add(hnd.last,hnd.Data{j}.name); end
		guidata(h0,hnd);
		c_cal_gui('replot')
		hnd = guidata(h0);
		
		% Axes labels
		labs = ['x' 'y' 'z'];
		if hnd.mode, u_s = 'E'; u_u = 'mV/m';
		else, u_s = 'V'; u_u = 'km/s';
		end
		for ax=1:3
			ylabel(h(ax),[u_s '_' labs(ax) ' [' u_u ']'])
			grid(h(ax),'on')
		end
		ylabel(h(4),'AUX')
		axes(h(4)); add_timeaxis; grid on
		
		% Time span
		av_zoom(hnd.tlim,'x',h);
		
		guidata(h0,hnd);
	else
		% Replotting NOT ALL data (hnd.last is set)
		c_cal_gui('replot')
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'update_legend'
	hnd = guidata(h0);
	if isempty(hnd.DataList)
		set(hnd.DLaxes,'Visible','off')
		legend(hnd.DLaxes,'off')
	else
		cla(hnd.DLaxes)
		set(hnd.DLaxes,'Visible','on','XTick',[],'YTick',[],'Box','on')
		
		% Make fake plots
		for j=1:length(hnd.DataList)
			ii = D_findByName(hnd.Data,hnd.DataList{j});
			hold(hnd.DLaxes,'on')
			plot(hnd.DLaxes,1,1,...
				hnd.Data{ii}.plot_style,...
				'Color',hnd.Data{ii}.plot_color,...
				'Visible','off')
			hold(hnd.DLaxes,'off')
		end
		
		% Make legend
		ii = D_findByNameList(hnd.Data,hnd.DataList);
		l_s = ['''' hnd.Data{ii(1)}.label ''''];
		if length(ii)>1
			for j=2:length(ii)
				l_s = [l_s ',''' hnd.Data{ii(j)}.label ''''];
			end
		end
		eval(['hxxx=legend(hnd.DLaxes,' l_s ');'])
		set(hxxx,'FontSize',7);
		legend(hnd.DLaxes,'boxoff')
		set(hnd.DLaxes,'Visible','off')
	end
	if isempty(hnd.AUXList)
		set(hnd.ALaxes,'Visible','off')
		legend(hnd.ALaxes,'off')
	else
		cla(hnd.ALaxes)
		set(hnd.ALaxes,'Visible','on','XTick',[],'YTick',[],'Box','on')
		
		% Make fake plots
		for j=1:length(hnd.AUXList)
			ii = D_findByName(hnd.Data,hnd.AUXList{j});
			hold(hnd.ALaxes,'on')
			plot(hnd.ALaxes,1,1,...
				hnd.Data{ii}.plot_style,...
				'Color',hnd.Data{ii}.plot_color,...
				'Visible','off')
			hold(hnd.ALaxes,'off')
		end
		
		% Make legend
		ii = D_findByNameList(hnd.Data,hnd.AUXList);
		l_s = ['''' hnd.Data{ii(1)}.label ''''];
		if length(ii)>1
			for j=2:length(ii)
				l_s = [l_s ',''' hnd.Data{ii(j)}.label ''''];
			end
		end
		eval(['hxxx=legend(hnd.ALaxes,' l_s ');'])
		set(hxxx,'FontSize',7);
		legend(hnd.ALaxes,'boxoff')
		set(hnd.ALaxes,'Visible','off')
	end
	guidata(h0,hnd);
	c_cal_gui('fix_plot_pos')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_DXslider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'update_DXslider'
	hnd = guidata(h0);
	set(hnd.DXedit,'String',...
    num2str(get(hnd.DXslider,'Value')));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_DXedit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'update_DXedit'
	hnd = guidata(h0);
	val = str2double(get(hnd.DXedit,'String'));
	% Determine whether val is a number between 0 and 1
	if isnumeric(val) & length(val)==1 & ...
		val >= get(hnd.DXslider,'Min') & ...
		val <= get(hnd.DXslider,'Max')
		set(hnd.DXslider,'Value',val);
	else
		% Increment the error count, and display it
		set(hnd.DXedit,'String','You have entered an invalid entry ');
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_DXcheckbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'update_DXcheckbox'
	hnd = guidata(h0);
	if get(hnd.DXcheckbox,'Value')==1
		% Lock
		set(hnd.DXedit,'Enable','off','BackgroundColor',inactive_color)
		set(hnd.DXslider,'Enable','off')
	else
		% Unclock
		set(hnd.DXedit,'Enable','on','BackgroundColor',active_color)
		set(hnd.DXslider,'Enable','on')
		set(hnd.DXslider,'Value',str2double(get(hnd.DXedit,'String')))
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_Ccheckbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'update_Ccheckbox'
	hnd = guidata(h0);
	ii = D_findByCLID(hnd.Data,cl_id);
	if isempty(ii), return, end
	
	if get(eval(['hnd.C' num2str(cl_id) 'checkbox']),'Value')==1
		% Show
		for j=1:length(ii)
			hnd.Data{ii(j)}.visible = 1;
			set(eval(['hnd.DATA' hnd.Data{ii(j)}.name 'checkbox']),...
				'Value', 1);
		end
		hnd.last = D_listNames(hnd.Data,ii);
		guidata(h0,hnd);
		c_cal_gui('replot');
	else
		% Hide
		kk = [];
		for j=1:length(ii)
			if hnd.Data{ii(j)}.visible == 0, kk = [kk j]; end
		end
		ii(kk) =[];
		if isempty(ii), return, end
		
		for j=1:length(ii)
			%disp(['hide ' hnd.Data{ii(j)}.name])
			delete(hnd.Data{ii(j)}.ploth)
			hnd.Data{ii(j)}.visible = 0;
			set(eval(['hnd.DATA' hnd.Data{ii(j)}.name 'checkbox']),...
				'Value', 0);
		end
		
		% remove from the plotting lists
		ii = cell(1,2);
		ii{1} = sort(D_findByNameList(hnd.Data,hnd.DataList));
		ii{2} = sort(D_findByNameList(hnd.Data,hnd.AUXList));
		for k=1:2
			for j=length(ii{k}):-1:1
				if hnd.Data{ii{k}(j)}.cl_id == cl_id
					if k==1, hnd.DataList(j) = [];
					else, hnd.AUXList(j) = [];
					end
				end
			end
		end
		
		guidata(h0,hnd);
		c_cal_gui('update_legend')
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update_DATAcheckbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'update_DATAcheckbox'
	hnd = guidata(h0);
	if get(eval(['hnd.DATA' vs 'checkbox']),'Value')==1
		% Plot the varible
		j = D_findByName(hnd.Data,vs);
		%disp(['plotting ' hnd.Data{j}.name])
		hnd.Data{j}.visible = 1;
		hnd.last = {vs};
		guidata(h0,hnd);
		c_cal_gui('replot')
		
		% Check if we need to show C# checkBox
		if get(eval(['hnd.C' num2str(hnd.Data{j}.cl_id) 'checkbox']),'Value')==0
			set(eval(['hnd.C' num2str(hnd.Data{j}.cl_id) 'checkbox']),...
				'Value',1)
		end
	else
		% Hide the varible
		j = D_findByName(hnd.Data,vs);
		hnd.Data{j}.visible = 0;
		delete(hnd.Data{j}.ploth)

		% Check if we need to hide C# checkBox
		cl_id = hnd.Data{j}.cl_id;
		ii = D_findByCLID(hnd.Data,cl_id);
		if isempty(ii), return, end
		hide_ok = 1;
		for j=1:length(ii)
			if hnd.Data{ii(j)}.visible == 1, hide_ok = 0; break, end
		end
		if hide_ok
			set(eval(['hnd.C' num2str(cl_id) 'checkbox']),...
				'Value',0)
		end
		
		% Remove from the plotting lists	
		j = L_find(hnd.DataList,vs);
		if ~isempty(j), hnd.DataList = L_rm_ii(hnd.DataList,j);
		else, hnd.AUXList = L_rm(hnd.AUXList,vs);
		end
		
		guidata(h0,hnd);
		c_cal_gui('update_legend')
	end
otherwise 
	disp('wrong action')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function corr_v_velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = corr_v_velocity(v,offset)
out = data;
out(:,2) = v(:,2) - offset(1);
out(:,3) = v(:,3) - offset(2);
out(:,4) = v(:,4) - offset(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function p_style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out, color] = p_style(cl_id, inst, sen)
out = '';
colrs = [0 0 0; 1 0 0; 0 .5 0; 0 0 1]; % black, red, dark green, blue
switch inst
case 'EFW'
	switch sen
	case 'all'
		out = '-';
	case 'p12'
		out = '--';
	case 'p34'
		out = '-.';	
	otherwise
		disp('unknown sensor')
	end
case 'EDI'
	out = '.';
case 'CIS'
	switch sen
	case 'HIA'
		out = '-o';
	case 'COD'
		out = '-+';
	otherwise
		disp('unknown sensor')
	end
otherwise
	disp('unknown instrument')
end
color = colrs(cl_id,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function get_plot_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p_data = get_plot_data(data, hnd)
% correct offsets, compute Ez, compute ExB (VxB).

p_data = [];

if data.visible
	p_data = data.data;
	
	% Correct offsets
	if data.editable 
		switch data.type
		case 'E'
			ofs = hnd.EFWoffset{data.cl_id};
			p_data = corrDSIOffsets(p_data,...
				real(ofs(1)),imag(ofs(1)),ofs(2));
		case 'V'
			if strcmp(data.sen,'COD')
				offset = hnd.CISCoffset{data.cl_id};
			elseif strcmp(data.sen,'HIA')
				offset = hnd.CISHoffset{data.cl_id};
			else, offset = [0 0 0];
			end
			p_data = corr_v_velocity(p_data,...
				hnd.EFWoffset{data.cl_id});
		otherwise
			c_log('proc','Unknown data type.')
		end
	end
	
	% Calculate Ez, convert V->E and V->E
	switch data.type
	case 'E'
		if any(p_data(:,4))==0,
			if isempty(data.B)
				c_log('proc','B is empty. No Ez')
			else
				[p_data,angle] = av_ed(p_data,...
					data.B,hnd.ang_limit);
				ii = find(abs(angle) < hnd.ang_limit);
				if length(ii) > 1
					p_data(ii,4) = p_data(ii,4)*NaN; 
				end
			end
		end
		% E->V if in V mode
		if ~hnd.mode
			if ~isempty(data.B)
				p_data = av_e_vxb(p_data,hnd.Data{j}.B,-1);
			else
				data.visible = 0; 
				p_data = [];
				c_log('proc','B is empty. No ExB')
			end
		end
	case 'V'
		% V->E if in E mode
		if hnd.mode
			if ~isempty(data.B)
				p_data = av_t_appl(...
					av_cross(p_data,data.B),'*(-1e-3)');
			else
				data.visible = 0;
				p_data = [];
				c_log('proc','B is empty. No VxB')
			end
		end
	otherwise
		% AUX data, do nothing
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function D_findByName
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_i = D_findByName(data_list,name_s)
%data_list is a cell array
%name_s is a string
data_i = [];
for j=1:length(data_list)
	if strcmp(data_list{j}.name,name_s)
		data_i = j;
		return
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function D_findByNameList
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_ii = D_findByNameList(data_list,name_list)
%data_list is a cell array
%name_list is a cell array of strings
data_ii = [];
for j=1:length(name_list)
	data_rec = D_findByName(data_list,name_list{j});
	if ~isempty(data_rec), data_ii = [data_ii data_rec]; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function D_findByCLID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_ii = D_findByCLID(data_list,cl_id)
%data_list is a cell array
%cl_id is integer [1-4]
data_ii = [];
for j=1:length(data_list)
	if data_list{j}.cl_id == cl_id
		data_ii = [data_ii j];
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function D_listNames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list = D_listNames(data_list,ii)
%data_list is a cell array
%ii is index array
list = cell(size(ii));
for j=1:length(ii)
	list{j} = data_list{ii(j)}.name;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L_add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newlist = L_add(list,s)
newlist = [list {s}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L_find
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ii = L_find(list,s_list)
ii = [];
if isstr(s_list)
	% fast search
	for j=1:length(list)
		if strcmp(list{j},s_list), ii = j; return, end
	end
else
	for k=1:length(s_list)
		for j=1:length(list)
			if strcmp(list{j},s_list{k}), ii = [ii j]; break, end
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L_rm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newlist = L_rm(list,s_list)
newlist = list;
if isstr(s_list)
	j = L_find(list,s_list);
	if ~isempty(j)
		newlist(j) = [];
	end
else
	ii = L_find(list,s_list);
	ii = sort(ii);
	if ~isempty(ii)
		for j=length(ii):-1:1, newlist(j) = []; end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L_rm_ii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newlist = L_rm_ii(list,ii)
newlist = list;
if length(ii)==1
	newlist(ii) = [];
else
	ii = sort(ii);
	for j=length(ii):-1:1, newlist(j) = []; end
end
