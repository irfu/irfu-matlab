function varargout = c_cal_gui(varargin)
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%
persistent h0;
%disp(h0)

sp = '.';

inactive_color = [0.831373 0.815686 0.784314]; %gray
active_color = [1 1 1]; %white

if nargin, action = varargin{1};
else, action = 'init';
end

if regexp(action,'^update_DATA.+radiobutton$')
	vs = action(12:end-11);
	action = 'update_DATAradiobutton';
elseif regexp(action,'^update_DATA.+checkbox$')
	vs = action(12:end-8);
	action = 'update_DATAcheckbox';
end

switch action
case 'init'
	%create figure
	h0 = figure(23);
	clf
	set(h0,'Position', [25 40 990 640])
	handles = guihandles(h0);
	pxa = .07; pya = .1; wa = .55; ha = .2; dya = .03;
	handles.DATApanel = uipanel('Position',[.85 .01 .14 .98]);
	set(handles.DATApanel,'Title','Data')
	
	%load data
	
	c_Edata = {'diE?p1234','diEs?p12','diEs?p34'};
	c_Ddata = {'diEDI?'};
	c_Vdata = {'diVCp?','diVCh?'};
	c_AUXdata = {'P?','NCp?','NCh?'};
	dd={c_Edata c_Ddata c_Vdata c_AUXdata};
	%d_s = {'Edata', 'Ddata', 'Vdata', 'AUXdata'};
	%for d=1:4, eval(['handles.' d_s{d} '= {};']),end
	handles.Data = {};
	handles.AUXData = {};
	handles.DataList = {};
	handles.DataLegList = {};
	handles.AUXList = {};
	handles.AUXLegList = {};
	handles.BPPData = {};
	handles.BData = {};
	handles.EFWoffset = {};
	handles.CISHoffset = {};
	handles.CISCoffset = {};
	handles.mode = 1; % 0 is for V, 1 is for E.
	handles.ang_limit = 15; % 15 degrees.
	handles.tlim = [0 0];
	no_active = 1;
	
	for cl_id=1:4
		s = av_ssub('C?',cl_id);
		pxd = .1; pyd = .95-(cl_id-1)*.24;
		hhh = uicontrol(handles.DATApanel,'Style','checkbox',...
			'Units','normalized','Position',[pxd pyd 0.4 .05],...
			'String',s,...
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
		handles.BData = [handles.BData {BFGM}];
		handles.BPPData = [handles.BPPData {BPP}];
		
		% load offsets
		old_pwd = pwd; cd(sp);
		clear Delta_off;
		c_eval('load mEDSI D?p12p34',cl_id)
		if exist(av_ssub('D?p12p34',cl_id),'var')
			c_eval('Delta_off=D?p12p34;',cl_id)
		else
			c_log('load',...
			'Cannot load D#p12p34. Probably we have only one probe pair')
			Delta_off = 0;
		end
		cd(old_pwd)
		
		% load data
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
					
					%tlim
					c_eval(['tlxxx(1) = ' vs '(1,1); tlxxx(2) = ' vs '(end,1);'],cl_id)
					if ~handles.tlim(1), handles.tlim(1) = tlxxx(1);
					elseif handles.tlim(1) > tlxxx(1), handles.tlim(1) = tlxxx(1); 
					end
					if handles.tlim(2) < tlxxx(2), handles.tlim(2) = tlxxx(2); end
					clear tlxxx
					
					% correct delta offset for Es?p12 and Es?p34
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
						hhd = uicontrol(handles.DATApanel,'Style','radiobutton',...
							'Units','normalized',...
							'Position',[pxd pyd-.02*ndata 0.1 .02],...
							'String','',...
							'Callback',['c_cal_gui(''update_DATA' vs 'radiobutton'')'],...
							'Tag',['DATA' vs 'radiobutton']);
						if strcmp(dsc.sig,'V')|strcmp(dsc.sig,'Vp')
							set(hhd,'Enable','off')
						end
						eval(['handles.DATA' vs 'radiobutton=hhd;clear hhd'])
					end
					hhd = uicontrol(handles.DATApanel,'Style','checkbox',...
						'Units','normalized',...
						'Position',[pxd+.1 pyd-.02*ndata 0.8 .02],...
						'String',data.label,...
						'Callback',['c_cal_gui(''update_DATA' vs 'checkbox'')'],...
						'Tag',['DATA' vs 'checkbox']);
					if data.visible, set(hhd,'Value',1), end
					if no_active & data.visible & ...
						strcmp(dsc.inst,'EFW') & strcmp(dsc.sig,'E')
						eval(['set(handles.DATA' vs 'radiobutton,''Value'',1)'])
						no_active = 0;
						handles.ActiveVar = vs;
					end
					eval(['handles.DATA' vs 'checkbox=hhd;clear hhd'])
					if d < 4
						%E and V data
						% resample B
						if ~isempty(handles.BData{cl_id})
							c_log('proc','resampling B GFM')
							c_eval(...
								['data.B = c_resamp(handles.BData{cl_id},' vs ');'],...
								cl_id)
						elseif ~isempty(handles.BPPData{cl_id})
							c_log('proc','resampling B PP')
							c_eval(...
								['data.B = c_resamp(handles.BPPData{cl_id},' vs ');'],...
								cl_id)
						else, data.B = [];
						end
						handles.Data = [handles.Data, {data}];
					else
						%AUX data
						cur_aux = 1;
						handles.AUXData = [handles.AUXData, {data}];
					end
					eval(['clear ' vs])
				else
					c_log('load',['cannot load ' vs])
				end
			end
		end
		if ndata==0, set(hhh,'Enable','off')
		else, set(hhh,'Value',1)
		end
		eval(['handles.' s 'checkbox=hhh;clear hhh'])
		
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
		handles.EFWoffset = [handles.EFWoffset {offset}];
		
		% load CIS offsets, must be only in Z.
		offset = [0 0 0];
		handles.CISHoffset = [handles.CISHoffset {offset}];
		handles.CISHoffset = [handles.CISCoffset {offset}];
	end
	
	%create Axes
	handles.Xaxes = axes('Position',[pxa pya+(ha+dya)*3 wa ha],'Tag','Xaxes');
	handles.Yaxes = axes('Position',[pxa pya+(ha+dya)*2 wa ha],'Tag','Yaxes');
	handles.Zaxes = axes('Position',[pxa pya+(ha+dya)*1 wa ha],'Tag','Zaxes');
	handles.AUXaxes = axes('Position',[pxa pya wa ha],'Tag','AUXaxes');
	
	%create Calibrators
	wp = .15;
	handles.DXpanel = uipanel('Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*3 wp ha]);
	handles.DYpanel = uipanel('Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*2 wp ha]);
	handles.DZpanel = uipanel('Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*1 wp ha]);
	%handles.AUXpanel = uipanel('Position',[pxa+wa+dya*2 pya wp ha]);
	set(handles.DXpanel,'Title','dX')
	set(handles.DYpanel,'Title','dY')
	set(handles.DZpanel,'Title','dZ')
	handles.DXslider = uicontrol(handles.DXpanel,'Style','slider',...
		'Units','normalized','Position',[0.1 0.7 0.8 0.2],...
		'Max',5,'Min',-5,'Value',0,'Callback','c_cal_gui(''update_DXslider'')');
	handles.DXedit = uicontrol(handles.DXpanel,'Style','edit',...
		'Units','normalized','Position',[0.1 0.1 0.4 0.2],...
		'String','0.0','BackgroundColor',active_color,...
		'Callback','c_cal_gui(''update_DXedit'')','Tag','DXedit');
	handles.DXcheckbox = uicontrol(handles.DXpanel,'Style','checkbox',...
		'Units','normalized','Position',[0.6 0.1 0.4 0.2],...
		'String','Lock',...
		'Callback','c_cal_gui(''update_DXcheckbox'')','Tag','DXcheckbox');

	guidata(h0,handles);
	c_cal_gui('replot_all')
	av_figmenu
case 'replot_all'
	%plot data
	handles = guidata(h0);
	handles.DataList = {};
	handles.DataLegList = {};
	handles.AUXList = {};
	handles.AUXLegList = {};
	if ~length(handles.Data), error('no data to calibrate'),end
	h = [handles.Xaxes handles.Yaxes handles.Zaxes handles.AUXaxes];
	for ax=1:3, cla(h(ax)), end
	for j=1:length(handles.Data)
		data = handles.Data{j};
		if data.visible
			p_data = data.data;
			if data.editable 
				switch data.type
				case 'E'
					ofs = handles.EFWoffset{data.cl_id};
					p_data = corrDSIOffsets(p_data,real(ofs(1)),imag(ofs(1)),ofs(2));
				case 'V'
					if strcmp(data.sen,'COD')
						offset = handles.CISCoffset{data.cl_id};
					elseif strcmp(data.sen,'HIA')
						offset = handles.CISHoffset{data.cl_id};
					else, offset = [0 0 0];
					end
					p_data = corr_v_velocity(p_data,handles.EFWoffset{data.cl_id});
				otherwise
					disp('Unknown type.')
				end
			end
			% calculate Ez, convert V->E and V->E
			switch data.type
			case 'E'
				if any(p_data(:,4))==0,
					if ~isempty(data.B)
						[p_data,angle]=av_ed(p_data,data.B,handles.ang_limit);
						ii = find(abs(angle) < handles.ang_limit);
						if length(ii) > 1, p_data(ii,4) = p_data(ii,4)*NaN; end
					end
				end
				if ~handles.mode
					if ~isempty(data.B)
						p_data = av_e_vxb(p_data,data.B,-1);
					else
						handles.Data{j}.visible = 0;
						continue
					end
				end
			case 'V'
				if handles.mode
					if ~isempty(data.B)
						p_data = av_t_appl(av_cross(p_data,data.B),'*(-1e-3)');
					else
						handles.Data{j}.visible = 0;
						continue
					end
				end
			otherwise
				disp('Unknown type.')
			end
			handles.DataList = [handles.DataList {data.name}];
			handles.DataLegList = [handles.DataLegList {data.label}];
			for ax=1:3
                %disp([epoch2iso(p_data(1,1)) ' ' data.name])
				hold(h(ax),'on')
				handles.Data{j}.ploth(ax) = plot(h(ax),p_data(:,1),p_data(:,1+ax),...
					p_style(data.cl_id, data.inst, data.sen));
				hold(h(ax),'off')
			end
		end
	end
	labs = ['x' 'y' 'z'];
	if handles.mode
		u_s = 'E';
		u_u = 'mV/m';
	else
		u_s = 'V';
		u_u = 'km/s';
	end
	l_s = ['''' handles.DataLegList{1} ''''];
	if length(handles.DataLegList)>1
		for j=2:length(handles.DataLegList)
			l_s = [l_s ',''' handles.DataLegList{j} ''''];
		end
	end
	for ax=1:3
		ylabel(h(ax),[u_s '_' labs(ax) ' [' u_u ']'])
		grid(h(ax),'on')
		eval(['hxxx=legend(h(ax),' l_s ',''Location'',''NorthEastOutside'');'])
		set(hxxx,'FontSize',7);
	end
	for j=1:length(handles.AUXData)
		data = handles.AUXData{j};
		if data.visible
			p_data = data.data;
			handles.AUXList = [handles.AUXList {data.name}];
			handles.AUXLegList = [handles.AUXLegList {data.label}];
			hold(h(4),'on')
			handles.AUXData{j}.ploth = plot(h(4),p_data(:,1),p_data(:,2),...
				p_style(data.cl_id, data.inst, data.sen));
			hold(h(4),'off')
		end
	end
	l_s = ['''' handles.AUXLegList{1} ''''];
	if length(handles.AUXLegList)>1
		for j=2:length(handles.AUXLegList)
			l_s = [l_s ',''' handles.AUXLegList{j} ''''];
		end
	end
	eval(['hxxx=legend(h(4),' l_s ',''Location'',''NorthEastOutside'');'])
	set(hxxx,'FontSize',7);
	ylabel(h(4),'AUX')
	axis(h(4)); add_timeaxis; grid on
	av_zoom(handles.tlim,'x',h);
	guidata(h0,handles);
case 'update_DXslider'
	handles = guidata(h0);
	set(handles.DXedit,'String',...
    num2str(get(handles.DXslider,'Value')));
case 'update_DXedit'
	handles = guidata(h0);
	val = str2double(get(handles.DXedit,'String'));
	% Determine whether val is a number between 0 and 1
	if isnumeric(val) & length(val)==1 & ...
		val >= get(handles.DXslider,'Min') & ...
		val <= get(handles.DXslider,'Max')
		set(handles.DXslider,'Value',val);
	else
		% Increment the error count, and display it
		set(handles.DXedit,'String','You have entered an invalid entry ');
	end
case 'update_DXcheckbox'
	handles = guidata(h0);
	if get(handles.DXcheckbox,'Value')==1
		%lock
		set(handles.DXedit,'Enable','off','BackgroundColor',inactive_color)
		set(handles.DXslider,'Enable','off')
	else
		%unclock
		set(handles.DXedit,'Enable','on','BackgroundColor',active_color)
		set(handles.DXslider,'Enable','on')
		set(handles.DXslider,'Value',str2double(get(handles.DXedit,'String')))
	end
case 'update_DATAcheckbox'
	handles = guidata(h0);
	if get(eval(['handles.DATA' vs 'checkbox']),'Value')==1
		%need to plot the varible
		f_ok = 0;
		for j=1:length(handles.Data)
			if strcmp(handles.Data{j}.name,vs)
				f_ok = 1;
				handles.Data{j}.visible = 1;
			end
		end
		if ~f_ok
			for j=1:length(handles.AUXData)
				if strcmp(handles.AUXData{j}.name,vs)
					f_ok = 1;
					handles.AUXData{j}.visible = 1;
				end
			end
		end
		guidata(h0,handles);
		c_cal_gui('replot_all')
	else
		%need to hide the varible
		f_ok = 0;
		for j=1:length(handles.Data)
			if strcmp(handles.Data{j}.name,vs)
				f_ok = 1;
				handles.Data{j}.visible = 0;
				delete(handles.Data{j}.ploth)
			end
		end
		if ~f_ok
			for j=1:length(handles.AUXData)
				if strcmp(handles.AUXData{j}.name,vs)
					f_ok = 1;
					handles.AUXData{j}.visible = 0;
					delete(handles.AUXData{j}.ploth)
				end
			end
		end
		guidata(h0,handles);
	end
otherwise 
	disp('wrong action')
end

function out = corr_v_velocity(v,offset)
out = data;
out(:,2) = v(:,2) - offset(1);
out(:,3) = v(:,3) - offset(2);
out(:,4) = v(:,4) - offset(3);

function out = p_style(cl_id, inst, sen)
out = '';
colrs = ['k','r','g','b'];
switch inst
case 'EFW'
	switch sen
	case 'all'
		out = '-';
	case 'p12'
		out = '--';
	case 'p34'
		out = '--';	
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
out = [out colrs(cl_id)];
