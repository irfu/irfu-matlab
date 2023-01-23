function c_cal_gui(varargin)
%C_CAL_GUI EFW calibration GUI
%
% This program will load EFW, CIS and EDI data for present in the
% current directory and allow the user to interactively determine
% offsets in EFW and CIS Vz.
%
% Example:
%	cd /tmp/my_event/20020304_1000
%	c_cal_gui
%
% See also ClusterProc/corrSOffsetM
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

persistent h0;
persistent inprog_mtx; % this is a kinda mutex to hold "in progress" status

sp = '.';

% plotting constants
inactive_color = [0.831373 0.815686 0.784314]; % Gray
active_color = [1 1 1]; % White
active_p_color = 'magenta';
pxa = .07+.12; pya = .1; wa = .47; ha = .215; dya = .01; % Positions
MM = [-100 100; -5 5]; % Min and max value for sliders
MMZ = [-200 200; .5 2.];
main_fig_id = 23;
raw_fig_id = 24;
spect_fig_id = 25;
% positions for 1024x768, 1600x1200
pos_main_fig   = [[  2  40 990 640];[   2 159 990 945]];
pos_raw_fig    = [[994 160 560 420];[1002 160 560 420]];
pos_spect_fig  = [[464 282 560 420];[1002 655 560 420]];

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
elseif regexp(action,'^update_D[X-Z]slider$')
  comp_s = action(9);
  switch comp_s
    case 'X'
      comp = 1;
    case 'Y'
      comp = 2;
    case 'Z'
      comp = 3;
    otherwise
      disp('strange component...')
      comp = 0;
  end
  action = 'update_Dslider';
elseif regexp(action,'^update_D[X-Z]edit$')
  comp_s = action(9);
  switch comp_s
    case 'X'
      comp = 1;
    case 'Y'
      comp = 2;
    case 'Z'
      comp = 3;
    otherwise
      disp('strange component...')
      comp = 0;
  end
  action = 'update_Dedit';
elseif regexp(action,'^update_D[X-Z]checkbox$')
  comp_s = action(9);
  switch comp_s
    case 'X'
      comp = 1;
    case 'Y'
      comp = 2;
    case 'Z'
      comp = 3;
    otherwise
      disp('strange component...')
      comp = 0;
  end
  action = 'update_Dcheckbox';
elseif regexp(action,'^click_[W-Z]axes$')
  curr_ax = action(7);
  if curr_ax=='W', curr_ax = 'AUX'; end
  action = 'click_axes';
end

switch action
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % init
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'init'
    inprog_mtx = 0;
    
    % Create figure
    if find(get(0,'children')==main_fig_id)
      pos_old = get(main_fig_id,'Position');
    else, pos_old = [];
    end
    h0 = figure(main_fig_id);
    clf
    set(main_fig_id,'Name', 'CLUSTER CAL GUI')
    
    hnd = guihandles(h0);
    
    % Save the screen size
    sc_s = get(0,'ScreenSize');
    if sc_s(3)==1600 && sc_s(4)==1200, hnd.scrn_size = 2;
    else, hnd.scrn_size = 1;
    end
    
    if isempty(pos_old), set(h0,'Position', pos_main_fig(hnd.scrn_size,:)), end
    
    hnd.DATApanel = uipanel('Position',[.85 .01 .14 .98]);
    set(hnd.DATApanel,'Title','Data')
    
    % Load data
    
    c_Edata = {'diE?p1234','diEs?p12','diEs?p32','diEs?p34', 'diELXs?p42'};
    c_Ddata = {'diEDI?'};
    c_Vdata = {'diVCp?','diVCh?'};
    c_AUXdata = {'P?','NCp?','NCh?'};
    dd={c_Edata c_Ddata c_Vdata c_AUXdata};
    %d_s = {'Edata', 'Ddata', 'Vdata', 'AUXdata'};
    %for d=1:4, eval(['hnd.' d_s{d} '= {};']),end
    hnd.Data = {};
    hnd.DataList = {};
    hnd.AUXList = {};
    hnd.DATArbList = {};
    hnd.BPPData = {};
    hnd.BData = {};
    hnd.EFWoffset = zeros(4,2);
    hnd.CISHoffset = zeros(4,2);
    hnd.CISCoffset = zeros(4,2);
    hnd.EFWoffset_save = zeros(4,2);
    hnd.CISHoffset_save = zeros(4,2);
    hnd.CISCoffset_save = zeros(4,2);
    hnd.mode = 1; % 0 is for V, 1 is for E.
    % hnd.last
    % empty means replotting all,
    % otherwise contains a cell array of data in a plot queue
    hnd.last = [];
    hnd.off = [0+0i 0];
    hnd.ang_limit = 15; % 15 degrees.
    hnd.tlim = [0 0];
    hnd.ts_marker = [];
    hnd.te_marker = [];
    hnd.c_visible = [1 1 1 1]; % all sc are visible by default
    no_active = 1;
    ncdata = 0; % number of non-AUX variables
    
    for cl_id=1:4
      s = irf_ssub('C?',cl_id);
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
      if c_load(irf_ssub('diB?',cl_id))
        irf_log('load',irf_ssub('loaded  diB?',cl_id))
        c_eval('BFGM=diB?;clear diB?',cl_id)
      else
        if c_load(irf_ssub('diBPP?',cl_id))
          irf_log('load',irf_ssub('loaded  diB?',cl_id))
          c_eval('BFGM=diBPP?;clear diBPP?',cl_id)
        end
      end
      % Load resampled FGM data
      [ok,Br] = c_load('diBr?',cl_id);
      if ok, irf_log('load',irf_ssub('loaded  diBr?',cl_id)), end
      [ok,Brs] = c_load('diBrs?',cl_id);
      if ok, irf_log('load',irf_ssub('loaded  diBrs?',cl_id)), end
      
      hnd.BData = [hnd.BData {BFGM}];
      hnd.BPPData = [hnd.BPPData {BPP}];
      clear BFGM BPP
      
      % Load data
      for d=1:4
        for j=1:length(dd{d})
          vs = irf_ssub(dd{d}{j},cl_id);
          if c_load(vs)
            ndata = ndata + 1;
            irf_log('load',['loaded ' vs])
            %{
					if eval(['any(any(isnan(' vs '(:,2:end))))']),
						irf_log('load',[vs ' contains no data'])
						continue
					end
            %}
            dsc = c_desc(vs);
            clear data
            data.name = vs;
            data.cl_id = str2double(dsc.cl_id);
            data.inst = dsc.inst;
            data.sig = dsc.sig;
            data.label = [dsc.inst ' ' dsc.sig];
            senStr = dsc.sen;
            if length(vs)>5 && strcmpi(vs(4:5),'LX'), senStr = [senStr ' LX']; end
            if dsc.sen, data.label = [data.label ' (' senStr ')']; end
            data.sen = dsc.sen;
            [data.plot_style,data.plot_color] =...
              p_style(data.cl_id, data.inst,data.sen);
            
            % tlim
            eval(['tlxxx(1) = ' vs '(1,1); tlxxx(2) = ' vs '(end,1);'])
            if ~hnd.tlim(1), hnd.tlim(1) = tlxxx(1);
            elseif hnd.tlim(1) > tlxxx(1), hnd.tlim(1) = tlxxx(1);
            end
            if hnd.tlim(2) < tlxxx(2), hnd.tlim(2) = tlxxx(2); end
            clear tlxxx
            
            eval(['data.data =' vs ';'])
            data.p_data = [];
            if ((strcmp(dsc.sen,'1234') || strcmp(dsc.sen,'12') || strcmp(dsc.sen,'32')) && ...
                strcmp(dsc.sig,'E')) || strcmp(dsc.sen,'COD') || ...
                strcmp(dsc.sig,'N') || strcmp(dsc.sig,'Np')
              data.visible = 0;
            else
              data.visible = 1;
            end
            if strcmp(dsc.sig,'E'), data.type = 'E';
            elseif strcmp(dsc.sig,'V') || strcmp(dsc.sig,'Vp')
              data.type = 'V';
            else, data.type = 'A';
            end
            if (strcmp(dsc.inst,'EFW') && strcmp(data.type,'E')) || ...
                (strcmp(dsc.inst,'CIS') && (strcmp(data.type,'V')))
              data.editable = 1;
            else, data.editable = 0;
            end
            if d==1 || d==3
              hhd = uicontrol(hnd.DATApanel,'Style','radiobutton',...
                'Units','normalized',...
                'Position',[pxd pyd-.02*ndata 0.1 .02],...
                'String','',...
                'Callback',['c_cal_gui(''update_DATA' vs 'radiobutton'')'],...
                'Tag',['DATA' vs 'radiobutton']); %#ok<NASGU>
              hnd.DATArbList = [hnd.DATArbList {vs}];
              eval(['hnd.DATA' vs 'radiobutton=hhd;clear hhd'])
            end
            hhd = uicontrol(hnd.DATApanel,'Style','checkbox',...
              'Units','normalized',...
              'Position',[pxd+.1 pyd-.02*ndata 0.8 .02],...
              'String',data.label,...
              'Callback',['c_cal_gui(''update_DATA' vs 'checkbox'')'],...
              'Tag',['DATA' vs 'checkbox']);
            if data.visible, set(hhd,'Value',1), end
            if no_active && data.visible && ...
                strcmp(dsc.inst,'EFW') && strcmp(dsc.sig,'E')
              eval(['set(hnd.DATA' vs 'radiobutton,''Value'',1)'])
              no_active = 0;
              hnd.ActiveVar = vs;
              hnd.old_ActiveVar = '';
            end
            eval(['hnd.DATA' vs 'checkbox=hhd;clear hhd'])
            
            % Delta offsets: remove automatic and apply CAA
            if d==1 && strcmp(dsc.quant,'dies')
              eval(['Del_caa = c_efw_delta_off(' vs '(1,1),cl_id);'])
              if ~isempty(Del_caa)
                [ok,Delauto] = c_load('D?p12p34',cl_id);
                if ~ok || isempty(Delauto)
                  irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cl_id))
                else
                  eval([vs ' = caa_corof_delta(' vs ',str2double(dsc.sen),Delauto,''undo'');'])
                  eval([vs ' = caa_corof_delta(' vs ',str2double(dsc.sen),Del_caa,''apply'');'])
                end
              end
              clear Del_caa
            end
            
            % Assign magnetic field to the variable
            % For EFW we use already resampled filed
            % for the others we resample.
            data.B = [];
            if d==1
              if strcmp(vs(1:4),'diEs')
                if ~isempty(Brs)
                  irf_log('proc',irf_ssub('using Brs?',cl_id))
                  data.B = Brs;
                end
              else
                if ~isempty(Br)
                  irf_log('proc',irf_ssub('using Br?',cl_id))
                  data.B = Br;
                end
              end
            end
            if d < 4 && isempty(data.B)
              % E and V data
              % Resample B
              if ~isempty(Brs)
                irf_log('proc','interpolating Brs')
                warning('off','MATLAB:interp1:NaNinY')
                c_eval(...
                  ['data.B = irf_resamp(Brs,' vs ');'],...
                  cl_id)
                warning('on','MATLAB:interp1:NaNinY')
              elseif ~isempty(hnd.BPPData{cl_id})
                irf_log('proc','interpolating B PP')
                c_eval(...
                  ['data.B = irf_resamp(hnd.BPPData{cl_id},' vs ');'],...
                  cl_id)
              elseif ~isempty(hnd.BData{cl_id})
                irf_log('proc','resampling B GFM')
                c_eval(...
                  ['data.B = irf_resamp(hnd.BData{cl_id},' vs ');'],...
                  cl_id)
              else, data.B = [];
              end
            end
            
            % Assing AUX flag
            if d < 4
              data.aux = 0;
              ncdata = ncdata + 1;
            else
              % AUX data
              data.aux = 1;
            end
            
            hnd.Data = [hnd.Data, {data}];
            eval(['clear ' vs])
            
          else
            irf_log('load',['cannot load ' vs])
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
      offset = c_ctl(cl_id,'dsiof');
      
      if isempty(offset)
        [dsiof_def, dam_def] = c_efw_dsi_off(hnd.tlim(1),cl_id);
        offset(1) = dsiof_def; offset(2) = dam_def;
        
        [ok1,Ddsi] = c_load('Ddsi?',cl_id); if ok1, offset(1) = Ddsi; end
        [ok2,Damp] = c_load('Damp?',cl_id); if ok2, offset(2) = Damp; end
        
        if ok1 || ok2, irf_log('calb','Using saved DSI offsets')
        else, irf_log('calb','Using default DSI offsets')
        end
        clear dsiof_def dam_def Ddsi Damp
      else
        irf_log('calb','Using user specified DSI offsets')
      end
      
      cd(old_pwd)
      hnd.EFWoffset(cl_id,:) = offset;
      
      % Load CIS offsets, must be only in Z.
      warning('off','MATLAB:load:variableNotFound')
      if exist('mCIS.mat','file')
        c_eval('load -mat mCIS DHdsi? DCdsi? DHz? DCz?',cl_id)
      end
      warning('on','MATLAB:load:variableNotFound')
      
      offset = [0 0];
      if exist(irf_ssub('DHdsi?',cl_id),'var')
        c_eval('offset(1)=DHdsi?;clear DHdsi?',cl_id)
        irf_log('load',...
          sprintf('HIA xy offset DHdsi%d = %.2f %.2f*i',...
          cl_id,real(offset(1)),imag(offset(1))))
      end
      if exist(irf_ssub('DHz?',cl_id),'var')
        c_eval('offset(2)=DHz?;clear DHz?',cl_id)
        irf_log('load',...
          sprintf('HIA z offset DHz%d = %.2f',cl_id,offset(2)))
      end
      hnd.CISHoffset(cl_id,:) = offset;
      if exist(irf_ssub('DCdsi?',cl_id),'var')
        c_eval('offset(1)=DCdsi?;clear DCdsi?',cl_id)
        irf_log('load',...
          sprintf('CODIF xy offset DHdsi%d = %.2f %.2f*i',...
          cl_id,real(offset(1)),imag(offset(1))))
      end
      if exist(irf_ssub('DCz?',cl_id),'var')
        c_eval('offset(2)=DCz?;clear DCz?',cl_id)
        irf_log('load',...
          sprintf('CODIF z offset DCz%d = %.2f',cl_id,offset(2)))
      end
      hnd.CISCoffset(cl_id,:) = offset;
    end
    
    hnd.EFWoffset_save = hnd.EFWoffset;
    hnd.CISHoffset_save = hnd.CISHoffset;
    hnd.CISCoffset_save = hnd.CISCoffset;
    hnd.off_updated = 1;
    
    % Check if we have any data apart from AUX
    if ~ncdata, error('No usefull data loaded'), end
    
    % Create Data Axes
    hnd.Xaxes = axes('Position',[pxa pya+(ha+dya)*3 wa ha],'Tag','Xaxes',...
      'ButtonDownFcn','c_cal_gui(''click_Xaxes'')');
    hnd.Yaxes = axes('Position',[pxa pya+(ha+dya)*2 wa ha],'Tag','Yaxes',...
      'ButtonDownFcn','c_cal_gui(''click_Yaxes'')');
    hnd.Zaxes = axes('Position',[pxa pya+(ha+dya)*1 wa ha],'Tag','Zaxes',...
      'ButtonDownFcn','c_cal_gui(''click_Zaxes'')');
    hnd.AUXaxes = axes('Position',[pxa pya wa ha],'Tag','AUXaxes',...
      'ButtonDownFcn','c_cal_gui(''click_Waxes'')');
    
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
    set(hnd.DXpanel,'Title','dEX','TitlePosition','centertop')
    set(hnd.DYpanel,'Title','dEY','TitlePosition','centertop')
    set(hnd.DZpanel,'Title','dAMP','TitlePosition','centertop')
    hnd.DXslider = uicontrol(hnd.DXpanel,'Style','slider',...
      'Units','normalized','Position',[0.1 0.7 0.8 0.2],...
      'Min',MM(2,1),'Max',MM(2,2),'Value',0,'Callback','c_cal_gui(''update_DXslider'')');
    hnd.DXedit = uicontrol(hnd.DXpanel,'Style','edit',...
      'Units','normalized','Position',[0.1 0.1 0.4 0.2],...
      'String','0.0','BackgroundColor',active_color,...
      'Callback','c_cal_gui(''update_DXedit'')','Tag','DXedit');
    hnd.DXcheckbox = uicontrol(hnd.DXpanel,'Style','checkbox',...
      'Units','normalized','Position',[0.6 0.1 0.4 0.2],...
      'String','Lock','Enable','off',...
      'Callback','c_cal_gui(''update_DXcheckbox'')','Tag','DXcheckbox');
    hnd.DYslider = uicontrol(hnd.DYpanel,'Style','slider',...
      'Units','normalized','Position',[0.1 0.7 0.8 0.2],...
      'Min',MM(2,1),'Max',MM(2,2),'Value',0,'Callback','c_cal_gui(''update_DYslider'')');
    hnd.DYedit = uicontrol(hnd.DYpanel,'Style','edit',...
      'Units','normalized','Position',[0.1 0.1 0.4 0.2],...
      'String','0.0','BackgroundColor',active_color,...
      'Callback','c_cal_gui(''update_DYedit'')','Tag','DYedit');
    hnd.DYcheckbox = uicontrol(hnd.DYpanel,'Style','checkbox',...
      'Units','normalized','Position',[0.6 0.1 0.4 0.2],...
      'String','Lock','Enable','off',...
      'Callback','c_cal_gui(''update_DYcheckbox'')','Tag','DYcheckbox');
    hnd.DZslider = uicontrol(hnd.DZpanel,'Style','slider',...
      'Units','normalized','Position',[0.1 0.7 0.8 0.2],...
      'Min',MMZ(2,1),'Max',MMZ(2,2),'Value',1,'Callback','c_cal_gui(''update_DZslider'')');
    hnd.DZedit = uicontrol(hnd.DZpanel,'Style','edit',...
      'Units','normalized','Position',[0.1 0.1 0.4 0.2],...
      'String','0.0','BackgroundColor',active_color,...
      'Callback','c_cal_gui(''update_DZedit'')','Tag','DZedit');
    hnd.DZcheckbox = uicontrol(hnd.DZpanel,'Style','checkbox',...
      'Units','normalized','Position',[0.6 0.1 0.4 0.2],...
      'String','Lock','Enable','off',...
      'Callback','c_cal_gui(''update_DZcheckbox'')','Tag','DZcheckbox');
    hbut = .05;
    hnd.SaveALLbutton = uicontrol(h0,'Style','pushbutton',...
      'Units','normalized','Position',[pxa+wa+dya*2+.015 .01 wp hbut],...
      'String','Save all',...
      'Callback','c_cal_gui(''press_SaveALLbutton'')','Tag','SaveALLbutton');
    hnd.SAVEbutton = uicontrol(h0,'Style','pushbutton',...
      'Units','normalized',...
      'Position',[pxa+wa+dya*2+.015+.55*wp pya+(ha+dya)*1-.01-hbut 0.45*wp hbut],...
      'String','Save',...
      'Callback','c_cal_gui(''press_SAVEbutton'')','Tag','SAVEbutton');
    hnd.RESETbutton = uicontrol(h0,'Style','pushbutton',...
      'Units','normalized',...
      'Position',[pxa+wa+dya*2+.015 pya+(ha+dya)*1-.01-hbut 0.45*wp hbut],...
      'String','Reset',...
      'Callback','c_cal_gui(''press_RESETbutton'')','Tag','RESETbutton');
    hnd.FWDbutton = uicontrol(h0,'Style','pushbutton',...
      'Units','normalized','Position',[pxa+wa/2+.015 .01 .45*wp hbut*.7],...
      'String','>>',...
      'Callback','c_cal_gui(''press_FWDbutton'')','Tag','FWDbutton');
    hnd.RWDbutton = uicontrol(h0,'Style','pushbutton',...
      'Units','normalized','Position',[pxa+wa/2-.45*wp-.015 .01 .45*wp hbut*.7],...
      'String','<<',...
      'Callback','c_cal_gui(''press_RWDbutton'')','Tag','RWDbutton');
    
    % Disable buttons
    set(hnd.SaveALLbutton,'Enable','off');
    set(hnd.SAVEbutton,'Enable','off');
    set(hnd.RESETbutton,'Enable','off');
    
    hnd.EVbutton = uicontrol(h0,'Style','togglebutton',...
      'Units','normalized','Position',[pxa+wa+dya*2+.015 pya wp hbut],...
      'String','V, ExB',...
      'Callback','c_cal_gui(''update_EVbutton'')','Tag','EVbutton');
    
    % Menu
    hnd.menu_time = uimenu(h0,'Label','&Time');
    hnd.menu_zoom_in = uimenu(hnd.menu_time,'Label','&Zoom in',...
      'Callback','c_cal_gui(''zoom_in'')',...
      'Accelerator','t',...
      'Enable','off');
    hnd.menu_zoom_un = uimenu(hnd.menu_time,'Label','&Undo zoom',...
      'Callback','c_cal_gui(''zoom_un'')',...
      'Accelerator','z',...
      'Enable','off');
    hnd.menu_zoom_rs = uimenu(hnd.menu_time,'Label','&Reset',...
      'Callback','c_cal_gui(''zoom_rs'')',...
      'Accelerator','r',...
      'Enable','off');
    
    % Menu tools
    hnd.menu_tools = uimenu(h0,'Label','&Tools');
    hnd.menu_show_raw = uimenu(hnd.menu_tools,'Label','&Show raw data',...
      'Callback','c_cal_gui(''show_raw'')',...
      'Accelerator','d');
    hnd.menu_show_b = uimenu(hnd.menu_tools,'Label','&Show B',...
      'Callback','c_cal_gui(''show_b'')',...
      'Accelerator','b');
    hnd.menu_show_spect = uimenu(hnd.menu_tools,'Label','&Spectrum',...
      'Callback','c_cal_gui(''show_spect'')',...
      'Accelerator','f',...
      'Enable','on');
    hnd.menu_cut_int = uimenu(hnd.menu_tools,'Label','&Cut the interval',...
      'Callback','c_cal_gui(''cut_int'')',...
      'Accelerator','k',...
      'Enable','off');
    hnd.menu_splot = uimenu(hnd.menu_tools,'Label','&Summary plots',...
      'Callback','caa_pl_summary(''noexport'')',...
      'Accelerator','m',...
      'Enable','on');
    hnd.menu_splot_im = uimenu(hnd.menu_tools,'Label','&S-plots+images',...
      'Callback','caa_pl_summary',...
      'Accelerator','i',...
      'Enable','on');
    
    guidata(h0,hnd);
    
    c_cal_gui('replot_all')
    
    % Take care of the active variable
    c_cal_gui('ch_active_var')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ch_active_var
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'ch_active_var'
    hnd = guidata(h0);
    if isempty(hnd.old_ActiveVar)
      % Initialize
      d_m = 2; vs = 'E'; vsz = 'AMP';
      j = D_findByName(hnd.Data,hnd.ActiveVar);
      hnd.off = hnd.EFWoffset(hnd.Data{j}.cl_id,:);
    else
      j = D_findByName(hnd.Data,hnd.ActiveVar);
      
      if strcmp(hnd.Data{j}.type,'E')
        d_m = 2; vs = 'E'; vsz = 'AMP';
        hnd.off = hnd.EFWoffset(hnd.Data{j}.cl_id,:);
        set(hnd.DXslider,'Enable','on')
        set(hnd.DXedit,'Enable','on')
        set(hnd.DYslider,'Enable','on')
        set(hnd.DYedit,'Enable','on')
        set(hnd.menu_show_spect,'Enable','on')
      else % Display V, offsets in CIS
        d_m = 1; vs = 'V'; vsz = 'Vz';
        if strcmp(hnd.Data{j}.sen,'HIA')
          hnd.off = hnd.CISHoffset(hnd.Data{j}.cl_id,:);
        else, hnd.off = hnd.CISCoffset(hnd.Data{j}.cl_id,:);
        end
        set(hnd.DXslider,'Enable','off')
        set(hnd.DXedit,'Enable','off')
        set(hnd.DYslider,'Enable','off')
        set(hnd.DYedit,'Enable','off')
        set(hnd.menu_show_spect,'Enable','off')
      end
    end
    set(hnd.DXedit,'String',num2str(real(hnd.off(1))))
    set(hnd.DXslider,'Min',MM(d_m,1),'Max',MM(d_m,2),'Value',real(hnd.off(1)))
    set(hnd.DXpanel,'Title',[num2str(MM(d_m,1)) ' < d' vs 'x < ' num2str(MM(d_m,2))])
    set(hnd.DYedit,'String',num2str(imag(hnd.off(1))))
    set(hnd.DYslider,'Min',MM(d_m,1),'Max',MM(d_m,2),'Value',imag(hnd.off(1)))
    set(hnd.DYpanel,'Title',[num2str(MM(d_m,1)) ' < d' vs 'y < ' num2str(MM(d_m,2))])
    set(hnd.DZedit,'String',num2str(hnd.off(2)))
    set(hnd.DZslider,'Min',MMZ(d_m,1),'Max',MMZ(d_m,2),'Value',hnd.off(2))
    set(hnd.DZpanel,'Title',[num2str(MMZ(d_m,1)) ' < d' vsz ' < ' num2str(MMZ(d_m,2))])
    
    guidata(h0,hnd);
    c_cal_gui('update_SAVEbuttons')
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
    
    if inprog_mtx, return, end
    inprog_mtx = 1;
    
    % Plot data
    hnd = guidata(h0);
    h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
    if isempty(hnd.Data), error('no data to calibrate'),end
    if isempty(hnd.last), disp('List empty, will do nothing'), return, end
    
    % Update selected records in DataLegList
    % we set hnd.last to name of the variable we need to add
    d_ii = D_findByNameList(hnd.Data,hnd.last);
    
    % Sanity check
    if isempty(d_ii)
      error('cannot find variables in the data list')
    end
    
    % Plotting
    for j=1:length(d_ii)
      %disp([action ': plotting ' hnd.Data{d_ii(j)}.name])
      
      % Process only visible data
      if ~hnd.Data{d_ii(j)}.visible, continue, end
      
      if hnd.Data{d_ii(j)}.aux
        % Sanity check
        if isempty(L_find(hnd.AUXList,hnd.Data{d_ii(j)}.name))
          hnd.AUXList = L_add(hnd.AUXList,hnd.Data{d_ii(j)}.name);
        else
          disp([action ': ' hnd.Data{d_ii(j)}.name ' already in the list'])
        end
        
        % Plotting
        hold(h(4),'on')
        hnd.Data{d_ii(j)}.ploth = plot(hnd.AUXaxes,...
          hnd.Data{d_ii(j)}.data(:,1),hnd.Data{d_ii(j)}.data(:,2),...
          hnd.Data{d_ii(j)}.plot_style,...
          'Color',hnd.Data{d_ii(j)}.plot_color);
        hold(h(4),'off')
      else
        % Remove the variable if it is already in the list
        ii = L_find(hnd.DataList,hnd.Data{d_ii(j)}.name);
        if ii
          hnd.DataList = L_rm_ii(hnd.DataList,ii);
          try
            delete(hnd.Data{d_ii(j)}.ploth)
          catch err
            if ~strcmp(err.identifier,'MATLAB:hg:udd_interface:CannotDelete')
              rethrow(err)
            end
          end
        end
        
        % Update p_data only if calibrations were changed
        if isempty(hnd.Data{d_ii(j)}.p_data) || hnd.off_updated
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
            
            % Take care of the active variable
            if strcmp(hnd.Data{d_ii(j)}.name,hnd.ActiveVar)
              set(hnd.Data{d_ii(j)}.ploth(ax),...
                'LineWidth',get(hnd.Data{d_ii(j)}.ploth(ax),'LineWidth')*4,...
                'Color',active_p_color)
            end
            hold(h(ax),'off')
          end
          %{
				xx = hnd.Data{d_ii(j)}.ploth;
				disp(sprintf('%s after: %f %f %f',hnd.Data{d_ii(j)}.name,xx(1),xx(2),xx(3)))
          %}
        end
      end
    end
    hnd.last = [];
    hnd.off_updated = 0;
    
    irf_zoom(h,'x',hnd.tlim(end,:));
    
    % Hide the markers so thet they will not contribute YLim
    ts_tmp = hide_t_marker(hnd,hnd.ts_marker);
    te_tmp = hide_t_marker(hnd,hnd.te_marker);
    
    set(h,'YLimMode','auto')
    
    % Show the markers back
    if ~isempty(ts_tmp.t)
      hnd.ts_marker = replot_t_marker(hnd,ts_tmp);
    end
    if ~isempty(te_tmp.t)
      hnd.te_marker = replot_t_marker(hnd,te_tmp);
    end
    
    guidata(h0,hnd);
    c_cal_gui('update_legend')
    
    inprog_mtx = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % replot_all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'replot_all'
    % Replot all data
    hnd = guidata(h0);
    h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
    if isempty(hnd.Data), error('no data to calibrate'),end
    
    if isempty(hnd.last)
      % Clear everything
      for ax=1:4, cla(h(ax)), end
      hnd.DataList = {};
      hnd.AUXList = {};
      
      % Plot
      for j=1:length(hnd.Data)
        hnd.last = L_add(hnd.last,hnd.Data{j}.name);
        hnd.Data{j}.p_data = [];
      end
      guidata(h0,hnd);
      c_cal_gui('replot')
      hnd = guidata(h0);
      
      % Axes labels
      labs = ['x' 'y' 'z'];
      if hnd.mode, u_s = 'E'; u_u = 'mV/m'; cmp = '';
      else, u_s = 'V'; u_u = 'km/s'; cmp = '\perp ';
      end
      for ax=1:3
        ylabel(h(ax),[u_s '_{' cmp labs(ax) '} [' u_u ']'])
        grid(h(ax),'on')
      end
      ylabel(h(4),'AUX')
      irf_timeaxis(h);
      grid(h(4),'on')
      
      % Time span
      irf_zoom(h,'x',hnd.tlim(end,:));
      
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
          l_s = [l_s ',''' hnd.Data{ii(j)}.label '''']; %#ok<AGROW>
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
          l_s = [l_s ',''' hnd.Data{ii(j)}.label '''']; %#ok<AGROW>
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
    % update_SAVEbuttons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_SAVEbuttons'
    hnd = guidata(h0);
    ava = hnd.Data{D_findByName(hnd.Data,hnd.ActiveVar)};
    cl_id = ava.cl_id;
    
    % See if we have offsets different from the saved ones and
    % enable SAVE and RESET buttons
    if strcmp(ava.type,'E'), off_s = hnd.EFWoffset_save(cl_id,:);
    elseif strcmp(ava.type,'V') && strcmp(ava.sen,'HIA')
      off_s = hnd.CISHoffset_save(cl_id,:);
    elseif strcmp(ava.type,'V') && strcmp(ava.sen,'COD')
      off_s = hnd.CISCoffset_save(cl_id,:);
    else, disp([action ': must be something wrong']), off_s = [0 0];
    end
    hxx = [hnd.SAVEbutton hnd.RESETbutton hnd.SaveALLbutton];
    if any(off_s-hnd.off), for j=1:3, set(hxx(j),'Enable','on'), end
    else, for j=1:3, set(hxx(j),'Enable','off'), end
    end
    
    % SaveALL button
    if any(hnd.EFWoffset_save(:) - hnd.EFWoffset(:)) || ...
        any(hnd.CISHoffset_save(:) - hnd.CISHoffset(:)) ||...
        any(hnd.CISCoffset_save(:) - hnd.CISCoffset(:))
      set(hnd.SaveALLbutton,'Enable','on')
    else, set(hnd.SaveALLbutton,'Enable','off')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update_off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_off'
    
    if inprog_mtx, return, end
    inprog_mtx = 1;
    
    hnd = guidata(h0);
    k = D_findByName(hnd.Data,hnd.ActiveVar);
    
    if strcmp(hnd.Data{k}.type,'E')
      % We calibrate E
      
      % Sanity check
      if ~strcmp(hnd.Data{k}.inst,'EFW')
        disp('update_off: we are doing something wrong with INST')
        inprog_mtx = 0;
        return
      end
      
      % See if offsets were really changed and replot everything
      if any(hnd.EFWoffset(hnd.Data{k}.cl_id,:) - hnd.off)
        hnd.EFWoffset(hnd.Data{k}.cl_id,:) = hnd.off;
        %disp(sprintf('%s : offsets %f %f %f',action,real(hnd.off(1)),imag(hnd.off(1)),hnd.off(2)))
        hnd.off_updated = 1;
        
        % Create list of variables which need to be replotted
        ii = D_findByCLID(hnd.Data,hnd.Data{k}.cl_id);
        up_list = {};
        for j=1:length(ii)
          if strcmp(hnd.Data{ii(j)}.inst,'EFW') && strcmp(hnd.Data{ii(j)}.sig,'E')
            up_list = L_add(up_list,hnd.Data{ii(j)}.name);
          end
        end
        % Sanity check
        if isempty(up_list)
          disp('update_off: we are doing something wrong with UP_LIST')
          inprog_mtx = 0;
          return
        end
        hnd.last = up_list;
        guidata(h0,hnd);
        inprog_mtx = 0;
        c_cal_gui('replot')
      end
    else
      % We calibrate V
      
      % See if offsets were really changed and replot everything
      if strcmp(hnd.Data{k}.sen,'HIA')
        d_off = hnd.CISHoffset(hnd.Data{k}.cl_id,:);
      else, d_off = hnd.CISCoffset(hnd.Data{k}.cl_id,:);
      end
      if any(d_off - hnd.off)
        if strcmp(hnd.Data{k}.sen,'HIA')
          hnd.CISHoffset(hnd.Data{k}.cl_id,:) = hnd.off;
        else, hnd.CISCoffset(hnd.Data{k}.cl_id,:) = hnd.off;
        end
        %disp(sprintf('%s : offsets %f %f %f',action,real(hnd.off(1)),imag(hnd.off(1)),hnd.off(2)))
        hnd.off_updated = 1;
        
        % Create list of variables which need to be replotted
        ii = D_findByCLID(hnd.Data,hnd.Data{k}.cl_id);
        up_list = {};
        for j=1:length(ii)
          if strcmp(hnd.Data{ii(j)}.type,'V') && ...
              strcmp(hnd.Data{ii(j)}.sen,hnd.Data{k}.sen)
            up_list = L_add(up_list,hnd.Data{ii(j)}.name);
          end
        end
        % Sanity check
        if isempty(up_list)
          disp('update_off: we are doing something wrong with UP_LIST')
          inprog_mtx = 0;
          return
        end
        hnd.last = up_list;
        guidata(h0,hnd);
        inprog_mtx = 0;
        c_cal_gui('replot')
      end
    end
    c_cal_gui('update_SAVEbuttons')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update_EVbutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_EVbutton'
    hnd = guidata(h0);
    
    % Clear last for security
    hnd.last = [];
    if get(hnd.EVbutton,'Value')==1
      % Go to V mode
      hnd.mode = 0;
    else
      % Return to E mode
      hnd.mode = 1;
    end
    %{
	for j=1:length(hnd.DATArbList)
		hxxx = eval(['hnd.DATA' hnd.DATArbList{j} 'radiobutton']);
		if strcmp(get(hxxx,'Enable'),'off')
			set(hxxx,'Enable','on');
		else, set(hxxx,'Enable','off');
		end
	end
    %}
    guidata(h0,hnd);
    c_cal_gui('replot_all')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % press_SAVEbutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'press_SAVEbutton'
    hnd = guidata(h0);
    
    ava = hnd.Data{D_findByName(hnd.Data,hnd.ActiveVar)};
    cl_id = ava.cl_id;
    
    if strcmp(ava.inst,'EFW') && strcmp(ava.type,'E')
      f_name = './mEDSI.mat';
      c_eval('Ddsi?=hnd.off(1);Damp?=hnd.off(2);')
      hnd.EFWoffset_save(cl_id,:) = hnd.off;
      var_s = 'Ddsi? Damp?';
    elseif strcmp(ava.inst,'CIS')
      f_name = './mCIS.mat';
      if strcmp(ava.sen,'HIA')
        c_eval('DHdsi?=hnd.off(1);DHz?=hnd.off(2);')
        hnd.CISHoffset_save(cl_id,:) = hnd.off;
        var_s = 'DHdsi? DHz?';
      else
        c_eval('DCdsi?=hnd.off(1);DCz?=hnd.off(2);')
        hnd.CISCoffset_save(cl_id,:) = hnd.off;
        var_s = 'DCdsi? DCz?';
      end
    else
      disp(['Cannot save ' ava.name])
      return
    end
    if exist(f_name,'file'), apd = ' -append';
    else, apd = '';
    end
    c_eval(['save ' f_name ' ' var_s ' -mat' apd],cl_id)
    irf_log('save',irf_ssub([var_s ' -> ' f_name],cl_id))
    guidata(h0,hnd);
    c_cal_gui('update_SAVEbuttons')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % press_SaveALLbutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'press_SaveALLbutton'
    hnd = guidata(h0);
    
    f_name = {'./mEDSI.mat' './mCIS.mat' './mCIS.mat'};
    o_name = {'EFW' 'CISH' 'CISC'};
    v1 = {'Ddsi?' 'DHdsi?' 'DCdsi?'};
    v2 = {'Damp?' 'DHz?' 'DCz?'};
    
    for j=1:3
      c_eval([v1{j} '=hnd.' o_name{j} 'offset(?,1);' ...
        v2{j} '=hnd.' o_name{j} 'offset(?,2);'])
      c_eval(['hnd.' o_name{j} 'offset_save(?,:)=hnd.' o_name{j} 'offset(?,:);'])
      c_eval(['s?=[''' v1{j} ''' '' '' ''' v2{j} '''];']);
      
      var_s = [s1 ' ' s2 ' ' s3 ' ' s4];
      eval(['save ' f_name{j} ' ' var_s ' -mat -append'])
      irf_log('save',[var_s ' -> ' f_name{j}])
      eval(['clear ' var_s])
    end
    guidata(h0,hnd);
    c_cal_gui('update_SAVEbuttons')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % press_RESETbutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'press_RESETbutton'
    hnd = guidata(h0);
    
    k = D_findByName(hnd.Data,hnd.ActiveVar);
    if strcmp(hnd.Data{k}.type,'E')
      hnd.off = hnd.EFWoffset_save(hnd.Data{k}.cl_id,:);
    elseif strcmp(hnd.Data{k}.type,'V')
      if strcmp(hnd.Data{k}.sen,'HIA')
        hnd.off = hnd.CISHoffset_save(hnd.Data{k}.cl_id,:);
      else, hnd.off = hnd.CISCoffset_save(hnd.Data{k}.cl_id,:);
      end
    else, disp('Oj! Something is rong with data.type')
    end
    comp_s = 'XYZ';
    for comp=1:3
      switch comp
        case 1
          val = real(hnd.off(1));
        case 2
          val = imag(hnd.off(1));
        case 3
          val = hnd.off(2);
      end
      set(eval(['hnd.D' comp_s(comp) 'edit']),'String',num2str(val));
      set(eval(['hnd.D' comp_s(comp) 'slider']),'Value',val);
    end
    
    guidata(h0,hnd);
    c_cal_gui('update_off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update_Dslider
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_Dslider'
    hnd = guidata(h0);
    val = get(eval(['hnd.D' comp_s 'slider']),'Value');
    set(eval(['hnd.D' comp_s 'edit']),'String',num2str(val));
    switch comp
      case 1
        hnd.off(1) = val + 1i*imag(hnd.off(1));
      case 2
        hnd.off(1) = 1i*val + real(hnd.off(1));
      case 3
        hnd.off(2) = val;
    end
    %disp(sprintf('%s : offsets %f %f %f',action,real(hnd.off(1)),imag(hnd.off(1)),hnd.off(2)))
    guidata(h0,hnd);
    c_cal_gui('update_off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update_Dedit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_Dedit'
    hnd = guidata(h0);
    val = str2double(get(eval(['hnd.D' comp_s 'edit']),'String'));
    % Determine whether val is a number between 0 and 1
    if isnumeric(val) && length(val)==1 && ...
        val >= get(eval(['hnd.D' comp_s 'slider']),'Min') && ...
        val <= get(eval(['hnd.D' comp_s 'slider']),'Max')
      set(eval(['hnd.D' comp_s 'slider']),'Value',val);
      switch comp
        case 1
          hnd.off(1) = val + 1i*imag(hnd.off(1));
        case 2
          hnd.off(1) = 1i*val + real(hnd.off(1));
        case 3
          hnd.off(2) = val;
      end
      %disp(sprintf('%s : offsets %f %f %f',action,real(hnd.off(1)),imag(hnd.off(1)),hnd.off(2)))
      guidata(h0,hnd);
      c_cal_gui('update_off')
    else
      set(hnd.DXedit,'String','You have entered an invalid entry ');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update_Dcheckbox
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_Dcheckbox'
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
    
    % Check if we try to hide the active variable
    if ~isempty(D_findByName(hnd.Data(ii),hnd.ActiveVar))
      set(eval(['hnd.C' num2str(cl_id) 'checkbox']),'Value',1)
      disp([action ': cannot hide the active variable'])
      return
    end
    
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
      kk = zeros(size(ii));
      for j=1:length(ii)
        if hnd.Data{ii(j)}.visible == 0, kk(j) = j; end
      end
      kk(kk==0) = [];
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
      if inprog_mtx
        set(eval(['hnd.DATA' vs 'checkbox']),'Value',0)
        return
      end
      inprog_mtx = 1; %#ok<NASGU>
      
      % Plot the varible
      j = D_findByName(hnd.Data,vs);
      %disp(['plotting ' hnd.Data{j}.name])
      hnd.Data{j}.visible = 1;
      hnd.last = {vs};
      guidata(h0,hnd);
      
      inprog_mtx = 0;
      c_cal_gui('replot')
      
      hnd = guidata(h0);
      
      % Check if we need to show C# checkBox
      if get(eval(['hnd.C' num2str(hnd.Data{j}.cl_id) 'checkbox']),'Value')==0
        set(eval(['hnd.C' num2str(hnd.Data{j}.cl_id) 'checkbox']),...
          'Value',1)
      end
    else
      if inprog_mtx
        set(eval(['hnd.DATA' vs 'checkbox']),'Value',0)
        return
      end
      inprog_mtx = 1;
      
      % Check if we are hiding the active variable
      if strcmp(hnd.ActiveVar,vs)
        disp('you cannot hide the active variable')
        set(eval(['hnd.DATA' vs 'checkbox']),'Value',1)
        return
      end
      
      % Hide the varible
      j = D_findByName(hnd.Data,vs);
      hnd.Data{j}.visible = 0;
      
      % Hide the markers so thet they will not contribute YLim
      ts_tmp = hide_t_marker(hnd,hnd.ts_marker);
      te_tmp = hide_t_marker(hnd,hnd.te_marker);
      
      % Delete the variable
      delete(hnd.Data{j}.ploth)
      h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
      set(h,'YLimMode','auto')
      
      % Show the markers back
      if ~isempty(ts_tmp.t)
        hnd.ts_marker = replot_t_marker(hnd,ts_tmp);
      end
      if ~isempty(te_tmp.t)
        hnd.te_marker = replot_t_marker(hnd,te_tmp);
      end
      
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
      
      inprog_mtx = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update_DATAradiobutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'update_DATAradiobutton'
    hnd = guidata(h0);
    if get(eval(['hnd.DATA' vs 'radiobutton,']),'Value')==1
      vs_old = hnd.ActiveVar;
      
      j = D_findByName(hnd.Data,vs);
      % Show the new data if not visible
      if hnd.Data{j}.visible==0
        %disp('update_DATAradiobutton: replot')
        set(eval(['hnd.DATA' vs 'checkbox']),'Value',1)
        c_cal_gui(['update_DATA' vs 'checkbox'])
        hnd = guidata(h0);
      end
      
      j = D_findByName(hnd.Data,vs);
      % Give the new data 4*width and active color
      if hnd.Data{j}.visible
        for k=1:length(hnd.Data{j}.ploth)
          set(hnd.Data{j}.ploth(k),...
            'LineWidth',get(hnd.Data{j}.ploth(k),'LineWidth')*4,...
            'Color',active_p_color)
        end
        hnd.ActiveVar = vs;
        hnd.old_ActiveVar = vs_old;
      else
        return
      end
      
      % Hide the old active rb
      set(eval(['hnd.DATA' vs_old 'radiobutton,']),'Value',0)
      
      % Give the old data usual width and color
      j = D_findByName(hnd.Data,vs_old);
      for k=1:length(hnd.Data{j}.ploth)
        set(hnd.Data{j}.ploth(k),...
          'LineWidth',get(hnd.Data{j}.ploth(k),'LineWidth')/4,...
          'Color',hnd.Data{j}.plot_color)
      end
      
      guidata(h0,hnd);
      c_cal_gui('ch_active_var')
    else
      % Show it back
      eval(['set(hnd.DATA' vs 'radiobutton,''Value'',1)'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % click_axes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'click_axes'
    hnd = guidata(h0);
    
    replot_ts_marker = 0;
    replot_te_marker = 0;
    
    t = get(eval(['hnd.' curr_ax 'axes']),'CurrentPoint');
    t = t(1);
    
    % If we press the left mouse button we modify TS_MARKER
    % If we press the right mouse button of Control-click
    % we modify TE_MARKER
    % Otherwise re return
    if strcmp(get(h0,'SelectionType'),'normal')
      hnd.ts_marker.t = t;
      replot_ts_marker = 1;
      
      % We check the TS_MARKER should always precede TE_MARKER
      if ~isempty(hnd.te_marker)
        if hnd.ts_marker.t > hnd.te_marker.t
          hnd.ts_marker.t = hnd.te_marker.t;
        end
      end
    elseif strcmp(get(h0,'SelectionType'),'alt')
      hnd.te_marker.t = t;
      replot_te_marker = 1;
      
      % We check the TS_MARKER should always precede TE_MARKER
      if ~isempty(hnd.ts_marker)
        if hnd.ts_marker.t > hnd.te_marker.t
          hnd.te_marker.t = hnd.ts_marker.t;
        end
      end
    else, return
    end
    
    if replot_ts_marker
      hnd.ts_marker = replot_t_marker(hnd,hnd.ts_marker);
    elseif replot_te_marker
      hnd.te_marker = replot_t_marker(hnd,hnd.te_marker);
    end
    
    % Enable menu
    if ~isempty(hnd.ts_marker) && ~isempty(hnd.te_marker)
      set(hnd.menu_zoom_in,'Enable','on')
      set(hnd.menu_cut_int,'Enable','on')
    end
    
    guidata(h0,hnd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % press_RWDbutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'press_RWDbutton'
    hnd = guidata(h0);
    
    dt_tmp = hnd.tlim(end,2) - hnd.tlim(end,1);
    hnd.ts_marker.t = hnd.tlim(end,1) - dt_tmp;
    hnd.te_marker.t = hnd.tlim(end,2) - dt_tmp;
    clear dt_tmp
    
    guidata(h0,hnd);
    c_cal_gui('zoom_in')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % press_FWDbutton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'press_FWDbutton'
    hnd = guidata(h0);
    
    dt_tmp = hnd.tlim(end,2) - hnd.tlim(end,1);
    hnd.ts_marker.t = hnd.tlim(end,1) + dt_tmp;
    hnd.te_marker.t = hnd.tlim(end,2) + dt_tmp;
    clear dt_tmp
    
    guidata(h0,hnd);
    c_cal_gui('zoom_in')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % zoom_in
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'zoom_in'
    hnd = guidata(h0);
    
    % Check if we really have update tlim
    % This check is probably unnecessary
    if hnd.ts_marker.t==hnd.tlim(end,1) && hnd.te_marker.t==hnd.tlim(end,2)
      return
    end
    
    hnd.tlim(end+1,:) = [hnd.ts_marker.t hnd.te_marker.t];
    hide_t_marker(hnd,hnd.ts_marker); hnd.ts_marker = [];
    hide_t_marker(hnd,hnd.te_marker); hnd.te_marker = [];
    
    h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
    irf_zoom(h,'x',hnd.tlim(end,:));
    set(h,'YLimMode','auto')
    irf_timeaxis(h)
    
    % Enable/disable menus
    set(hnd.menu_zoom_in,'Enable','off')
    set(hnd.menu_cut_int,'Enable','off')
    set(hnd.menu_zoom_un,'Enable','on')
    set(hnd.menu_zoom_rs,'Enable','on')
    
    guidata(h0,hnd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % zoom_un
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'zoom_un'
    hnd = guidata(h0);
    
    % Check if we really can undo
    % This check is probably unnecessary
    if size(hnd.tlim,1)<=1
      disp('Cannot undo. Already at the the original time span.')
      return
    end
    
    hnd.tlim(end,:) = [];
    hide_t_marker(hnd,hnd.ts_marker); hnd.ts_marker = [];
    hide_t_marker(hnd,hnd.te_marker); hnd.te_marker = [];
    
    h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
    irf_zoom(h,'x',hnd.tlim(end,:));
    set(h,'YLimMode','auto')
    irf_timeaxis(h)
    
    % Enable/disable menus
    set(hnd.menu_zoom_in,'Enable','off')
    set(hnd.menu_cut_int,'Enable','off')
    if size(hnd.tlim,1)==1
      set(hnd.menu_zoom_un,'Enable','off')
      set(hnd.menu_zoom_rs,'Enable','off')
    end
    
    guidata(h0,hnd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % zoom_un
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'zoom_rs'
    hnd = guidata(h0);
    
    % Check if we really can undo
    % This check is probably unnecessary
    if size(hnd.tlim,1)<=1
      disp('Cannot reset. Already at the the original time span.')
      return
    end
    
    hnd.tlim(2:end,:) = [];
    hide_t_marker(hnd,hnd.ts_marker); hnd.ts_marker = [];
    hide_t_marker(hnd,hnd.te_marker); hnd.te_marker = [];
    
    h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];
    irf_zoom(h,'x',hnd.tlim(end,:));
    set(h,'YLimMode','auto')
    irf_timeaxis(h)
    
    % Enable/disable menus
    set(hnd.menu_zoom_in,'Enable','off')
    set(hnd.menu_cut_int,'Enable','off')
    set(hnd.menu_zoom_un,'Enable','off')
    set(hnd.menu_zoom_rs,'Enable','off')
    
    guidata(h0,hnd);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % show_raw
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'show_raw'
    hnd = guidata(h0);
    
    j = D_findByName(hnd.Data,hnd.ActiveVar);
    
    % Create figure
    if find(get(0,'children')==raw_fig_id)
      pos_old = get(raw_fig_id,'Position');
    else, pos_old = [];
    end
    fig = figure(raw_fig_id);
    clf
    set(raw_fig_id,'Name', 'RAW DATA')
    if isempty(pos_old), set(fig,'Position', pos_raw_fig(hnd.scrn_size,:)), end
    
    if strcmp(hnd.Data{j}.type,'E')
      d_tmp = {};
      leg_tmp = {};
      dp_tmp = {};
      dp_sid = [];
      legp_tmp = {};
      % Load E
      if strcmp(hnd.Data{j}.sen,'1234')
        s = {'12','32','34'};
        probes = 1:4;
      else
        s = {num2str(hnd.Data{j}.sen)};
        sss = num2str(hnd.Data{j}.sen);
        probes = [str2double(sss(1)) str2double(sss(2))];
        clear sss
      end
      for sid = 1:length(s)
        [ok,E_tmp] = c_load(['wE?p' s{sid}],hnd.Data{j}.cl_id);
        if ok
          E_tmp = irf_tlim(E_tmp,hnd.tlim(end,:));
          d_tmp = [d_tmp {E_tmp}]; clear E_tmp %#ok<AGROW>
          leg_tmp = [leg_tmp {irf_ssub(['wE?p' s{sid}],hnd.Data{j}.cl_id)}]; %#ok<AGROW>
        end
      end
      % Load P
      for sid = probes
        [ok,P_tmp] = c_load(['P10Hz?p' num2str(sid)],hnd.Data{j}.cl_id);
        if ok
          P_tmp = irf_tlim(P_tmp,hnd.tlim(end,:));
          dp_tmp = [dp_tmp {P_tmp}]; clear P_tmp %#ok<AGROW>
          dp_sid = [dp_sid sid]; %#ok<AGROW>
          legp_tmp = [legp_tmp {irf_ssub(['P10Hz?p' num2str(sid)],hnd.Data{j}.cl_id)}]; %#ok<AGROW>
        end
      end
      if isempty(d_tmp) && isempty(dp_tmp), return, end
      figure(raw_fig_id), clf
      if ~isempty(d_tmp)
        h1 = irf_subplot(2,1,1);
        irf_plot(d_tmp,'comp')
        irf_zoom(h1,'x',hnd.tlim(end,:));
        set(h1,'YLimMode','auto')
        ylabel('E [mV/m]')
        legend(leg_tmp{:},'Location','NorthEastOutside')
      end
      if ~isempty(dp_tmp)
        h2 = irf_subplot(2,1,2);
        irf_plot(dp_tmp,'comp')
        irf_zoom(h1,'x',hnd.tlim(end,:));
        set(h2,'YLimMode','auto')
        ylabel(h2,'V [V]')
        legend(h2,legp_tmp{:},'Location','NorthEastOutside')
      end
      irf_plot_axis_align;
      clear d_tmp dp_tmp leg_tmp legp_tmp s sid dp_sid ok
    else % Display V CIS
      if strcmp(hnd.Data{j}.sen,'HIA'), v_s = irf_ssub('VCh?',hnd.Data{j}.cl_id);
      else, v_s = irf_ssub('VCp?',hnd.Data{j}.cl_id);
      end
      [ok,V_tmp] = c_load(v_s);
      if ~ok, return, end
      V_tmp = irf_tlim(V_tmp,hnd.tlim(end,:));
      if ~isempty(V_tmp)
        figure(raw_fig_id), clf
        irf_plot(V_tmp)
        ylabel([v_s ' GSE [km/s]'])
        irf_zoom(gca,'x',hnd.tlim(end,:));
        set(gca,'YLimMode','auto')
        legend('Vx','Vy','Vz')
      end
      clear V_tmp ok v_s
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % show_b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'show_b'
    hnd = guidata(h0);
    
    j = D_findByName(hnd.Data,hnd.ActiveVar);
    
    % Create figure
    if find(get(0,'children')==raw_fig_id)
      pos_old = get(raw_fig_id,'Position');
    else, pos_old = [];
    end
    fig = figure(raw_fig_id);
    clf
    set(raw_fig_id,'Name', 'Magnetic field')
    if isempty(pos_old), set(fig,'Position', pos_raw_fig(hnd.scrn_size,:)), end
    h1 = irf_plot(irf_abs(hnd.Data{j}.B));
    irf_zoom(h1,'x',hnd.tlim(end,:));
    ylabel(['B SC' num2str(hnd.Data{j}.cl_id) ' [nT]'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % show_spect
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'show_spect'
    hnd = guidata(h0);
    
    j = D_findByName(hnd.Data,hnd.ActiveVar);
    
    if ~strcmp(hnd.Data{j}.type,'E')
      disp('should not see this: no need for cis')
      return,
    end
    
    % Create figure
    if find(get(0,'children')==spect_fig_id)
      pos_old = get(spect_fig_id,'Position');
    else, pos_old = [];
    end
    fig = figure(spect_fig_id);
    clf
    set(spect_fig_id,'Name', 'SPECTRUM')
    if isempty(pos_old), set(fig,'Position', pos_spect_fig(hnd.scrn_size,:)), end
    
    if strcmp(hnd.Data{j}.sen,'1234'), E_tmp = hnd.Data{j};
    else
      ii = D_findByName(hnd.Data,irf_ssub('diE?p1234',hnd.Data{j}.cl_id));
      if ~isempty(ii)
        E_tmp = hnd.Data{ii};
      else
        c_log('proc','cannot make spectra: no full res E loaded')
        return
      end
    end
    
    ndata = length(E_tmp.data(:,1));
    % Guess the sampling frequency
    sf = c_efw_fsample(E_tmp.data(:,1));
    if sf==25, nfft = 512;
    else, nfft = 4096;
    end
    
    % Exclude NaNs
    ii = find(~isnan(E_tmp.data(:,2)));
    
    % Check if we have enough data
    if length(ii) < nfft/2
      irf_log('proc','not enough points for FFT')
      return
    elseif length(ii) < nfft, nfft = nfft/2;
    end
    
    % Compute FFT
    [Px,Freqx] = irf_psd(E_tmp.data(ii,[1 2]),nfft,sf);
    [Py,Freqy] = irf_psd(E_tmp.data(ii,[1 3]),nfft,sf);
    
    figure(spect_fig_id), clf
    loglog(Freqx,Px,Freqy,Py)
    set(gca,'Xlim',[.05 5]);
    set(gca,'YLimMode','auto')
    freqs = [.25 .5 1]; % harmonics of spin freq
    dy = get(gca,'YLim');
    hold on
    for j=1:length(freqs)
      plot(freqs(j)*[1 1],dy,'r')
    end
    hold off
    ylabel('PSD [(mV/m)^2/Hz]')
    title(...
      sprintf('Spectrum of %s %d data points, nfft=%d',E_tmp.name,ndata,nfft))
    legend('X DSI','Y DSI')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cut_int
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'cut_int'
    hnd = guidata(h0);
    
    % Check if we really have update tlim
    % This check is probably unnecessary
    if hnd.ts_marker.t==hnd.tlim(end,1) && hnd.te_marker.t==hnd.tlim(end,2)
      return
    end
    
    tlim_tmp = [hnd.ts_marker.t hnd.te_marker.t];
    j = D_findByName(hnd.Data,hnd.ActiveVar);
    
    btn = questdlg([hnd.ActiveVar ' : Remove ' num2str(tlim_tmp(2)-tlim_tmp(1)) ...
      ' sec starting from ' epoch2iso(tlim_tmp(1),1)], ...
      'Confirmation', ...
      'Yes','No','No');
    if strcmp(btn,'No'), return, end
    
    if strcmp(hnd.Data{j}.type,'E')
      btn = questdlg('Remove corresponding single probe potentials?', ...
        'Question', ...
        'Yes','No','No');
      rm_pot = 1;
      if strcmp(btn,'No'), rm_pot = 0; end
      if strcmp(hnd.Data{j}.sen,'1234')
        var_list = {'wE?p12','wE?p32','wE?p34'};
        if rm_pot, var_list = [var_list, ...
            {'P10Hz?p1','P10Hz?p2','P10Hz?p3','P10Hz?p4'}];
        end
      else
        var_list = {['wE?p' hnd.Data{j}.sen]};
        if rm_pot, var_list = [var_list,...
            {['P10Hz?p' hnd.Data{j}.sen(1)],['P10Hz?p' hnd.Data{j}.sen(2)]}];
        end
      end
      
      % Remove the desired interval from the raw data
      n_ok = 0;
      tint_tmp = [];
      for vi = 1:length(var_list)
        v_s = irf_ssub(var_list{vi},hnd.Data{j}.cl_id);
        ok = c_load(v_s);
        if ok
          irf_log('proc',...
            ['cutting ' v_s ' ' num2str(tlim_tmp(2)-tlim_tmp(1)) 'sec from ' epoch2iso(tlim_tmp(1),1)])
          
          eval(['ii=find(' v_s '(:,1)>tlim_tmp(1) & ' v_s '(:,1)<tlim_tmp(2));']);
          dsc = c_desc(v_s);
          eval([v_s '(ii,:)=[];save ' dsc.file ' ' v_s ' -append']);
          if isempty(tint_tmp), eval(['tint_tmp=[' v_s '(1,1) '  v_s '(end,1)];']);
          else
            eval(['tint1_tmp=[' v_s '(1,1) '  v_s '(end,1)];']);
            if tint1_tmp(1)<tint_tmp(1), tint_tmp(1) = tint1_tmp(1); end
            if tint1_tmp(2)>tint_tmp(2), tint_tmp(2) = tint1_tmp(2); end
            clear tint1_tmp
          end
          n_ok = n_ok + 1;
          eval(['clear ' v_s]);
        end
      end
      clear ii v_s
      if n_ok<1, disp('was no raw data'), return, end
      
      % Reprocess the data using the updated raw data
      c_get_batch(tint_tmp(1),tint_tmp(2)-tint_tmp(1),hnd.Data{j}.cl_id,...
        'vars','e','nosrc');
      if rm_pot
        c_get_batch(tint_tmp(1),tint_tmp(2)-tint_tmp(1),hnd.Data{j}.cl_id,...
          'vars','p','nosrc');
      end
      % Variables which we need to reload
      var_list = {'diEs?p12','diEs?p32','diEs?p34','diE?p1234'};
      if rm_pot, var_list = [var_list, {'P?'}]; end
      
    else % Cut V CIS
      if strcmp(hnd.Data{j}.sen,'HIA'), v_s = irf_ssub('VCh?',hnd.Data{j}.cl_id);
      else, v_s = irf_ssub('VCp?',hnd.Data{j}.cl_id);
      end
      % Variables which we need to reload
      var_list = {['di' v_s],v_s}; n_ok = 0;
      for j=1:length(var_list)
        v_s = var_list{j};
        ok = c_load(v_s);
        if ok
          irf_log('proc',...
            ['cutting ' v_s ' ' num2str(tlim_tmp(2)-tlim_tmp(1)) 'sec from ' epoch2iso(tlim_tmp(1),1)])
          eval(['ii=find(' v_s '(:,1)>tlim_tmp(1) & ' v_s '(:,1)<tlim_tmp(2));']);
          dsc = c_desc(v_s);
          eval([v_s '(ii,:)=[];save ' dsc.file ' ' v_s ' -append; clear v_s']);
          n_ok = n_ok + 1;
        end
      end
      clear ii v_s
      if n_ok<1, disp('was no raw data'), return, end
      
      % We need only diV*
      var_list = var_list(1);
    end
    
    % Reload updated variables
    n_ok = 0;
    for vi = 1:length(var_list)
      v_s = irf_ssub(var_list{vi},hnd.Data{j}.cl_id);
      jj = D_findByName(hnd.Data,v_s);
      if ~isempty(jj)
        ok = c_load(v_s);
        if ok
          eval(['hnd.Data{jj}.data=' v_s ';clear ' v_s])
          n_ok = n_ok + 1;
        else
          irf_log('load',...
            ['Cannot load ' v_s '. Removing it from the list.']);
          hnd.Data{jj} = [];
        end
      end
    end
    if n_ok<1
      disp('could not load any updated data. this is strange')
      return
    end
    
    guidata(h0,hnd);
    c_cal_gui('replot_all')
  otherwise
    disp('wrong action')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                        H E L P   F U N C T I O N S
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function replot_t_marker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmarker = hide_t_marker(hnd,marker) %#ok<INUSL>
newmarker.t = [];
if ~isempty(marker)
  newmarker.t = marker.t;
  % Hide the marker if marker.h exists
  if isfield(marker,'h')
    try
      for j=4:-1:1, delete(marker.h(j)),end
    catch err
      if ~strcmp(err.identifier,'MATLAB:hg:udd_interface:CannotDelete')
        rethrow(err)
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function replot_t_marker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmarker = replot_t_marker(hnd,marker)
h = [hnd.Xaxes hnd.Yaxes hnd.Zaxes hnd.AUXaxes];

% Try to hide the marker
newmarker = hide_t_marker(hnd,marker);

for j=1:length(h)
  yy = get(h(j),'YLim');
  hold(h(j),'on')
  newmarker.h(j) = plot(h(j),[newmarker.t newmarker.t],yy,...
    'LineWidth',1,'Color','red');
  hold(h(j),'off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function corr_v_velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = corr_v_velocity(v,offset)
out = v;
out(:,2) = v(:,2) - real(offset(1));
out(:,3) = v(:,3) - imag(offset(1));
out(:,4) = v(:,4) - offset(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function p_style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out, color] = p_style(cl_id, inst, sen)
out = '';
colrs = [0 0 0; 1 0 0; 0 .5 0; 0 0 1]; % black, red, dark green, blue
switch inst
  case 'EFW'
    switch sen
      case '1234'
        out = '-';
      case '12'
        out = '--';
      case '32'
        out = '--';
      case '34'
        out = '-.';
      otherwise
        disp('unknown sensor')
    end
  case 'EDI'
    out = '.';
  case 'CIS'
    switch sen
      case 'HIA'
        out = '.-';
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
  %disp(['get_plot_data : ' data.name])
  
  % Correct offsets
  if data.editable
    %disp(['get_plot_data : correcting offsets in ' data.name])
    switch data.type
      case 'E'
        ofs = hnd.EFWoffset(data.cl_id,:);
        %disp(sprintf('offsets are: %f %f %f',real(ofs(1)),imag(ofs(1)),ofs(2)))
        p_data = caa_corof_dsi(p_data,...
          real(ofs(1)),imag(ofs(1)),ofs(2));
      case 'V'
        if strcmp(data.sen,'HIA')
          ofs = hnd.CISHoffset(data.cl_id,:);
        else, ofs = hnd.CISCoffset(data.cl_id,:);
        end
        %disp(sprintf('offsets are: %f %f %f',real(ofs(1)),imag(ofs(1)),ofs(2)))
        if isempty(data.B), p_data = [];
        else
          p_data = corr_v_velocity(p_data,ofs);
          [xxx,p_data]=irf_dec_parperp(data.B,p_data); %#ok<ASGLU>
          clear xxx
        end
      otherwise
        irf_log('proc','Unknown data type.')
    end
  end
  
  % Calculate Ez, convert V->E and V->E
  switch data.type
    case 'E'
      if any(p_data(:,4))==0
        if isempty(data.B)
          irf_log('proc','B is empty. No Ez')
        else
          warning('off','MATLAB:interp1:NaNinY')
          [p_data,angle] = irf_edb(p_data,...
            data.B,hnd.ang_limit);
          warning('on','MATLAB:interp1:NaNinY')
          ii = find(abs(angle) < hnd.ang_limit);
          if length(ii) > 1
            p_data(ii,4) = p_data(ii,4)*NaN;
          end
        end
      end
      % E->V if in V mode
      if ~hnd.mode
        if ~isempty(data.B)
          warning('off','MATLAB:interp1:NaNinY')
          p_data = irf_e_vxb(p_data,data.B,-1);
          warning('on','MATLAB:interp1:NaNinY')
        else
          data.visible = 0;
          p_data = [];
          irf_log('proc','B is empty. No ExB')
        end
      end
    case 'V'
      % V->E if in E mode
      if hnd.mode
        if ~isempty(data.B)
          p_data = irf_tappl(...
            irf_cross(p_data,data.B),'*(-1e-3)');
        else
          data.visible = 0;
          p_data = [];
          irf_log('proc','B is empty. No VxB')
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
  if ~isempty(data_rec), data_ii = [data_ii data_rec]; end %#ok<AGROW>
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
    data_ii = [data_ii j]; %#ok<AGROW>
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
if isempty(list), return, end
if ischar(s_list)
  % fast search
  for j=1:length(list)
    if strcmp(list{j},s_list), ii = j; return, end
  end
else
  for k=1:length(s_list)
    for j=1:length(list)
      if strcmp(list{j},s_list{k}), ii = [ii j]; break, end %#ok<AGROW>
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L_rm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newlist = L_rm(list,s_list)
newlist = list;
if ischar(s_list)
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
