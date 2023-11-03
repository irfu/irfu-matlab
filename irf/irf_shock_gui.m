function [out] = irf_shock_gui(scd,varName)
%IRF_SHOCK_GUI GUI for determining shock parameters.
%
%   THIS FUNCTION IS IN DEVELOPMENT AND HAS NOT BEEN PROPERLY TESTED!
%
%   IRF_SHOCK_GUI(scd) Starts a GUI. scd is a struct with fields containing
%   TSeries objects of spacecraft data. Select up- and downstream intervals
%   of the shock by clicking in the plot. Optionally select the shock foot.
%   Then select methods to calculate shock normal and shock speed. The
%   methods are described in irf_shock_normal. Click "Calculate" to display
%   the results.
%
%   scd field names:
%       B       -   3D magnetic field (nT)
%       V       -   3D ion or electron bulk velocity (km/s)
%       n       -   ion or electron number density (cm^-3)
%       Ti/Te   -   ion/electron temperature (eV) (optional)
%       R       -   Spacecraft position as given by R =
%                   mms.get_data('R_gse',tint) (optional)
%
%   Results displayed:
%       n       -   Shock normal vector
%       Vsh     -   Shock speed along normal vector
%       thBn    -   Angle between normal and upstream magnetic field
%       thVn    -   Angle between normal and upstream flow
%       nd/nu   -   Density compression rate
%       Bd/Bu   -   Magnetic compression rate, norm(Bd)/norm(Bu)
%       Ma      -   Alfven Mach #, sc frame or NIF
%       Mf      -   Fast Magnetosonic Mach #, sc frame or NIF
%       Ms      -   Sonic Mach #, sc frame or NIF
%       beta_i  -   Upstream ion beta
%       beta_e  -   Upstream electron beta
%
%   IRF_SHOCK_GUI(scd,varName) Also saves a variable, containing normal
%   vectors and plasma parameters, to the Workspace. Variable contains:
%       nvec    -   As returned by irf_shock_normal
%       par     -   As returned by irf_shock_parameters
%       data    -   Up- and downstream values used
%
%   See also:
%       IRF_SHOCK_NORMAL, IRF_SHOCK_PARAMETERS, IRF_4_V_GUI, IRF_MINVAR_GUI
%

%   Written by: Andreas Johlander, andreasj@irfu.se
%
%   TODO: Fix all velocity methods
%       Replace uicontrols with text objects

%   Update 1:
%       By Ahmad Lalti,  on 5-2-2021.
%       Update descrition:
%       - Give the option for manually inputing the
%         upstream and downstream parameters.
%       - Display the magnetic field and velocity as vectors and not as
%         norms
%        01-07-2021
%       - Use nanmean instead of mean for omnidata averaging
%       - removed correction for Vy abberation in omni data since its
%       already done from SDC
%       - use tlim instead of resample to cut the time interval for omni
%       data averaging



%% Structure of ud:
%   Nin             -   Number of data inputs
%   uih             -   Handles in panels:
%       up  - Upstream panel
%           panel - The panel
%           pb - Push-button
%       dw  - Downstream panel
%           panel - The panel
%           pb - Push-button
%       mt  - Method panel
%           panel - The panel
%           ntx - Normal text
%           npu - Normal pop-up menu
%           vtx - Velocity text
%           vpu - Normal pop-up menu
%       cl  - Shock parameter panel
%           panel - The panel
%           pb - Push-button
%           nvec
%           Vsh
%           thBn
%           thVn
%           r_n
%           r_b
%           Ma
%           Mms
%           Ms
%           beta_i
%           beta_e
%
%   ax              -   Axis handles, 1xNin array
%   params          -   Parameters to plot (maybe more)
%           Bd
%           Bu
%           Vd
%           Vu
%           nd
%           nu
%   scd             -   Structure containing data, can also be string for
%                       internal use
%   tu              -   Upstream time, interval 1x2 array
%   td              -   Downstream time, interval 1x2 array
%   varName         -   Name of variable in workspace
%   shp             -   Dunno
%   normal_method   -   Method for normal vector
%   vel_method      -   Method for shock velocity
%   mach_method     -   Method for Mach numbers
%   use_omni        -   Boolean for if omni data is used as upstream
%
%   %Others
%   t_start_epoch
%   zoomStack

%% handle input
if ischar(scd)
  % not first call
  ud = get(gcf,'userdata');
  % switch for action
  switch scd
    case 'clu' % click upstream
      ud = clickt(ud,{'u'});
      ud = mark_times(ud);
      ud = get_avg_field(ud,ud.scd,{'u'});
      % Initiate ugly fields to replace omni if needed
      fn = fieldnames(ud.scd);
      for k = 1:length(fn)
        ud.sc_up.([fn{k},'u']) = ud.params.([fn{k},'u']);
      end
      ud = set_omni(ud); % to ensure omni is used if chosen
      ud = display_vals(ud);
      set(gcf,'userdata',ud)
    case 'cld' % click downstream
      ud = clickt(ud,{'d'});
      ud = get_avg_field(ud,ud.scd,{'d'});
      ud = mark_times(ud);
      ud = display_vals(ud);
      set(gcf,'userdata',ud)
    case 'clf' % click foot
      ud = clickt(ud,{'f'});
      ud = get_avg_field(ud,ud.scd,{'f'});
      ud = mark_times(ud);
      ud = display_vals(ud);
      set(gcf,'userdata',ud)
    case 'set_met' % click calculate
      ud = set_methods(ud);
      set(gcf,'userdata',ud)
    case 'set_omni'
      ud = set_omni(ud);
      ud = display_vals(ud);
      set(gcf,'userdata',ud)
    case 'plot_omni'
      plot_omni(ud)
    case 'manual_input'
      ud = manual_input(ud);
      set(gcf,'userdata',ud)
    case 'calc' % click calculate
      % first update methods
      ud = set_methods(ud);

      % time is set between up- and downstream intervals
      ud.params.t = irf_time(mean([ud.tu(2),ud.td(1)]),'epoch>epochtt');
      % return up and downstream tints for output
      ud.params.tintu = irf_time(ud.tu,'epoch>epochtt');
      ud.params.tintd = irf_time(ud.td,'epoch>epochtt');
      ud.params.tintf = irf_time(ud.tf,'epoch>epochtt');

      % check for manual input to overwrite other inputs
      ud = manual_input(ud);
      % get shock parameters (Mach #, beta, Fcp,...)

      ud.shp.par = irf_shock_parameters(ud.params);

      % set parameters for shock foot width methods
      if isfield(ud.shp.par,'Fcpf')
        ud.params.Fcp = ud.shp.par.Fcpf; % foot ion cycl. freq.
        ud.params.dTf = diff(ud.tf);
        ud.params.d2u = sign(mean(ud.tu)-mean(ud.td));
      end

      % calculate shock normals and speeds
      ud.shp.nvec = irf_shock_normal(ud.params);

      ud.shp.data = ud.params;

      ud = display_prop(ud);



      set(gcf,'userdata',ud)
  end
else

  % input names
  fn = fieldnames(scd);
  % sc position is moved
  if ismember('R',fn)
    R = scd.R;
    scd = rmfield(scd,'R');
    fn = fieldnames(scd);
    inpR = 1;
  else
    inpR = 0;
  end
  % number of data inputs
  Nin = numel(fn);

  % possible inputs
  poss_inp = {'B','V','n','Ti','Te','R'};
  for k = 1:Nin % remove all fields that are not allowed
    if ~ismember(fn{k},poss_inp)
      irf.log('w',['Removes field ',fn{k},'.']);
      scd = rmfield(scd,fn{k});
      Nin = Nin-1;
    end
  end

  % User data that is used everywhere
  ud = [];
  % number of inputs
  ud.Nin = Nin;
  % spacecraft data
  ud.scd = scd;
  % handles to UI elements
  ud.uih = [];

  if nargin == 1 % do not write to workspace if no filename given
    ud.doSave = 0;
  else % if filename is given, save to workspace
    ud.doSave = 1;
    ud.varName = varName;
  end
  % parameters used in calculation
  ud.params = [];
  if inpR
    ud.params.R = R;
  end
  % initiate parameters
  ud.params.Bu = NaN*ones(1,3); ud.params.Bd = NaN*ones(1,3); ud.params.Bf = NaN*ones(1,3);
  ud.params.nu = NaN; ud.params.nd = NaN; ud.params.nf = NaN;
  ud.params.Vu = NaN*ones(1,3); ud.params.Vd = NaN*ones(1,3); ud.params.Vf = NaN*ones(1,3);
  ud.params.Tiu = NaN; ud.params.Tid = NaN; ud.params.Tif = NaN;
  ud.params.Teu = NaN; ud.params.Ted = NaN; ud.params.Tef = NaN;
  % to make omni stuff work
  ud.sc_up.Bu = ud.params.Bu;
  ud.sc_up.nu = ud.params.nu;
  ud.sc_up.Vu = ud.params.Vu;
  ud.sc_up.Tiu = ud.params.Tiu;

  % initiate GUI
  ud = init_gui(ud);
  % set default to not use OMNI, must change in GUI as well if changed
  ud.use_omni.B = 0;
  ud.use_omni.n = 0;
  ud.use_omni.V = 0;
  ud.use_omni.Ti = 0;
  % plot data in panels
  ud = get_avg_field(ud,scd,[]);
  % align time axis
  fnp = fieldnames(scd);
  tint = scd.(fnp{1}).time([1,end]); % not optimal
  irf_zoom(ud.ax,'x',tint)
  % fix labels
  ud = set_labels(ud);
  % t_start_epoch
  gfud = get(gcf,'userdata');
  ud.t_start_epoch = gfud.t_start_epoch;
  % up/downstream time intervals
  ud.tu = ud.scd.B.time([1,end]).epochUnix;
  ud.td = [ud.tu(1);ud.tu(1)]; % so no interval is shown first
  ud.tf = ud.td;
  % normal vector struct
  ud.shp.nvec = [];
  % parameter struct
  ud.shp.par = [];
  % default normal method (if changed, must also change in pop-up menu)
  ud.normal_method = 'mx3';
  % default velocity method (if changed, must also change in pop-up menu)
  ud.vel_method = 'mf';
  % set figure userdata
  set(gcf,'userdata',ud)
end

if ud.doSave % if filename is given, save to workspace
  assignin('base',ud.varName, ud.shp)
end

if nargout == 1 % output
  out = ud.shp;
end

end


function [ud] = init_gui(ud)

% initiate figure
ax = irf_plot(ud.Nin,'newfigure');

ax(1).Position(3) = 0.45;
for k = 1:ud.Nin
  ax(k).Position(1) = 0.1;
end
pause(0.001)
irf_plot_axis_align(ax);
% axis handle array
ud.ax = ax;

%% Upstream panel
% panel
ud.uih.up.panel = uipanel('Units', 'normalized',...
  'position',[0.56 0.91 0.1467 0.075],...
  'fontsize',14,...
  'Title','Upstream');
% push button for time
ud.uih.up.pb = uicontrol('style','push',...
  'Units', 'normalized',...
  'Parent',ud.uih.up.panel,...
  'position',[0.05 0.2 0.8 0.6],...
  'fontsize',14,...
  'string','Set times',...
  'callback','irf_shock_gui(''clu'')');

%% Downstream panel
% panel
ud.uih.dw.panel = uipanel('Units', 'normalized',...
  'position',[0.7067 0.91 0.1467 0.075],...
  'fontsize',14,...
  'Title','Downstream');
% push button for time
ud.uih.dw.pb = uicontrol('style','push',...
  'Units', 'normalized',...
  'Parent',ud.uih.dw.panel,...
  'position',[0.05 0.2 0.8 0.6],...
  'fontsize',14,...
  'string','Set times',...
  'callback','irf_shock_gui(''cld'')');

%% Foot panel
% panel
ud.uih.ft.panel = uipanel('Units', 'normalized',...
  'position',[0.8534 0.91 0.1467 0.075],...
  'fontsize',14,...
  'Title','Shock foot');
% push button for time
ud.uih.ft.pb = uicontrol('style','push',...
  'Units', 'normalized',...
  'Parent',ud.uih.ft.panel,...
  'position',[0.05 0.2 0.8 0.6],...
  'fontsize',14,...
  'string','Set times',...
  'callback','irf_shock_gui(''clf'')');

%% OMNI panel
% panel
ud.uih.omni.panel = uipanel('Units', 'normalized',...
  'position',[0.56 0.8 0.44 0.11],...
  'fontsize',14,...
  'Title','Use OMNI as upstream');

% B-field checkbox
ud.uih.omni.pbB = uicontrol('style','checkbox',...
  'Units', 'normalized',...
  'Parent',ud.uih.omni.panel,...
  'position',[0.05 0.7 0.8 0.3],...
  'fontsize',14,...
  'string','B',...
  'callback','irf_shock_gui(''set_omni'')');
% density checkbox
ud.uih.omni.pbN = uicontrol('style','checkbox',...
  'Units', 'normalized',...
  'Parent',ud.uih.omni.panel,...
  'position',[0.25 0.7 0.8 0.3],...
  'fontsize',14,...
  'string','n',...
  'callback','irf_shock_gui(''set_omni'')');
% velocity checkbox
ud.uih.omni.pbV = uicontrol('style','checkbox',...
  'Units', 'normalized',...
  'Parent',ud.uih.omni.panel,...
  'position',[0.50 0.7 0.8 0.3],...
  'fontsize',14,...
  'string','V',...
  'callback','irf_shock_gui(''set_omni'')');
% temperature checkbox
ud.uih.omni.pbTi = uicontrol('style','checkbox',...
  'Units', 'normalized',...
  'Parent',ud.uih.omni.panel,...
  'position',[0.75 0.7 0.8 0.3],...
  'fontsize',14,...
  'string','Ti',...
  'callback','irf_shock_gui(''set_omni'')');

% push button to plot omni data
ud.uih.dw.pb = uicontrol('style','push',...
  'Units', 'normalized',...
  'Parent',ud.uih.omni.panel,...
  'position',[0.05 0.15 0.8 0.5],...
  'fontsize',14,...
  'string','Plot OMNI data',...
  'callback','irf_shock_gui(''plot_omni'')');

%% Method panel
% panel
ud.uih.mt.panel = uipanel('Units', 'normalized',...
  'position',[0.56 0.64 0.44 0.16],...
  'fontsize',14,...
  'Title','Methods');

% normal text
ud.uih.mt.ntx = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.uih.mt.panel,...
  'position',[0.05 0.7 0.3 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Normal: ');
% pop-up menu for normal method
ud.uih.mt.npu = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.uih.mt.panel,...
  'position',[0.3 0.7 0.65 0.2],...
  'fontsize',14,...
  'string',{'Magnetic coplanarity','Velocity coplanarity',...
  'Mixed mode 1','Mixed mode 2','Mixed mode 3',...
  '-------','Farris et al.','Slavin & Holzer mean','Peredo et al.',...
  'Fairfield Meridian 4o','Fairfield Meridian No 4o',...
  'Formisano Unnorm. z = 0'},...
  'Value',5,...
  'Callback','irf_shock_gui(''set_met'')');

% velocity text
ud.uih.mt.vtx = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.uih.mt.panel,...
  'position',[0.05 0.4 0.6 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Velocity: ');
% pop-up menu for velocity method
ud.uih.mt.vpu = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.uih.mt.panel,...
  'position',[0.3 0.4 0.65 0.2],...
  'fontsize',14,...
  'string',{'Gosling & Thomsen','Mass flux','Smith & Burton','Moses et al. (Experimental)'},...
  'Value',2,...
  'Callback','irf_shock_gui(''set_met'')');

% velocity text
ud.uih.mt.mtx = uicontrol('style','text',...
  'Units', 'normalized',...
  'Parent',ud.uih.mt.panel,...
  'position',[0.05 0.1 0.6 0.2],...
  'fontsize',14,...
  'HorizontalAlignment','left',...
  'string','Mach #: ');
% pop-up menu for velocity method
ud.uih.mt.mpu = uicontrol('style','popupmenu',...
  'Units', 'normalized',...
  'Parent',ud.uih.mt.panel,...
  'position',[0.3 0.1 0.65 0.2],...
  'fontsize',14,...
  'string',{'Spacecraft','NIF','NIF, Vsh = 0'},...
  'Value',1,...
  'Callback','irf_shock_gui(''set_met'')');

%% Parameter values panel
% panel
ud.uih.par.panel = uipanel('Units', 'normalized',...
  'position',[0.56 0.37 0.44 0.27],...
  'fontsize',14,...
  'Title','Up/downstream values');

% parameters to show
par_name = {'Bu','nu','Vu','Teu','Tiu','Bd','nd','Vd','Ted','Tid'};
par_str = par_name;

% positions of text boxes
nt = length(par_name);
post = zeros(nt,4);
% width
post(:,3) = 0.4;
% height
post(:,4) = 0.15;
% first column left
post(1:ceil(nt/2),1) = 0.15;
% second column left
post(ceil(nt/2)+1:end,1) = 0.65;
% first column up
post(1:ceil(nt/2),2) = fliplr(linspace(0.05,0.8,ceil(nt/2)));
% second column up
post(ceil(nt/2)+1:end,2) = fliplr(linspace(0.05,0.8,nt-ceil(nt/2)));

% make the text boxes
for k = 1:nt

  ud.uih.cl.(par_name{k}) = uicontrol('style','text',...
    'Units', 'normalized',...
    'position',post(k,:),...
    'fontsize',14,...
    'Parent',ud.uih.par.panel,...
    'HorizontalAlignment','left',...
    'string',[par_str{k},' = ']);
  ud.uih.cl.([par_name{k} 'c'])=uicontrol('Style','togglebutton',...
    'Units','normalized',...
    'Position',[post(k,1)-0.1 post(k,2)+0.05 0.07 0.07],...
    'Parent',ud.uih.par.panel,...
    'Callback','irf_shock_gui(''manual_input'')');

end

%% Calculate panel
% panel
ud.uih.cl.panel = uipanel('Units', 'normalized',...
  'position',[0.56 0.02 0.44 0.35],...
  'fontsize',14,...
  'Title','Shock & upstream parameters');

% push button
ud.uih.cl.pb = uicontrol('style','push',...
  'Units', 'normalized',...
  'position',[0.2 0.9 0.6 0.1],...
  'fontsize',14,...
  'BackgroundColor',[.85,1,.9],...
  'Parent',ud.uih.cl.panel,...
  'HorizontalAlignment','left',...
  'string','Calculate',...
  'callback','irf_shock_gui(''calc'')');

% parameters to show
par_name = {'n','Vsh','thBn','thVn','r_n','r_b','Ma','Mf','Ms','beta_i','beta_e'};
par_str = par_name;
% "/" is not allowed in variable name
par_str{5} = 'nd/nu'; par_str{6} = 'Bd/Bu';

% positions of text boxes
nt = length(par_name);
post = zeros(nt,4);
% width
post(:,3) = 0.45;
% height
post(:,4) = 0.15;
% first column left
post(1:ceil(nt/2),1) = 0.05;
% second column left
post(ceil(nt/2)+1:end,1) = 0.555;
% first column up
post(1:ceil(nt/2),2) = fliplr(linspace(0.05,0.7,ceil(nt/2)));
% second column up
post(ceil(nt/2)+1:end,2) = fliplr(linspace(0.05,0.7,nt-ceil(nt/2)));

% make the text boxes
for k = 1:nt
  ud.uih.cl.(par_name{k}) = uicontrol('style','text',...
    'Units', 'normalized',...
    'position',post(k,:),...
    'fontsize',14,...
    'Parent',ud.uih.cl.panel,...
    'HorizontalAlignment','left',...
    'string',[par_str{k},' = ']);
end

%ud.params = [];

end


function ud = set_labels(ud)
% set ylabels and legends in panels, aslo sets ylabels
for k = 1:ud.Nin
  ud.ax(k).YLabel.Interpreter = 'tex';
  switch ud.ax(k).Tag
    case 'B'
      ud.ax(k).YLabel.String = 'B (nT)';
      irf_legend(ud.ax(k),{'B_x','B_y','B_z'},[0.98,0.98])
    case 'V'
      ud.ax(k).YLabel.String = 'V (km/s)';
      irf_legend(ud.ax(k),{'V_x','V_y','V_z'},[0.98,0.98])
    case 'n'
      ud.ax(k).YLabel.String = 'n (cm^{-3})';
      ud.ax(k).YLim(1) = 0;
    case 'Ti'
      ud.ax(k).YLabel.String = 'Ti (eV)';
      ud.ax(k).YLim(1) = 0;
    case 'Te'
      ud.ax(k).YLabel.String = 'Te (eV)';
      ud.ax(k).YLim(1) = 0;
    otherwise
      ud.ax(k).YLabel.String = ':P';
  end
end
end


function [ud] = clickt(ud,str) % click time, str is "u" or "d" or "f"
% Give some log, should probably be "mark upstream/downstream"
irf.log('w',['Mark ',num2str(nargout), ' time interval for averaging.'])
% Click times
[t,~] = ginput(2);
t = sort(t);
fig = gcf;
t = t+fig.UserData.t_start_epoch;
ud.(['t',str{1}]) = t;
end


function [ud] = mark_times(ud) % mark times, str is 'u' or 'd' or 'f'
ucol = [0.7,0.7,0];
dcol = [0,0.7,0.7];
fcol = [0.7,0.35,0.35];
for k = 1:ud.Nin
  delete(findall(ud.ax(k).Children,'Tag','irf_pl_mark'));
end
pause(0.001)

for k = 1:ud.Nin % new feature in axescheck, does not allow handle arrays
  irf_pl_mark(ud.ax(k),ud.tu',ucol)
  irf_pl_mark(ud.ax(k),ud.td',dcol)
  irf_pl_mark(ud.ax(k),ud.tf',fcol)
end
end


function [ud] = set_methods(ud)
% '--' will crash the program, the order is very important
nmet = {'mc','vc','mx1','mx2','mx3','--',...
  'farris','slho','per','fa4o','fan4o','foun'};
vmet = {'gt','mf','sb','mo'};

if ud.uih.mt.npu.Value <= 5 || isfield(ud.params,'R') % not a model or with position
  ud.normal_method = nmet{ud.uih.mt.npu.Value};
  ud.vel_method = vmet{ud.uih.mt.vpu.Value};
elseif ud.uih.mt.npu.Value == 6  % -- is not actually a method
  msgbox('Not a method')
else  % if no position and model give error
  msgbox('Model not applicable whith no spacecraft position data in input.',...
    'No position')
end

% special for mach method
if ud.uih.mt.mpu.Value == 1 % spacecraft frame
  ud.params.ref_sys = 'sc';
else% nif
  ud.params.ref_sys = 'nif';
  ud.params.nvec = ud.shp.nvec.n.(ud.normal_method);
  if ud.uih.mt.mpu.Value == 2% Vsh is determined from velocity method
    ud.params.Vsh = ud.shp.nvec.Vsh.(ud.vel_method).(ud.normal_method);
  else % Vsh = 0
    ud.params.Vsh = 0;
  end
end
end


function [ud] = display_vals(ud) % print up/downstream values
% Now displays norm of vectors because it is shorter.

% norm of magnetic field vector
if ~isfield(ud.uih.cl,'Bum')
  temp=round((ud.params.Bu),2);
  ud.uih.cl.Bu.String = ['Bu = (',[num2str(temp(1)) ',' num2str(temp(2)) ',' num2str(temp(3))],') nT'];
  clear temp
end
if ~isfield(ud.uih.cl,'Bud')
  temp=round((ud.params.Bd),2);
  ud.uih.cl.Bd.String = ['Bd = (',[num2str(temp(1)) ',' num2str(temp(2)) ',' num2str(temp(3))],') nT'];
  clear temp
end
% destiny
if ~isfield(ud.uih.cl,'num')
  ud.uih.cl.nu.String = ['nu = ',num2str(round(ud.params.nu,2)),' /cc'];
end
if ~isfield(ud.uih.cl,'ndm')
  ud.uih.cl.nd.String = ['nd = ',num2str(round(ud.params.nd,2)),' /cc'];
end
% velocity
if ~isfield(ud.uih.cl,'Vum')
  temp=round((ud.params.Vu));
  ud.uih.cl.Vu.String = ['Vu = (',[num2str(temp(1)) ',' num2str(temp(2)) ',' num2str(temp(3))],') km/s'];
  clear temp
end
if ~isfield(ud.uih.cl,'Vdm')
  temp=round((ud.params.Vd));
  ud.uih.cl.Vd.String = ['Vd = (',[num2str(temp(1)) ',' num2str(temp(2)) ',' num2str(temp(3))],') km/s'];
  clear temp
end
% electron temperature
if ~isfield(ud.uih.cl,'Teum')
  ud.uih.cl.Teu.String = ['Teu = ',num2str(round(ud.params.Teu,2)),' eV'];
end
if ~isfield(ud.uih.cl,'Tedm')
  ud.uih.cl.Ted.String = ['Ted = ',num2str(round(ud.params.Ted,2)),' eV'];
end
% ion temperature
if ~isfield(ud.uih.cl,'Tium')
  ud.uih.cl.Tiu.String = ['Tiu = ',num2str(round(ud.params.Tiu,2)),' eV'];
end
if ~isfield(ud.uih.cl,'Tidm')
  ud.uih.cl.Tid.String = ['Tid = ',num2str(round(ud.params.Tid,2)),' eV'];
end
end


function [ud] = display_prop(ud) % write out results

% normal vector
nstr = ['[',num2str(round(ud.shp.nvec.n.(ud.normal_method)(1),2)),',',...
  num2str(round(ud.shp.nvec.n.(ud.normal_method)(2),2)),',',...
  num2str(round(ud.shp.nvec.n.(ud.normal_method)(3),2)),']'];
ud.uih.cl.n.String = ['n = ',nstr];
% shock angle
ud.uih.cl.thBn.String = ['thBn = ',num2str(round(ud.shp.nvec.thBn.(ud.normal_method),2))];
% flow incidence angle
ud.uih.cl.thVn.String = ['thVn = ',num2str(round(ud.shp.nvec.thVn.(ud.normal_method),2))];
% flow incidence angle
ud.uih.cl.thVn.String = ['thVn = ',num2str(round(ud.shp.nvec.thVn.(ud.normal_method),2))];
% density compression ratio
ud.uih.cl.r_n.String = ['nd/nu = ',num2str(round(ud.shp.data.nd/ud.shp.data.nu,2))];
% magnetic compression ratio
ud.uih.cl.r_b.String = ['Bd/Bu = ',num2str(round(norm(ud.shp.data.Bd)/norm(ud.shp.data.Bu),2))];
% Alfven mach
if isfield(ud.shp.par,'Mau')
  ud.uih.cl.Ma.String = ['Ma = ',num2str(round(ud.shp.par.Mau,2))];
end
% MS mach
if isfield(ud.shp.par,'Mfu')
  ud.uih.cl.Mf.String = ['Mf = ',num2str(round(ud.shp.par.Mfu,2))];
end
% sound mach
if isfield(ud.shp.par,'Mfu')
  ud.uih.cl.Ms.String = ['Ms = ',num2str(round(ud.shp.par.Msu,2))];
end
% ion beta
if isfield(ud.shp.par,'biu')
  ud.uih.cl.beta_i.String = ['beta_i = ',num2str(round(ud.shp.par.biu,2))];
end
% electron beta
if isfield(ud.shp.par,'beu')
  ud.uih.cl.beta_e.String = ['beta_e = ',num2str(round(ud.shp.par.beu,2))];
end
% shock speed
ud.uih.cl.Vsh.String = ['Vsh = ',num2str(round(ud.shp.nvec.Vsh.(ud.vel_method).(ud.normal_method),0)),' km/s'];
end


function [ud] = get_avg_field(ud,par,fin)
%AVG_FIELD Calculates average value in a time interval
%   out = AVG_FIELD(ud,par,fin)
%
%   Examples:
%           par =
%               B: [1x1 TSeries]
%
%       >> avg = avg_field(ud,par,{'avg'});
%           out =
%               Bavg: [1.52,2.53,-0.22]
%
% -------------------------------------
%           par =
%               B: [1x1 TSeries]
%               V: [1x1 TSeries]
%
%       >> avg = avg_field(ud,par,{'u','d'});
%           out =
%               Bu: ...
%               Bd: ...
%               Vu: ...
%               Vd: ...

% number of parameter inputs
fnp = fieldnames(par);
nin = numel(fnp);

N = numel(fin);

% Plot all inputs only first time
if isempty(ud.ax(1).Tag)
  for k = 1:nin
    if ~strcmpi(fnp{k},'r')
      irf_plot(ud.ax(k),par.(fnp{k}));
      ud.ax(k).Tag = fnp{k};
    end
  end
end

% Set color order (this is done twice for some reason)
color_order = zeros(nargout*3,3);
for i = 1:3:nargout*3
  color_order(i:i+2,:) = [0,0,0; 0,0,1; 1,0,0];
end
for k = 1:nin % set color order
  ud.ax(k).ColorOrder = color_order;
end

% Nested loop for calculating averages
for i = 1:N
  for k = 1:nin
    % variable
    X = par.(fnp{k});
    % find time indicies
    idt = interp1(X.time.epochUnix,1:length(X.time.epochUnix),ud.(['t',fin{i}]),'nearest');
    % handle out-of-panel clicks
    if isnan(idt(1)) && t(1) < X.time(1).epochUnix
      idt(1) = 1;
    end
    if isnan(idt(2)) && t(2) > X.time(end).epochUnix
      idt(2) = length(X.time.epochUnix);
    end
    tinti = X.time(idt).epochUnix';
    hold(ud.ax(k),'on')
    % Calculate average
    x_avg = double(nanmean(X.data(idt(1):idt(2),:)));
    % Plot average as solid lines lines
    delete(findall(ud.ax(k).Children,'Tag',['avg_mark',fin{i}]))
    irf_plot(ud.ax(k),[tinti',[x_avg;x_avg]],'Tag',['avg_mark',fin{i}])
    ud.ax(k).ColorOrder = color_order;
    ud.params.([fnp{k},fin{i}]) = x_avg;
  end
end

for k = 1:nin % set various things
  % remove XTickLabels for all but the last panel
  if k<nin; ud.ax(k).XTickLabels = ''; end
  % thicker lines in plots
  ud.ax(k).LineWidth = 1.15;
end

if ~isempty(ud.params)
  ud.params = orderfields(ud.params);
end
end


function [ud] = set_omni(ud)
% sets OMNI-valu as the average value in the upstream time interval (entire
% interval if not set)
ud.use_omni.B = ud.uih.omni.pbB.Value;
ud.use_omni.n = ud.uih.omni.pbN.Value;
ud.use_omni.V = ud.uih.omni.pbV.Value;
ud.use_omni.Ti = ud.uih.omni.pbTi.Value;


% Read omni data if not already done and if a box is checked
if ~isfield(ud,'omnidata') && (ud.use_omni.B || ud.use_omni.n || ud.use_omni.V || ud.use_omni.Ti)
  irf.log('w','Reading OMNI data from internet...')
  ud.omnidata = [];
  tint = ud.scd.B.time([1,end])+[-60,60]*5; % set time interval +- 5 mins
  ud.omnidata = irf_get_data(tint,'bx,by,bz,n,vx,vy,vz,t','omni_min');
  % Re-correct Vy for abberation
  %   ud.omnidata(:,6) = ud.omnidata(:,7)+29.8;
  % change temperature units from K to eV
  u = irf_units;
  ud.omnidata(:,9) = ud.omnidata(:,9)*u.kB/u.e;
  irf.log('w','Done reading OMNI data.')
end

if ud.use_omni.B
  %   Bomni = nanmean(irf_resamp(ud.omnidata(:,1:4),ud.tu),1);
  %   ud.params.Bu = Bomni(2:4);

  Bomni = nanmean(irf.ts_vec_xyz(EpochUnix(ud.omnidata(:,1)),ud.omnidata(:,2:4)).tlim(EpochUnix(ud.tu)+[-60 60]).data,1);
  ud.params.Bu = Bomni;
else; ud.params.Bu = ud.sc_up.Bu;
end
if ud.use_omni.n
  %   nomni = nanmean(irf_resamp(ud.omnidata(:,[1,5]),ud.tu),1);
  %   ud.params.nu = nomni(2);

  nomni = nanmean(irf.ts_scalar(EpochUnix(ud.omnidata(:,1)),ud.omnidata(:,5)).tlim(EpochUnix(ud.tu)+[-60 60]).data,1);
  ud.params.nu = nomni;
else; ud.params.nu = ud.sc_up.nu;
end
if ud.use_omni.V
  %   Vomni = nanmean(irf_resamp(ud.omnidata(:,[1,6:8]),ud.tu),1);
  %   ud.params.Vu = Vomni(2:4);

  Vomni = nanmean(irf.ts_vec_xyz(EpochUnix(ud.omnidata(:,1)),ud.omnidata(:,6:8)).tlim(EpochUnix(ud.tu)+[-60 60]).data,1);
  ud.params.Vu = Vomni;
else; ud.params.Vu = ud.sc_up.Vu;
end
if ud.use_omni.Ti
  %   Tomni = nanmean(irf_resamp(ud.omnidata(:,[1,9]),ud.tu),1);
  %   ud.params.Tiu = Tomni(2);

  Tomni = nanmean(irf.ts_scalar(EpochUnix(ud.omnidata(:,1)),ud.omnidata(:,9)).tlim(EpochUnix(ud.tu)+[-60 60]).data,1);
  ud.params.Tiu = Tomni;
else; ud.params.Tiu = ud.sc_up.Tiu;
end
end


function [] = plot_omni(ud)
h = irf_plot(4,'newfigure');

hca = irf_panel(h,'B');
irf_plot(hca,ud.omnidata(:,1:4))
ylabel(hca,'B (nT)')
irf_legend(hca,{'B_x','B_y','B_z'},[0.98,0.98])
hold(hca,'on')
irf_plot(hca,[ud.tu+[-60; 60],[ud.params.Bu;ud.params.Bu]])

hca = irf_panel(h,'n_i');
irf_plot(hca,ud.omnidata(:,[1,5]))
ylabel(hca,'n (cm^{-3})')
hold(hca,'on')
irf_plot(hca,[ud.tu+[-60; 60],[ud.params.nu;ud.params.nu]])

hca = irf_panel(h,'V');
irf_plot(hca,ud.omnidata(:,[1,6:8]))
ylabel(hca,'V (km/s)')
irf_legend(hca,{'V_x','V_y','V_z'},[0.98,0.98])
hold(hca,'on')
irf_plot(hca,[ud.tu+[-60; 60],[ud.params.Vu;ud.params.Vu]])

hca = irf_panel(h,'T');
irf_plot(hca,ud.omnidata(:,[1,9]))% needs SI
ylabel(hca,'Ti (eV)')
hold(hca,'on')
irf_plot(hca,[ud.tu+[-60; 60],[ud.params.Tiu;ud.params.Tiu]])

for i = 1:4
  irf_pl_mark(h(i),ud.tu'+[-60 60],[0.7,0.7,0])
end
irf_zoom(h(1:end),'x',EpochUnix(ud.omnidata(:,1)));
end

function [ud] = manual_input(ud)

par_name = {'Bu','nu','Vu','Teu','Tiu','Bd','nd','Vd','Ted','Tid'};
Uns = {'nT','/cc','Km/s','eV','eV','nT','/cc','Km/s','eV','eV'};
par_str = par_name;

% positions of text boxes
nt = length(par_name);
post = zeros(nt,4);
% width
post(:,3) = 0.4;
% height
post(:,4) = 0.15;
% first column left
post(1:ceil(nt/2),1) = 0.15;
% second column left
post(ceil(nt/2)+1:end,1) = 0.65;
% first column up
post(1:ceil(nt/2),2) = fliplr(linspace(0.05,0.8,ceil(nt/2)));
% second column up
post(ceil(nt/2)+1:end,2) = fliplr(linspace(0.05,0.8,nt-ceil(nt/2)));

% make the text boxes
for k = 1:nt

  if ud.uih.cl.([par_name{k} 'c']).Value && ~isfield(ud.uih.cl,[par_name{k} 'm'])
    ud.uih.cl.(par_name{k}).String=[par_str{k},' = '];


    ud.uih.cl.([par_name{k} 'm'])=uicontrol('Style','edit',...
      'Units','normalized',...
      'Position',[post(k,1)+0.09 post(k,2)+0.01 0.15 0.15],...
      'Parent',ud.uih.par.panel,...
      'String',num2str(round(100*ud.params.(par_name{k}))/100),'Callback','irf_shock_gui(''calc'')');

    ud.uih.cl.([par_name{k} 't2'])=uicontrol('Style','text',...
      'Units','normalized',...
      'Position',[post(k,1)+0.24 post(k,2)+0.02 0.15 0.1],...
      'Parent',ud.uih.par.panel,...
      'String',Uns{k},'FontSize',12,'HorizontalAlignment','left');

  elseif ud.uih.cl.([par_name{k} 'c']).Value && isfield(ud.uih.cl,[par_name{k} 'm'])
    ud.params.(par_name{k}) = cell2mat(textscan(ud.uih.cl.([par_name{k} 'm']).String,'%f'))';

  elseif isfield(ud.uih.cl,[par_name{k} 'm'])
    delete(ud.uih.cl.([par_name{k} 'm']));
    ud.uih.cl = rmfield(ud.uih.cl, [par_name{k} 'm']);
    delete(ud.uih.cl.([par_name{k} 't2']));
    ud.uih.cl = rmfield(ud.uih.cl, [par_name{k} 't2']);
    ud = get_avg_field(ud,ud.scd,{'u'});
    ud = set_omni(ud);
    ud = display_vals(ud);
  end

end


end

