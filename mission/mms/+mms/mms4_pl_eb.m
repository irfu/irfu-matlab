function h = mms4_pl_eb(Tint, dMode)
%MMS.MMS4_PL_EB  Summary plot - E & B at 4 MS S/C
%
%  h = MMS.MMS4_PL_EB(Tint, [MODE])
%
%  MODE - one of 'fast' (default) or 'brst'

%Tint = irf.tint('2015-05-24T02:10:00Z/2015-05-24T02:30:00Z');

if nargin<2, dMode = 'fast'; end
switch dMode
  case {'brst','fast'}
  otherwise
    error('MODE must be one of: ''brst'', ''fast''')
end

%% Load data
tic
for mmsId = 1:4
  fprintf('Loading MMS%d\n',mmsId);
  %% FGM & EDP
  B_dmpa_fgm_srvy = mms.get_data('B_dmpa_fgm_srvy_l2',Tint,mmsId);
  if isempty(B_dmpa_fgm_srvy)
    irf.log('warning','loading L2pre DFG')
    B_dmpa_fgm_srvy = mms.get_data('B_dmpa_dfg_srvy_l2pre',Tint,mmsId);
    if isempty(B_dmpa_fgm_srvy)
      irf.log('warning','loading QL DFG')
      B_dmpa_fgm_srvy = mms.get_data('B_dmpa_dfg_srvy_ql',Tint,mmsId); %#ok<NASGU>
    end
  end
  E_dsl_edp = mms.get_data(['E_dsl_edp_' dMode '_l2'],Tint,mmsId);
  if isempty(E_dsl_edp)
    irf.log('warning','loading QL DCE')
    E_dsl_edp = mms.get_data(['E_dsl_edp_' dMode '_ql'],Tint,mmsId); %#ok<NASGU>
  end

  %{
  E2d_dsl_edp = mms.get_data(['E2d_dsl_edp_' dMode '_l2pre'],Tint,mmsId);
  if isempty(E2d_dsl_edp)
    irf.log('warning','loading QL DCE2d')
    E2d_dsl_edp = mms.get_data(['E2d_dsl_edp_' dMode '_ql'],Tint,mmsId); %#ok<NASGU>
  end
  %}

  V_edp = mms.get_data(['V_edp_' dMode '_l2'],Tint,mmsId);
  if isempty(V_edp)
    irf.log('warning','loading SITL DCV')
    V_edp = mms.get_data('V_edp_fast_sitl',Tint,mmsId); %#ok<NASGU>
  end

  R_gse = mms.get_data('R_gse',Tint,mmsId); %#ok<NASGU>

  c_eval([...
    'E? = E_dsl_edp;'...
    'P? = V_edp;'...
    'B? = B_dmpa_fgm_srvy;'...
    'R? = R_gse;'],...
    mmsId)
  clear B_dmpa_fgm_srvy E_dsl_edp V_edp R_gse
end
fprintf('Data loaded\n');
if ~isempty(R1), gseR = [R1.time.epochUnix double(R1.data(:,1:3))];
elseif ~isempty(R2), gseR = [R2.time.epochUnix double(R2.data(:,1:3))];
elseif ~isempty(R3), gseR = [R3.time.epochUnix double(R3.data(:,1:3))];
elseif ~isempty(R4), gseR = [R4.time.epochUnix double(R4.data(:,1:3))];
else, gseR = [];
end
toc
%% Plot
tic
% define Cluster colors
mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
% Official MMS colors
%mmsColors=[0 0 0; .8 .4 0 ; 0 0.6 0.5 ; 0.35 0.7 .9];

h = irf_figure(8);

if 1
  hca = irf_panel('B'); set(hca,'ColorOrder',mmsColors)
  irf_pl_tx(hca,'abs(B?)',1)
  ylabel(hca,'|B| [nT]')
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98, 0.1],'color','cluster');

  hca = irf_panel('Bx'); set(hca,'ColorOrder',mmsColors)
  irf_pl_tx(hca,'B?',1)
  ylabel(hca,'Bx [nT]')

  hca = irf_panel('By'); set(hca,'ColorOrder',mmsColors)
  irf_pl_tx(hca,'B?',2)
  ylabel(hca,'By [nT]')
end

hca = irf_panel('Bz'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',3)
ylabel(hca,'Bz [nT]')

hca = irf_panel('Ex'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'E?',1)
ylabel(hca,'Ex [mV/m]')

hca = irf_panel('Ey'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'E?',2)
ylabel(hca,'Ey [mV/m]')

if 1
  hca = irf_panel('Ez'); set(hca,'ColorOrder',mmsColors)
  irf_pl_tx(hca,'E?',3)
  ylabel(hca,'Ez [mV/m]')
end

if 1
  hca = irf_panel('ScPot'); set(hca,'ColorOrder',mmsColors)
  irf_pl_tx(hca,'P?')
  ylabel(hca,'ScPot [V]')
end


irf_zoom(h,'x',Tint)
irf_plot_axis_align(h)
add_position(h(end),gseR, 'obj','earth')
xlabel(h(end),'')
title(h(1),Tint.start.utc)

toc
return

%% B/V
c_eval('if ~isempty(E?) && ~isempty(B?), EdB? = irf_edb(irf.ts_vec_xy(E?.time,E?.data(:,1:2)),B?); EdB?.units = ''mV/m''; VExB? = irf_e_vxb(EdB?,B?,-1); else VExB?=[]; end')  %#ok<UNRCH>

h = irf_figure(8);

hca = irf_panel('B'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'abs(B?)',1)
ylabel(hca,'|B| [nT]')
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98, 0.1],'color','cluster');

hca = irf_panel('Bx'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',1)
ylabel(hca,'Bx [nT]')

hca = irf_panel('By'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',2)
ylabel(hca,'By [nT]')

hca = irf_panel('Bz'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',3)
ylabel(hca,'Bz [nT]')

hca = irf_panel('Vx'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'VExB?',1)
ylabel(hca,'Vx [km/s]')

hca = irf_panel('Vy'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'VExB?',2)
ylabel(hca,'Vy [km/s]')

hca = irf_panel('Vz'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'VExB?',3)
ylabel(hca,'Vz [km/s]')

if 1
  hca = irf_panel('ScPot'); set(hca,'ColorOrder',mmsColors)
  irf_pl_tx(hca,'P?')
  ylabel(hca,'-ScPot [V]')
end

irf_zoom(h,'x',Tint)
irf_plot_axis_align(h)
add_position(h(end),gseR)
xlabel(h(end),'')
title(h(1),Tint.start.utc)