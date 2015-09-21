function h = mms4_pl_eb(Tint)
%MMS.MMS4_PL_EB  Summary plot - E & B at 4 MS S/C
%
%  H = MMS.MMS4_PL_EB(Tint)

%Tint = irf.tint('2015-05-24T02:10:00Z/2015-05-24T02:30:00Z');

%% Load data
%load /data/mms/irfu/mmsR.mat
%epoTmp = EpochTT(R.time);
%gsmR1 = [epoTmp.epochUnix R.gsmR1];
%gsmR1(isnan(gsmR1(:,4)),:) = [];

for scId = 1:4
  fprintf('Loading MMS%d\n',scId);
  c_eval([...
    'E? = mms.db_get_ts(''mms?_edp_fast_ql_dce'',''mms?_edp_dce_xyz_dsl'',Tint);'...
    'P? = mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_psp'',Tint);'...
    'B? = mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_gsm_dmpa'',Tint);'...
    'R? = mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_ql_pos_gse'',Tint);'],...
    scId)
end
fprintf('Data loaded\n');
gseR1 = [R1.time.epochUnix double(R1.data(:,1:3))];

%% Plot
% define Cluster colors
mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
% Official MMS colors
%mmsColors=[0 0 0; .8 .4 0 ; 0 0.6 0.5 ; 0.35 0.7 .9];

h = irf_plot(7,'newfigure');

hca = irf_panel('Bx'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',1)
ylabel(hca,'Bx [nT]')

hca = irf_panel('By'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',2)
ylabel(hca,'By [nT]')

hca = irf_panel('Bz'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'B?',3)
ylabel(hca,'Bz [nT]')

hca = irf_panel('Ex'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'E?',1)
ylabel(hca,'Ex [mV/m]')

hca = irf_panel('Ey'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'E?',2)
ylabel(hca,'Ey [mV/m]')

hca = irf_panel('Ez'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'E?',3)
ylabel(hca,'Ez [mV/m]')

if 1
hca = irf_panel('ScPot'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'P?')
ylabel(hca,'-ScPot [V]')
end


irf_zoom(h,'x',Tint)
irf_plot_axis_align(h)
add_position(h(end),gseR1)
xlabel(h(end),'')
title(h(1),Tint.start.utc)

return

%% B/V
c_eval('if ~isempty(E?) && ~isempty(B?), EdB? = irf_edb(irf.ts_vec_xy(E?.time,E?.data(:,1:2)),B?); EdB?.units = ''mV/m''; VExB? = irf_e_vxb(EdB?,B?,-1); else VExB?=[]; end') 

h = irf_plot(7,'newfigure');

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
add_position(h(end),gsmR1)
xlabel(h(end),'')
title(h(1),Tint.start.utc)