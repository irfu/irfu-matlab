%MMS.MMS4_PL_EB  Summary plot - E & B at 4 MS S/C

tint = irf.tint('2015-05-24T02:10:00Z/2015-05-24T02:30:00Z');

%% Load data
load /data/mms/irfu/mmsR.mat
epoTmp = EpochTT(R.time);
gsmR1 = [epoTmp.epochUnix R.gsmR1];

%E1 = mms.db_get_ts('mms1_edp_comm_ql_dce2d','mms1_edp_dce_xyz_dsl',tint);
%P1 = mms.db_get_ts('mms1_edp_comm_l2_scpot','mms1_edp_psp',tint);
%B1 = mms.db_get_ts('mms1_dfg_srvy_ql','mms1_dfg_srvy_gsm_dmpa',tint);

for scId = 1:4
  fprintf('Loading MMS%d\n',scId);
  c_eval([...
    'E? = mms.db_get_ts(''mms?_edp_comm_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);'...
    'P? = mms.db_get_ts(''mms?_edp_comm_l2_scpot'',''mms?_edp_psp'',tint);'...
    'B? = mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_gsm_dmpa'',tint);'],...
    scId)
end
fprintf('Data loaded\n');

%% Plot
% define Cluster colors
mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
% Official MMS colors
%mmsColors=[0 0 0; .8 .4 0 ; 0 0.6 0.5 ; 0.35 0.7 .9];

h = irf_plot(7,'newfigure');

hca = irf_panel('Bx'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{B1.x,B2.x,B3.x,B4.x},'comp')
ylabel(hca,'Bx [nT]')

hca = irf_panel('By'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{B1.y,B2.y,B3.y,B4.y},'comp')
ylabel(hca,'By [nT]')

hca = irf_panel('Bz'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{B1.z,B2.z,B3.z,B4.z},'comp')
ylabel(hca,'Bz [nT]')

if 1
hca = irf_panel('ScPot'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{P1,P2,P3,P4},'comp')
ylabel(hca,'-ScPot [V]')
end

hca = irf_panel('Ex'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{E1.x,E2.x,E3.x,E4.x},'comp')
ylabel(hca,'Ex [mV/m]')

hca = irf_panel('Ey'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{E1.y,E2.y,E3.y,E4.y},'comp')
ylabel(hca,'Ey [mV/m]')

hca = irf_panel('Ez'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,{E1.z,E2.z,E3.z,E4.z},'comp')
ylabel(hca,'Ez [mV/m]')


irf_zoom(h,'x',tint)
irf_plot_axis_align(h)
add_position(h(end),gsmR1)
xlabel(h(end),'')
title(h(1),tint.start.utc)