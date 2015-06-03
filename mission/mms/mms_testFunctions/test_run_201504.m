%% Init
data_root='/data/mms';
if ~exist(data_root,'dir')
  error('DATA_ROOT(%s) does not exist',data_root)
end
if ismac,
  setenv('CDF_BASE','/Applications/cdf35_0-dist/')
  outDir = '/Users/yuri/Documents/MATLAB/MMS/test_proc/';
else outDir = data_root;
end
cd(outDir)
if ~exist('log','dir'), mkdir('log'), end
if ~exist('out','dir'), mkdir('out'), end
setenv('DROPBOX_ROOT', [outDir filesep 'out'])
setenv('DATA_PATH_ROOT', [outDir filesep 'out'])
setenv('LOG_PATH_ROOT', [outDir filesep 'log'])
MMS_CONST=mms_constants;

%% Define time
%tint = irf.tint('2015-04-16T00:00:00Z/2015-04-16T06:00:00Z');
%tint = irf.tint('2015-04-16T18:00:00Z/2015-04-16T23:59:59Z');
tint = irf.tint('2015-05-15T00:00:00Z/2015-05-15T05:59:59Z');
%tint = irf.tint('2015-04-20T18:00:00Z/2015-04-20T23:59:59Z');
mmsId = 'mms4'; 

prf = [data_root filesep mmsId]; utc = tint.start.toUtc(); 
mo = utc(6:7); yyyy=utc(1:4); day=utc(9:10); hh=utc(12:13); mm=utc(15:16);

li = mms.db_list_files([mmsId '_fields_hk_l1b_101'],tint); if length(li)>1, error('li>1'), end
HK_101_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_fields_hk_l1b_105'],tint); if length(li)>1, error('li>1'), end
HK_105_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_fields_hk_l1b_10e'],tint); if length(li)>1, error('li>1'), end
HK_10E_File = [li.path filesep li.name];
DCE_File  = [prf '/edp/comm/l1b/dce128/' yyyy '/' mo '/' mmsId ...
  '_edp_comm_l1b_dce128_' yyyy mo day hh mm '00_v0.8.0.cdf'];
DCV_File  = [prf '/edp/comm/l1b/dcv128/' yyyy '/' mo '/' mmsId ...
  '_edp_comm_l1b_dcv128_' yyyy mo day hh mm '00_v0.8.0.cdf'];

%% Test QL
irf.log('log_out','screen'), irf.log('notice')
procId = MMS_CONST.SDCProc.ql; procName='QL'; scId=str2double(mmsId(end));
tmMode=MMS_CONST.TmMode.comm; samplerate = MMS_CONST.Samplerate.comm_128;

mms_sdp_datamanager('init',...
  struct('scId',scId,'tmMode',tmMode,'procId',procId,...
  'samplerate',samplerate));
mms_sdp_load(HK_10E_File,'hk_10e');
mms_sdp_load(HK_105_File,'hk_105');
mms_sdp_load(HK_101_File,'hk_101');
mms_sdp_load(DCE_File,'dce');
mms_sdp_load(DCV_File,'dcv');

dce = mms_sdp_datamanager('dce');
probe2sc_pot = mms_sdp_datamanager('probe2sc_pot');
phase = mms_sdp_datamanager('phase');
spinfits = mms_sdp_datamanager('spinfits');
delta_off = mms_sdp_datamanager('delta_off');
dce_xyz_dsl = mms_sdp_datamanager('dce_xyz_dsl');

epochE = EpochTT2000(dce.time).toEpochUnix().epoch;
epochS = EpochTT2000(spinfits.time).toEpochUnix().epoch;

%% Summary plot
E_YLIM = 7;

figure(71), clf
h = irf_plot(4);

hca = irf_panel('E');
irf_plot(hca,[epochE double(dce.e12.data) double(dce.e34.data)])
hold(hca,'on'), comp = 1;
irf_plot(hca,[epochS double(spinfits.sfit.e12(:,comp)) double(spinfits.sfit.e34(:,comp))])
ylabel(hca,'E spin [mV/m]')
title(hca,mmsId), set(hca,'YLim',49*[-1 1])

if 0
hca = irf_panel('Phase');
irf_plot(hca,[epochE double(phase.data)])
ylabel(hca,'Phase [deg]'), set(hca,'YLim',[0 360])
end

hca = irf_panel('Ex'); comp = 1;
irf_plot(hca,{[epochE double(dce_xyz_dsl.data(:,comp))],...
  [epochS double(spinfits.sfit.e12(:,comp+1))-real(delta_off)],...
  [epochS double(spinfits.sfit.e34(:,comp+1))]...
  },'comp')
ylabel(hca,'Ex DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('Ey'); comp = 2;
irf_plot(hca,{[epochE double(dce_xyz_dsl.data(:,comp))],...
  [epochS double(spinfits.sfit.e12(:,comp+1))-imag(delta_off)],...
  [epochS double(spinfits.sfit.e34(:,comp+1))]...
  },'comp')
ylabel(hca,'Ey DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('V');
irf_plot(hca,[epochE double(probe2sc_pot.data)])
ylabel(hca,'P2ScPot [V]'), set(hca,'YLim',[-14 0])

irf_plot_ylabels_align(h), irf_zoom(h,'x',epochE([1 end])')

%% Delta offsets
Delta_p12_p34 = double(spinfits.sfit.e12(:,2:3)) - ...
  double(spinfits.sfit.e34(:,2:3));

figure(73), clf
h = irf_plot(3);
hca = irf_panel('Ex'); 
irf_plot(hca,[epochS Delta_p12_p34(:,1)])
ylabel(hca,'\Delta_{p12,p34}_x [mV/m]')
set(hca,'YLim',[-1.9 0]), title(hca,mmsId)

hca = irf_panel('Ey'); 
irf_plot(hca,[epochS Delta_p12_p34(:,2)])
ylabel(hca,'\Delta_{p12,p34}_y [mV/m]')
set(hca,'YLim',[-.9 .9])

hca = irf_panel('V');
irf_plot(hca,[epochE double(probe2sc_pot.data)])
ylabel(hca,'P2ScPot [V]'), set(hca,'YLim',[-14 0])

irf_plot_ylabels_align(h), irf_zoom(h,'x',epochE([1 end])')


%% Run
mms_sdc_sdp_proc('ql', DCE_File,  DCV_File, HK_10E_File, HK_101_File);
mms_sdc_sdp_proc('scpot', DCE_File,  DCV_File, HK_10E_File, HK_101_File);

mms_sdc_sdp_proc('l2pre', DCE_File,  DCV_File, HK_10E_File, DEFATT_File);

%%
mms_sdc_sdp_proc('l2a','out/mms4_edp_comm_l2pre_dce2d_20150405000000_v0.0.0.cdf')

%% Plot
dce2d=dataobj('out/mms4_edp_comm_l2pre_dce2d_20150405000000_v2.0.0.cdf');
e12=getmat(dce2d,'mms4_edp_dce_spinfit_e12');
e34=getmat(dce2d,'mms4_edp_dce_spinfit_e34');
h = irf_plot({e12,e34},'comp');
irf_zoom(h,'x',irf_time([2015 4 5 12 40 0])+[0 60*45])
ylabel(h(1),'Ex [mV/m]')
ylabel(h(2),'Ey [mV/m]')
ylabel(h(3),'std')
ylabel(h(4),'adc')
irf_plot_ylabels_align(h)
title(h(1),'MMS4')

%% Load B
li = mms.db_list_files([mmsId '_dfg_srvy_ql'],tint);
if length(li)>1, error('li>1'), end
DFG_File = [li.path filesep li.name];
dfg = dataobj(DFG_File);
B = getmat(dfg,[mmsId '_dfg_srvy_gsm_dmpa']);

c_eval('B?=B;',mmsId);
c_eval('diEs?p12=[epochS double(spinfits.sfit.e12)];',mmsId);
c_eval('diEs?p34=[epochS double(spinfits.sfit.e34)];',mmsId);
c_eval('P?=[epochE double(probe2sc_pot.data)];',mmsId);

%% Plot with B

figure(75), clf
h = irf_plot(4,'newfigure');

hca = irf_panel('B');
hTmp = irf_plot(hca,B);
hTmp(1).Color=[0 0 0]; hTmp(2).Color=[0 0.5 0]; hTmp(3).Color=[1 0 0];
hTmp(4).Color=[.5 .5 .5];
set(hca,'ColorOrder',[[0 0 0];[0 0.5 0];[1 0 0];[.5 .5 .5]])
irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.1])
ylabel(hca,'B [nT]')
title(hca,mmsId), set(hca,'YLim',64*[-1 1])

hca = irf_panel('Ex'); comp = 1;
irf_plot(hca,{[epochE double(dce_xyz_dsl.data(:,comp))],...
  [epochS double(spinfits.sfit.e12(:,comp+1))],...
  [epochS double(spinfits.sfit.e34(:,comp+1))]...
  },'comp')
ylabel(hca,'Ex DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('Ey'); comp = 2;
irf_plot(hca,{[epochE double(dce_xyz_dsl.data(:,comp))],...
  [epochS double(spinfits.sfit.e12(:,comp+1))],...
  [epochS double(spinfits.sfit.e34(:,comp+1))]...
  },'comp')
ylabel(hca,'Ey DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('V');
irf_plot(hca,[epochE double(probe2sc_pot.data)])
ylabel(hca,'P2ScPot [V]'), set(hca,'YLim',[-14 0])

irf_plot_ylabels_align(h), irf_zoom(h,'x',epochE([1 end])')
