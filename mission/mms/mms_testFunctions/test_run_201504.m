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

if 0
%% Define files 20150408
mmsId = 'mms1'; prf = [data_root filesep mmsId];  %#ok<UNRCH>
DCE_File = [prf '/edp/comm/l1b/dce128/2015/04/mms1_edp_comm_l1b_dce128_20150408060000_v0.8.0.cdf'];
DCV_File = [prf '/edp/comm/l1b/dcv128/2015/04/mms1_edp_comm_l1b_dcv128_20150408060000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms1_fields_hk_l1b_101_20150408_v0.3.0.cdf'];
HK_105_File = [prf '/fields/hk/l1b/105/2015/04/mms1_fields_hk_l1b_105_20150408_v0.1.0.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms1_fields_hk_l1b_10e_20150408_v0.1.0.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';
end
if 0
%% Define files 20150411-00:00
mmsId = 'mms3'; prf = [data_root filesep mmsId];  %#ok<UNRCH>
DCE_File = [prf '/edp/comm/l1b/dce128/2015/04/mms3_edp_comm_l1b_dce128_20150411000000_v0.8.0.cdf'];
DCV_File = [prf '/edp/comm/l1b/dcv128/2015/04/mms3_edp_comm_l1b_dcv128_20150411000000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms3_fields_hk_l1b_101_20150411_v0.3.1.cdf'];
HK_105_File = [prf '/fields/hk/l1b/105/2015/04/mms3_fields_hk_l1b_105_20150411_v0.1.1.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms3_fields_hk_l1b_10e_20150411_v0.1.1.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';
end
if 0
%% Define files 20150413-00:00
mmsId = 'mms2'; prf = [data_root filesep mmsId]; %#ok<UNRCH>
DCE_File  = [prf '/edp/comm/l1b/dce128/2015/04/mms2_edp_comm_l1b_dce128_20150413000000_v0.8.0.cdf'];
DCV_File  = [prf '/edp/comm/l1b/dcv128/2015/04/mms2_edp_comm_l1b_dcv128_20150413000000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms2_fields_hk_l1b_101_20150413_v0.3.1.cdf'];
HK_105_File = [prf '/fields/hk/l1b/105/2015/04/mms2_fields_hk_l1b_105_20150413_v0.1.1.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms2_fields_hk_l1b_10e_20150413_v0.1.1.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';
end
if 0
%% Define files 20150413-06:00
mmsId = 'mms2'; prf = [data_root filesep mmsId]; %#ok<UNRCH>
DCE_File  = [prf '/edp/comm/l1b/dce128/2015/04/mms2_edp_comm_l1b_dce128_20150413060000_v0.8.0.cdf'];
DCV_File  = [prf '/edp/comm/l1b/dcv128/2015/04/mms2_edp_comm_l1b_dcv128_20150413060000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms2_fields_hk_l1b_101_20150413_v0.3.1.cdf'];
HK_105_File = [prf '/fields/hk/l1b/105/2015/04/mms2_fields_hk_l1b_105_20150413_v0.1.1.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms2_fields_hk_l1b_10e_20150413_v0.1.1.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';
end
if 0
%% Define files 20150413-12:00
mmsId = 'mms2'; prf = [data_root filesep mmsId]; %#ok<UNRCH>
DCE_File  = [prf '/edp/comm/l1b/dce128/2015/04/mms2_edp_comm_l1b_dce128_20150413120000_v0.8.0.cdf'];
DCV_File  = [prf '/edp/comm/l1b/dcv128/2015/04/mms2_edp_comm_l1b_dcv128_20150413120000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms2_fields_hk_l1b_101_20150413_v0.3.1.cdf'];
HK_105_File = [prf '/fields/hk/l1b/105/2015/04/mms2_fields_hk_l1b_105_20150413_v0.1.1.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms2_fields_hk_l1b_10e_20150413_v0.1.1.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';
end
if 0
%% Define files 20150413-18:00
mmsId = 'mms2'; prf = [data_root filesep mmsId]; %#ok<UNRCH>
DCE_File  = [prf '/edp/comm/l1b/dce128/2015/04/mms2_edp_comm_l1b_dce128_20150413180000_v0.8.0.cdf'];
DCV_File  = [prf '/edp/comm/l1b/dcv128/2015/04/mms2_edp_comm_l1b_dcv128_20150413180000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms2_fields_hk_l1b_101_20150413_v0.3.1.cdf'];
HK_105_File = [prf '/fields/hk/l1b/105/2015/04/mms2_fields_hk_l1b_105_20150413_v0.1.1.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms2_fields_hk_l1b_10e_20150413_v0.1.1.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';
end
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

%% Additional processing
dce = mms_sdp_datamanager('dce');
probe2sc_pot = mms_sdp_datamanager('probe2sc_pot');
phase = mms_sdp_datamanager('phase');
spinfits = mms_sdp_datamanager('spinfits');
dce_xyz_dsl = mms_sdp_datamanager('dce_xyz_dsl');

%% Summary plot
E_YLIM = 7;
epochE = EpochTT2000(dce.time).toEpochUnix().epoch;
epochS = EpochTT2000(spinfits.time).toEpochUnix().epoch;

figure(71), clf
h = irf_plot(5);

hca = irf_panel('E');
irf_plot(hca,[epochE double(dce.e12.data) double(dce.e34.data)])
hold(hca,'on'), comp = 1;
irf_plot(hca,[epochS double(spinfits.sfit.e12(:,comp)) double(spinfits.sfit.e34(:,comp))])
ylabel(hca,'E spin [mV/m]')
title(hca,mmsId), set(hca,'YLim',49*[-1 1])

hca = irf_panel('Phase');
irf_plot(hca,[epochE double(phase.data)])
ylabel(hca,'Phase [deg]'), set(hca,'YLim',[0 360])

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


%% test files
dce = dataobj([data_root DCE_File],'KeepTT2000');
dcv = dataobj([data_root DCV_File],'KeepTT2000');

tE=dce.data.Epoch.data; tV=dcv.data.Epoch.data;
dce_ind = find(tE>tV(1)); dcv_ind = find(tE(end)>tV);

irf_plot([EpochTT2000(tV(dcv_ind)).toEpochUnix().epoch ...
  double(tV(dcv_ind)-tE(dce_ind))*1e-6])
ylabel('\Delta t [ms]')
ylabel('tDCV-tDCE [ms]')
title(mmsId)

%% Run
mms_sdc_sdp_proc('ql', [data_root DCE_File], [data_root DCV_File],...
  [data_root HK_10E_File], [data_root HK_101_File]);

mms_sdc_sdp_proc('l2pre', [data_root DCE_File], [data_root DCV_File],...
  [data_root HK_10E_File], [data_root DEFATT_File]);


%%
mms_sdc_sdp_proc('l2a','out/mms4_edp_comm_l2pre_dce2d_20150405000000_v0.0.0.cdf')

%% Plot
dce2d=dataobj('out/mms4_edp_comm_l2pre_dce2d_20150405000000_v2.0.0.cdf','KeepTT2000');
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
