%% Init
data_root='/data/mms';
  
if ~exist(data_root,'dir')
  error('DATA_ROOT(%s) does not exist',data_root)
end

if ismac,
  setenv('CDF_BASE','/Applications/cdf35_0-dist/')
  outDir = '/Users/yuri/Documents/MATLAB/MMS/test_proc/';
else
  outDir = data_root;
end

cd(outDir)
if ~exist('log','dir'), mkdir('log'), end
if ~exist('out','dir'), mkdir('out'), end
setenv('DROPBOX_ROOT', [outDir filesep 'out'])
setenv('DATA_PATH_ROOT', [outDir filesep 'out'])
setenv('LOG_PATH_ROOT', [outDir filesep 'log'])

MMS_CONST=mms_constants;

%% Define files 20150408
mmsId = 'mms1'; prf = [data_root filesep mmsId];
DCE_File = [prf '/edp/comm/l1b/dce128/2015/04/mms1_edp_comm_l1b_dce128_20150408060000_v0.8.0.cdf'];
DCV_File = [prf '/edp/comm/l1b/dcv128/2015/04/mms1_edp_comm_l1b_dcv128_20150408060000_v0.8.0.cdf'];
HK_101_File = [prf '/fields/hk/l1b/101/2015/04/mms1_fields_hk_l1b_101_20150408_v0.2.0.cdf'];
HK_10E_File = [prf '/fields/hk/l1b/10e/2015/04/mms1_fields_hk_l1b_10e_20150408_v0.1.0.cdf'];
DEFATT_File = []; %'ancillary/mms4/defatt/MMS4_DEFATT_2015097_2015098.V00';

%% Test QL
irf.log('log_out','screen'), irf.log('notice')

procId = MMS_CONST.SDCProc.ql; procName = 'QL';
tmMode=MMS_CONST.TmMode.comm; samplerate = MMS_CONST.Samplerate.comm_128;
scId = str2double(mmsId(end));

mms_sdp_datamanager('init',...
  struct('scId',scId,'tmMode',tmMode,'procId',procId,...
  'samplerate',samplerate));
mms_sdp_load(HK_10E_File,'hk_10e');
mms_sdp_load(HK_101_File,'hk_101');
mms_sdp_load(DCE_File,'dce');
mms_sdp_load(DCV_File,'dcv');

%% Test plot
dce = mms_sdp_datamanager('dce');
dcv = mms_sdp_datamanager('dcv');
sc_pot = mms_sdp_datamanager('sc_pot');
phase = mms_sdp_datamanager('phase');
dce_xyz_dsl = mms_sdp_datamanager('dce_xyz_dsl');
spinfits = mms_sdp_datamanager('spinfits');

%%
ALPHA = 3;
figure(71), clf
h = irf_plot(5);
hca = irf_panel('E');
epochE = EpochTT2000(dce.time).toEpochUnix().epoch;
epochS = EpochTT2000(spinfits.time).toEpochUnix().epoch;
irf_plot(hca,[epochE double(dce.e12.data) double(dce.e34.data)])
hold(hca,'on'), comp = 1;
irf_plot(hca,[epochS double(spinfits.sfit.e12(:,comp)) double(spinfits.sfit.e34(:,comp))])
ylabel(hca,'E spin [mV/m]')
hca = irf_panel('Phase');
irf_plot(hca,[epochE double(phase.data)])
ylabel(hca,'Phase [deg]')
hca = irf_panel('Ex'); comp = 1;
%irf_plot(hca,[epochE double(dce_xyz_dsl.data(:,comp))*ALPHA])
irf_plot(hca,{...
  [epochE double(dce_xyz_dsl.data(:,comp))],...
  [epochS double(spinfits.sfit.e12(:,comp+1))],...
  [epochS double(spinfits.sfit.e34(:,comp+1))]...
  },'comp')
ylabel(hca,'Ex DSL [mV/m]')
hca = irf_panel('Ey'); comp = 2;
%irf_plot(hca,[epochE double(dce_xyz_dsl.data(:,comp))*ALPHA])
irf_plot(hca,{...
  [epochE double(dce_xyz_dsl.data(:,comp))],...
  [epochS double(spinfits.sfit.e12(:,comp+1))],...
  [epochS double(spinfits.sfit.e34(:,comp+1))]...
  },'comp')
ylabel(hca,'Ey DSL [mV/m]')
hca = irf_panel('V');
irf_plot(hca,[epochE double(sc_pot.data)])
ylabel(hca,'ScPot [V]')
irf_plot_ylabels_align(h)
irf_zoom(h,'x',epochE([1 end])')

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
