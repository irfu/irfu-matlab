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

load /data/mms/irfu/mmsR.mat
epocRTmp = EpochTT(R.time);

%% Define time
flagComm = true;
%tint = irf.tint('2015-04-16T00:00:00Z/2015-04-16T06:00:00Z');
%tint = irf.tint('2015-04-16T18:00:00Z/2015-04-16T23:59:59Z');
%tint = irf.tint('2015-05-15T00:00:00Z/2015-05-15T05:59:59Z');
%tint = irf.tint('2015-04-20T18:00:00Z/2015-04-20T23:59:59Z');
%tint = irf.tint('2015-05-06T12:00:00Z/2015-05-06T17:59:59Z');
tint = irf.tint('2015-06-21T00:00:00Z/2015-06-21T05:59:59Z'); 
%tint = irf.tint('2015-06-22T00:00:00Z/2015-06-22T23:59:59Z'); flagComm = false;
mmsId = 'mms2'; 

prf = [data_root filesep mmsId]; utc = tint.start.toUtc(); 
mo = utc(6:7); yyyy=utc(1:4); day=utc(9:10); hh=utc(12:13); mm=utc(15:16);

li = mms.db_list_files([mmsId '_fields_hk_l1b_101'],tint); if length(li)>1, error('li>1'), end
HK_101_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_fields_hk_l1b_105'],tint); if length(li)>1, error('li>1'), end
HK_105_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_fields_hk_l1b_10e'],tint); if length(li)>1, error('li>1'), end
HK_10E_File = [li.path filesep li.name];
if flagComm
  DCE_File  = [prf '/edp/comm/l1b/dce128/' yyyy '/' mo '/' mmsId ...
    '_edp_comm_l1b_dce128_' yyyy mo day hh mm '00_v0.8.0.cdf'];
  DCV_File  = [prf '/edp/comm/l1b/dcv128/' yyyy '/' mo '/' mmsId ...
    '_edp_comm_l1b_dcv128_' yyyy mo day hh mm '00_v0.8.0.cdf'];
else
  DCE_File  = [prf '/edp/fast/l1b/dce/' yyyy '/' mo '/' mmsId ...
    '_edp_fast_l1b_dce_' yyyy mo day '_v1.1.0.cdf'];
  DCV_File = [];
end

gsmR = [epocRTmp.epochUnix R.(['gsmR' mmsId(end)])];

%% Test QL - DMNGR
irf.log('log_out','screen'), irf.log('notice')
procId = MMS_CONST.SDCProc.ql; procName='QL'; scId=str2double(mmsId(end));
tmMode=MMS_CONST.TmMode.comm; samplerate = MMS_CONST.Samplerate.comm_128;
% Initialize DMGR
Dmgr = mms_sdp_dmgr(scId,procId,tmMode,samplerate);
Dmgr.set_param('hk_10e',HK_10E_File);
Dmgr.set_param('hk_105',HK_105_File);
Dmgr.set_param('hk_101',HK_101_File);
Dmgr.set_param('dce',DCE_File);
if ~isempty(DCV_File),Dmgr.set_param('dcv',DCV_File); end
% Process
dce = Dmgr.dce;
probe2sc_pot = Dmgr.probe2sc_pot;
phase = Dmgr.phase;
spinfits = Dmgr.spinfits;
delta_off = Dmgr.delta_off;
dce_xyz_dsl = Dmgr.dce_xyz_dsl;

% Construct TSeries
DceSL = irf.ts_vec_xy(dce_xyz_dsl.time,[dce.e12.data dce.e34.data]);
DceDSL = irf.ts_vec_xyz(dce_xyz_dsl.time,dce_xyz_dsl.data);
Phase = irf.ts_scalar(dce_xyz_dsl.time,phase.data);
AdcOff12 = irf.ts_scalar(spinfits.time,spinfits.sfit.e12(:,1));
AdcOff34 = irf.ts_scalar(spinfits.time,spinfits.sfit.e34(:,1));
Es12 = irf.ts_vec_xy(spinfits.time,spinfits.sfit.e12(:,2:3));
Es34 = irf.ts_vec_xy(spinfits.time,spinfits.sfit.e34(:,2:3));
P2scPot = irf.ts_scalar(probe2sc_pot.time,probe2sc_pot.data);

%% Summary plot
E_YLIM = 7;

figure(71), clf
h = irf_plot(4);

hca = irf_panel('E');
irf_plot(hca,DceSL)
hold(hca,'on') 
irf_plot(hca,{AdcOff12,AdcOff34},'comp')
ylabel(hca,'E SL [mV/m]')
title(hca,mmsId), set(hca,'YLim',49*[-1 1])

if 0
hca = irf_panel('Phase');
irf_plot(hca,Phase)
ylabel(hca,'Phase [deg]'), set(hca,'YLim',[0 360])
end

hca = irf_panel('Ex');
irf_plot(hca,{DceDSL.x,Es12.x-real(delta_off),Es34.x},'comp')
ylabel(hca,'Ex DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('Ey');
irf_plot(hca,{DceDSL.y,Es12.y-imag(delta_off),Es34.y},'comp')
ylabel(hca,'Ey DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('V');
irf_plot(hca,P2scPot)
ylabel(hca,'P2ScPot [V]'), set(hca,'YLim',[-14 0])

%irf_plot_ylabels_align(h), 
irf_zoom(h,'x',DceDSL.time)
%add_position(h(end),gsmR), xlabel(h(end),'')

%% Delta offsets
Delta_p12_p34 = double(spinfits.sfit.e12(:,2:3)) - ...
  double(spinfits.sfit.e34(:,2:3));
epochS = EpochTT(spinfits.time).epochUnix;
epochE = EpochTT(dce.time).epochUnix;

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

%irf_plot_ylabels_align(h), 
irf_zoom(h,'x',epochE([1 end])')


%% Run
mms_sdc_sdp_proc('ql', DCE_File,  DCV_File, HK_10E_File, HK_101_File);
mms_sdc_sdp_proc('scpot', DCE_File,  DCV_File, HK_10E_File, HK_101_File);

%% L2Pre
tt = irf_time(tint.start.utc,'utc>doy');
DEFATT_File = [data_root filesep 'ancillary' filesep mmsId filesep 'defatt'...
  filesep 'MMS' mmsId(end) '_DEFATT_' ...
  sprintf('%d%d_%d%d',tt(1),tt(2)-1,tt(1),tt(2)) '.V00'];
mms_sdc_sdp_proc('l2pre', DCE_File,  DCV_File, HK_10E_File, DEFATT_File);

%% L2a
mms_sdc_sdp_proc('l2a','out/mms4_edp_comm_l2pre_dce2d_20150506120000_v0.1.0.cdf')

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
B = mms.db_get_ts([mmsId '_dfg_srvy_ql'],[mmsId '_dfg_srvy_gsm_dmpa'],tint);

%% Plot with B

figure(75), clf
h = irf_plot(4,'newfigure');

hca = irf_panel('B');
hTmp = irf_plot(hca,B);
hTmp(1).Color=[0 0 0]; hTmp(2).Color=[0 0.5 0]; hTmp(3).Color=[1 0 0];
%hTmp(4).Color=[.5 .5 .5];
set(hca,'ColorOrder',[[0 0 0];[0 0.5 0];[1 0 0];[.5 .5 .5]])
irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.1])
ylabel(hca,'B [nT]')
title(hca,mmsId), set(hca,'YLim',64*[-1 1])

hca = irf_panel('Ex');
irf_plot(hca,{DceDSL.x,Es12.x,Es34.x},'comp')
ylabel(hca,'Ex DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('Ey');
irf_plot(hca,{DceDSL.y,Es12.y,Es34.y},'comp')
ylabel(hca,'Ey DSL [mV/m]'), set(hca,'YLim',E_YLIM*[-1 1])

hca = irf_panel('V');
irf_plot(hca,P2scPot)
ylabel(hca,'P2ScPot [V]'), set(hca,'YLim',[-14 0])

irf_plot_ylabels_align(h), irf_zoom(h,'x',DceDSL.time)
