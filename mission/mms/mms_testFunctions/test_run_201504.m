%% Init
data_root='/data/mms';
if ~exist(data_root,'dir')
  error('DATA_ROOT(%s) does not exist',data_root)
end
if ismac
  setenv('CDF_BASE','/Applications/cdf35_0-dist/')
  outDir = '/Users/yuri/Documents/MATLAB/MMS/test_proc/';
else, outDir = data_root;
end
cd(outDir)
if ~exist('log','dir'), mkdir('log'), end
if ~exist('out','dir'), mkdir('out'), end
setenv('DROPBOX_ROOT', [outDir filesep 'out'])
setenv('DATA_PATH_ROOT', [outDir filesep 'out'])
setenv('LOG_PATH_ROOT', '')
%setenv('LOG_PATH_ROOT', [outDir filesep 'log'])
global MMS_CONST;
MMS_CONST=mms_constants;
global ENVIR
ENVIR.CAL_PATH_ROOT = '/data/mms/irfu/cal/';
setenv('CAL_PATH_ROOT',ENVIR.CAL_PATH_ROOT)

%load /data/mms/irfu/mmsR.mat
%epocRTmp = EpochTT(R.time);

%% Define time
flagComm = 0;
%tint = irf.tint('2015-04-16T00:00:00Z/2015-04-16T06:00:00Z');
%tint = irf.tint('2015-04-16T18:00:00Z/2015-04-16T23:59:59Z');
%tint = irf.tint('2015-05-15T00:00:00Z/2015-05-15T05:59:59Z');
%tint = irf.tint('2015-04-20T18:00:00Z/2015-04-20T23:59:59Z');
%tint = irf.tint('2015-05-06T12:00:00Z/2015-05-06T17:59:59Z');
%tint = irf.tint('2015-06-21T00:00:00Z/2015-06-21T05:59:59Z'); 
%tint = irf.tint('2015-06-22T00:00:00Z/2015-06-22T23:59:59Z'); flagComm = false;
%tint = irf.tint('2015-08-15T13:00:00Z/2015-08-15T13:59:59Z'); flagComm = 2;
%tint = irf.tint('2015-09-11T09:30:00Z/2015-09-11T09:59:59Z'); flagComm = 2;
%tint = irf.tint('2015-10-07T11:00:00Z/2015-10-07T13:59:59Z'); flagComm = 2;
%tint = irf.tint('2015-10-16T05:02:34Z/2015-10-16T16:34:04Z'); flagComm = 2;
%tint = irf.tint('2015-10-22T04:17:34Z/2015-10-22T16:42:34Z'); flagComm = 2;
%tint = irf.tint('2015-12-18T00:00:00Z/2015-12-18T11:59:59Z'); flagComm = 2;
%tint = irf.tint('2015-12-21T11:55:00Z/2015-12-21T21:55:00Z'); flagComm = 3;
%tint = irf.tint('2015-12-23T00:00:00Z/2015-12-23T08:30:00Z'); flagComm = 2;
%tint = irf.tint('2015-11-16T01:25:04Z/2015-11-16T14:41:04Z'); flagComm = 2;
%tint = irf.tint('2015-11-29T00:01:34Z/2015-11-29T13:46:44Z'); flagComm = 2;
%tint = irf.tint('2015-11-13T01:44:14Z/2015-11-13T15:00:34Z'); flagComm = 2;
%tint = irf.tint('2015-12-14T00:00:00Z/2015-12-14T02:59:59Z'); flagComm = 2;
%tint = irf.tint('2016-01-03T00:00:00Z/2016-01-03T07:59:59Z'); flagComm = 2;
%tint = irf.tint('2016-01-31T00:00:00Z/2016-01-31T03:59:59Z'); flagComm = 2;
%tint = irf.tint('2016-01-03T23:00:00Z/2016-01-03T23:59:59Z'); flagComm = 2;
%tint = irf.tint('2016-02-25T20:00:00Z/2016-02-25T23:59:59Z'); flagComm = 2;
%tint = irf.tint('2016-06-12T05:00:00Z/2016-06-12T06:00:00Z'); flagComm = 2;
%tint = irf.tint('2016-09-17T00:00:00Z/2016-09-17T23:59:59Z'); flagComm = 2;
tint = irf.tint('2017-10-21T00:00:00Z/2017-10-21T23:59:59Z'); flagComm = 2;
mmsId = 'mms1'; 

prf = [data_root filesep mmsId]; utc = tint.start.toUtc(); 
mo = utc(6:7); yyyy=utc(1:4); day=utc(9:10); hh=utc(12:13); mm=utc(15:16);%
li = mms.db_list_files([mmsId '_fields_hk_l1b_101'],tint); if length(li)>1, error('li>1'), end
HK_101_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_fields_hk_l1b_105'],tint); if length(li)>1, error('li>1'), end
HK_105_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_fields_hk_l1b_10e'],tint); if length(li)>1, error('li>1'), end
HK_10E_File = [li.path filesep li.name];
li = mms.db_list_files([mmsId '_aspoc_srvy_l2'],tint); if length(li)>1, error('li>1'), end
ASPOC_File = [li.path filesep li.name];
if flagComm==1
  DCE_File  = [prf '/edp/comm/l1b/dce128/' yyyy '/' mo '/' mmsId ...
    '_edp_comm_l1b_dce128_' yyyy mo day hh mm '00_v0.8.0.cdf'];
  DCV_File  = [prf '/edp/comm/l1b/dcv128/' yyyy '/' mo '/' mmsId ...
    '_edp_comm_l1b_dcv128_' yyyy mo day hh mm '00_v0.8.0.cdf'];
elseif flagComm==2
  li = mms.db_list_files([mmsId '_edp_fast_l1b_dce'],tint); if length(li)>1, error('li>1'), end
  DCE_File = [li.path filesep li.name];
  DCV_File = [];
else
  li = mms.db_list_files([mmsId '_edp_slow_l1b_dce'],tint); if length(li)>1, error('li>1'), end
  DCE_File = [li.path filesep li.name];
  DCV_File = [];
end

%gseR = [epocRTmp.epochUnix R.(['gseR' mmsId(end)])];

%% Test QL - DMNGR
irf.log('log_out','screen'), irf.log('notice')
procId = MMS_CONST.SDCProc.ql; procName='QL'; scId=str2double(mmsId(end));
tmMode=MMS_CONST.TmMode.comm;
switch flagComm
  case 0, samplerate = MMS_CONST.Samplerate.slow;
  case 1, samplerate = MMS_CONST.Samplerate.comm_128;
  case 2, samplerate = MMS_CONST.Samplerate.fast;
  otherwise
    error('Bad flagComm')
end
% Initialize DMGR
Dmgr = mms_sdp_dmgr(scId,procId,tmMode,samplerate);
Dmgr.set_param('hk_10e',HK_10E_File);
Dmgr.set_param('hk_105',HK_105_File);
Dmgr.set_param('hk_101',HK_101_File);
Dmgr.set_param('aspoc',ASPOC_File);
list = mms.db_list_files([mmsId '_ancillary_defatt'],tint); 
for ii=1:length(list)
  irf.log('notice', [procName ' proc using: ',list(ii).name]);
  [dataTmp, src_fileData] = mms_load_ancillary([list(ii).path, filesep, ...
    list(ii).name], 'defatt');
  Dmgr.set_param('defatt', dataTmp);
end
Dmgr.set_param('dce',DCE_File);
if ~isempty(DCV_File),Dmgr.set_param('dcv',DCV_File); end
% Process
dce = Dmgr.dce;
probe2sc_pot = Dmgr.probe2sc_pot;
phase = Dmgr.phase;
spinfits = Dmgr.spinfits;
delta_off = Dmgr.delta_off;
adc_off = Dmgr.adc_off;
dce_xyz_dsl = Dmgr.dce_xyz_dsl;
dcv = Dmgr.dcv;

sampleRate=Dmgr.samplerate;
%phase_an

%%
spinEpoch = irf_spin_epoch( ...
  irf.ts_vec_xyz(dce_xyz_dsl.time, dce_xyz_dsl.data), ...
  phase, 'fCut', 1/40, 'nspins', 31, ...
  'samplefreq', samplerate);
dce_xyz_dsl.data = dce_xyz_dsl.data - spinEpoch.data;

%% Construct TSeries
DceSL = irf.ts_vec_xy(dce_xyz_dsl.time,[dce.e12.data dce.e34.data]);
DceDSL = irf.ts_vec_xyz(dce_xyz_dsl.time,dce_xyz_dsl.data);
Phase = irf.ts_scalar(dce_xyz_dsl.time,phase.data);
AdcOff12 = irf.ts_scalar(dce_xyz_dsl.time,adc_off.e12(:,1));
AdcOff34 = irf.ts_scalar(dce_xyz_dsl.time,adc_off.e34(:,1));
Es12 = irf.ts_vec_xy(spinfits.time,spinfits.sfit.e12(:,2:3));
Es34 = irf.ts_vec_xy(spinfits.time,spinfits.sfit.e34(:,2:3));
P2scPot = irf.ts_scalar(probe2sc_pot.time,probe2sc_pot.data);
Dcv = irf.ts_scalar(dcv.time,[dcv.v1.data dcv.v2.data dcv.v3.data dcv.v4.data]);

%% Summary plot
E_YLIM = 7;

h = irf_figure(73,5,'reset');

hca = irf_panel('E');
irf_plot(hca,DceSL)
hold(hca,'on') 
irf_plot(hca,{AdcOff12,AdcOff34},'comp')
ylabel(hca,'E SL [mV/m]')
title(hca,mmsId), set(hca,'YLim',49*[-1 1])

if 0
hca = irf_panel('Phase'); %#ok<UNRCH>
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

hca = irf_panel('Vs');
irf_plot(hca,Dcv)
ylabel(hca,'PPot [V]'), set(hca,'YLim',[-14 0])
irf_legend(hca,{'P1','P2','P3','P4'},[0.9,0.02])

%irf_plot_ylabels_align(h), 
irf_zoom(h,'x',DceDSL.time)
%add_position(h(end),gseR), xlabel(h(end),'')

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

%% L2a
tt = irf_time(tint.start.utc,'utc>doy');
DEFATT_File = [data_root filesep 'ancillary' filesep mmsId filesep 'defatt'...
  filesep 'MMS' mmsId(end) '_DEFATT_' ...
  sprintf('%d%d_%d%d',tt(1),tt(2)-1,tt(1),tt(2)) '.V00'];
mms_sdc_sdp_proc('l2a', DCE_File,  DCV_File, HK_10E_File);

%% L2pre
li = mms.db_list_files([mmsId '_dfg_srvy_l2pre'],tint); if length(li)>1, error('li>1'), end
DFG_File = [li.path filesep li.name];
mms_sdc_sdp_proc('l2pre','out/mms2_edp_slow_l2a_dce2d_20151221000000_v1.0.0.cdf',DFG_File)

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
B = mms.db_get_ts([mmsId '_afg_srvy_ql'],[mmsId '_afg_srvy_dmpa'],tint);

%% Plot with B
E_YLIM = 7;

figure(75), clf
h = irf_plot(4,'newfigure');

hca = irf_panel('B');
hTmp = irf_plot(hca,B);
hTmp(1).Color=[0 0 0]; hTmp(2).Color=[0 0.5 0]; hTmp(3).Color=[1 0 0];
hold(hca,'on'), hTmp = irf_plot(hca,B.abs());
hTmp.Color=[.5 .5 .5];
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
