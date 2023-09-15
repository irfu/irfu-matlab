Tint = irf.tint('2015-08-12T00:00:00Z/2015-08-12T23:59:59Z');

%% Load SDP E & B
if 0
  c_eval('B? = mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',Tint);') %#ok<UNRCH>
  c_eval('E? = mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',Tint);')
end

%% Load EDI
%c_eval('do=dataobj(''mms?_edi_srvy_ql_efield_20150801_v0.2.*.cdf''); Edi?=get_ts(do,''mms?_edi_E_dmpa''); EdiBc?=get_ts(do,''mms?_edi_E_bc_dmpa''); EdiQ?=get_ts(do,''mms?_edi_quality_bc''); clear do')
% mms?_edi_E_dmpa is from Bestarg and has a bug now. Maybe OK in next
% version
%c_eval('Edi?=mms.db_get_ts(''mms?_edi_srvy_ql_efield'',''mms?_edi_E_dmpa'',Tint);')
c_eval('EdiBc?=mms.db_get_ts(''mms?_edi_srvy_ql_efield'',''mms?_edi_E_bc_dmpa'',Tint);')
c_eval('EdiQ?=mms.db_get_ts(''mms?_edi_srvy_ql_efield'',''mms?_edi_quality_bc'',Tint);')
c_eval('if ~isempty(EdiQ?), EdiBc?.data(EdiQ?.data<2,:) = NaN; end')

%% PSP
c_eval('Pfast? = mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_psp'',Tint);')
c_eval('Pslow? = mms.db_get_ts(''mms?_edp_slow_l2_scpot'',''mms?_edp_psp'',Tint);')

Tint1 = irf.tint('2015-08-01T00:10:00Z/2015-08-01T05:10:00Z');
% Trich to avoid lines diring gaps in slow data
c_eval('ii=find(diff(Pfast?.time.epochUnix)>100); if ~isempty(ii), Pfast?.data(ii,:) = NaN; end')
c_eval('ii=find(diff(Pslow1.time.epochUnix)>100); if ~isempty(ii), Pslow?.data(ii,:) = NaN; end')
clear ii

%% ASPOC
c_eval('AI?=mms.db_get_ts(''mms?_aspoc_srvy_l2'',''mms?_asp_ionc'',Tint);')
% incorrect fillval setting
c_eval('AI?.data(AI?.data==min(AI?.data))=NaN;')

%%
if 0
  h = irf_plot(2); %#ok<UNRCH>

  hca = irf_panel('Ex');
  hl = irf_plot(hca,{E1.x,Edi1.x},'comp');
  hl.Children(1).Marker = '.';

  hca = irf_panel('Ey');
  hl = irf_plot(hca,{E1.y,Edi1.y},'comp');
  hl.Children(1).Marker = '.';
end

%% l2pre - fast
c_eval(['tsTmp=mms.db_get_ts(''mms?_edp_fast_l2pre_dce2d'',''mms?_edp_dce_spinfit_e12'',Tint);Es12fast?=irf.ts_vec_xy(tsTmp.time,tsTmp.data(:,3:4));'...
  'tsTmp=mms.db_get_ts(''mms?_edp_fast_l2pre_dce2d'',''mms?_edp_dce_spinfit_e34'',Tint);Es34fast?=irf.ts_vec_xy(tsTmp.time,tsTmp.data(:,3:4));']);
%% l2pre - slow
c_eval(['tsTmp=mms.db_get_ts(''mms?_edp_slow_l2pre_dce2d'',''mms?_edp_dce_spinfit_e12'',Tint);Es12slow?=irf.ts_vec_xy(tsTmp.time,tsTmp.data(:,3:4));'...
  'tsTmp=mms.db_get_ts(''mms?_edp_slow_l2pre_dce2d'',''mms?_edp_dce_spinfit_e34'',Tint);Es34slow?=irf.ts_vec_xy(tsTmp.time,tsTmp.data(:,3:4));']);

%%
h = irf_plot(9,'newfigure');

irf_plot(h(1),{Es12fast1.x,Es34fast1.x,Es12slow1.x,Es34slow1.x},'comp')
if ~isempty(EdiBc1), hold(h(1),'on'), hl = irf_plot(h(1),EdiBc1.x,'m.'); end
irf_plot(h(2),{Es12fast2.x,Es34fast2.x,Es12slow2.x,Es34slow2.x},'comp')
if ~isempty(EdiBc2), hold(h(2),'on'), hl = irf_plot(h(2),EdiBc2.x,'m.'); end
irf_plot(h(3),{Es12fast3.x,Es34fast3.x,Es12slow3.x,Es34slow3.x},'comp')
if ~isempty(EdiBc3), hold(h(3),'on'), hl = irf_plot(h(3),EdiBc3.x,'m.'); end
irf_plot(h(4),{Es12fast4.x,Es34fast4.x,Es12slow4.x,Es34slow4.x},'comp')
if ~isempty(EdiBc4), hold(h(4),'on'), hl = irf_plot(h(4),EdiBc4.x,'m.'); end

irf_plot(h(5),{Es12fast1.y,Es34fast1.y,Es12slow1.y,Es34slow1.y},'comp')
if ~isempty(EdiBc1), hold(h(5),'on'), hl = irf_plot(h(5),EdiBc1.y,'m.'); end
irf_plot(h(6),{Es12fast2.y,Es34fast2.y,Es12slow2.y,Es34slow2.y},'comp')
if ~isempty(EdiBc2), hold(h(6),'on'), hl = irf_plot(h(6),EdiBc2.y,'m.'); end
irf_plot(h(7),{Es12fast3.y,Es34fast3.y,Es12slow3.y,Es34slow3.y},'comp')
if ~isempty(EdiBc3), hold(h(7),'on'), hl = irf_plot(h(7),EdiBc3.y,'m.'); end
irf_plot(h(8),{Es12fast4.y,Es34fast4.y,Es12slow4.y,Es34slow4.y},'comp')
if ~isempty(EdiBc4), hold(h(8),'on'), hl = irf_plot(h(8),EdiBc4.y,'m.'); end

irf_pl_tx(h(9),'Pfast?'); hold(h(9),'on'), irf_pl_tx(h(9),'Pslow?');

ylabel(h(1),'Ex sc1')
ylabel(h(2),'Ex sc2')
ylabel(h(3),'Ex sc3')
ylabel(h(4),'Ex sc4')

ylabel(h(5),'Ey sc1')
ylabel(h(6),'Ey sc2')
ylabel(h(7),'Ey sc3')
ylabel(h(8),'Ey sc4')

ylabel(h(9),'Probe2Sc [V]')

hl = legend(h(1),'12f','34f','12s','34s','EDI');
hl.Box = 'off';

%%
h = irf_plot(7,'newfigure');

hca = h(1);
irf_pl_tx(hca,'Es12fast?.x'), hold(hca,'on')
irf_pl_tx(hca,'Es12slow?.x')

hca = h(2);
irf_pl_tx(hca,'Es34fast?.x'), hold(hca,'on')
irf_pl_tx(hca,'Es34slow?.x')

hca = h(3);
irf_pl_tx(hca,'Es12fast?.y'), hold(hca,'on')
irf_pl_tx(hca,'Es12slow?.y')

hca = h(4);
irf_pl_tx(hca,'Es34fast?.y'), hold(hca,'on')
irf_pl_tx(hca,'Es34slow?.y')

c_eval('dEfast? = Es12fast? - Es34fast?; dEslow? = Es12slow? - Es34slow?;')

hca = h(5);
irf_pl_tx(hca,'abs(dEfast?.x)'), hold(hca,'on')
irf_pl_tx(hca,'abs(dEslow?.x)')

hca = h(6);
irf_pl_tx(hca,'abs(dEfast?.y)'), hold(hca,'on')
irf_pl_tx(hca,'abs(dEslow?.y)')

hca = h(7);
irf_pl_tx(hca,'Pfast?'); hold(hca,'on'), irf_pl_tx(hca,'Pslow?');

ylabel(h(1),'E12_x [mV/m]')
ylabel(h(2),'E34_x [mV/m]')
ylabel(h(3),'E12_y [mV/m]')
ylabel(h(4),'E34_y [mV/m]')
ylabel(h(5),'|E12_x-E34_x|')
ylabel(h(6),'|E12_y-E34_y|')
ylabel(h(7),'Probe2Sc [V]')
