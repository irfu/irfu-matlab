function caa_pl_efw_pea_hia(cl_id,vars,plot_range,plot_vz)
%CAA_PL_EFW_PEA_HIA  compare E from EFW, PEACE and CIS_HIA
%
% CAA_PL_EFW_PEA_HIA(CL_ID,[VARS,PLOT_RANGE,PLOT_VZ])
%
% VARS : instrument names (edi,pea,hia,cod) separated by '|'
%
% PLOT_RANGE : 0 - line plot, 1 - min and max averages, 2 - fill area
% between min and max. Can be set as one number for all VARS of for each
% variable separately
%
% PLOT_VZ : 0 - turn off plotting of Vz
%
% Example:
%     caa_pl_efw_pea_hia(1,'edi|pea|cod','022')
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin < 4 || plot_vz ~=0
	plot_vz = 1;
end

flag_edi = 0;
flag_pea = 0;
flag_hia = 0;
flag_cod = 0;

plot_range_edi = 0;
plot_range_hia = 0;
plot_range_pea = 0;
plot_range_cod = 0;

if nargin==1
	vars = 'edi|hia|cod|pea';
end

if nargin < 3, plot_range = 0; end

toks = tokenize(vars,'|');
if isempty(toks), error('No variables specified to compare with'), end
for i=1:length(toks)
	switch lower(toks{i})
		case {'cod','codif'}
			flag_cod = 1;
			plot_range_cod = get_plot_range(plot_range,i);
		case 'edi'
			flag_edi = 1;
			plot_range_edi = get_plot_range(plot_range,i);
		case 'hia'
			flag_hia = 1;
			plot_range_hia = get_plot_range(plot_range,i);
		case {'pea','peace'}
			flag_pea = 1;
			plot_range_pea = get_plot_range(plot_range,i);
		otherwise
			disp(['unknown variable : ' toks{i}])
	end
end


%% EFW
efw = my_load(cl_id,'C?_CP_EFW_L3_E');
if isempty(efw), return, end
E_Vec_xy_ISR2 = getmat(efw, irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L3_E',cl_id) );
tint = [E_Vec_xy_ISR2(1,1) E_Vec_xy_ISR2(end,1)];
efwp = my_load(cl_id,'C?_CP_EFW_L3_P');
ScPot = getmat(efwp, irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',cl_id) );

%% New time
STEP = 60;
dt = ceil(range(tint)/STEP)*STEP;
t = tint(1):STEP:tint(1)+dt;
t = t';

%% SAX
[ok,SAX] = c_load('SAX?',cl_id);
if ~ok
	getData(ClusterDB,tint(1),range(tint),cl_id,'sax')
	[ok,SAX] = c_load('SAX?',cl_id);
	if ~ok
		error('cannot load SAX')
	end
end
clear ok

%% R
[ok,R] = c_load('R?',cl_id);
if ~ok
	getData(ClusterDB,tint(1),range(tint),cl_id,'r')
	[ok,R] = c_load('R?',cl_id);
	if ~ok
		disp('cannot load R')
		R = [];
	end
end
clear ok

%% V
if flag_edi
	[ok,V] = c_load('V?',cl_id);
	if ~ok
		getData(ClusterDB,tint(1),range(tint),cl_id,'v')
		[ok,V] = c_load('V?',cl_id);
		if ~ok
			error('cannot load V')
		end
	end
	clear ok
end

%% FGM
fgm = my_load(cl_id,'C?_CP_FGM_SPIN');
if isempty(fgm), return, end
B_vec_xyz_gse = getmat(fgm, irf_ssub('B_vec_xyz_gse__C?_CP_FGM_SPIN',cl_id) );
B_vec_xyz_ISR2 = c_gse2dsi(B_vec_xyz_gse,SAX);

ttt = irf_e_vxb([E_Vec_xy_ISR2 E_Vec_xy_ISR2(:,1)*0],B_vec_xyz_ISR2,-1);
E_Vec_xy_ISR2(:,4) = ttt(:,4); % Vz in 3rd column
clear ttt

%% EDI
if flag_edi
	try
		edi = my_load(cl_id,'C?_CP_EDI_MP');
		iEDI_Vec_xyz_gse = getmat(edi,irf_ssub('E_xyz_gse__C?_CP_EDI_MP',cl_id) );

		B_EDI = irf_resamp(B_vec_xyz_gse,iEDI_Vec_xyz_gse);
		evxb = irf_tappl(irf_cross(B_EDI,irf_resamp(V,B_EDI)),'*1e-3*(-1)');

		EDI_Vec_xyz_gse = iEDI_Vec_xyz_gse;
		EDI_Vec_xyz_gse(:,2:4) = EDI_Vec_xyz_gse(:,2:4) + evxb(:,2:4); %#ok<NASGU>
		clear evxb

		EDI_Vec_xyz_ISR2 = c_gse2dsi(EDI_Vec_xyz_gse,SAX);
		
		% Vz in 3rd column
		ttt = irf_e_vxb(EDI_Vec_xyz_ISR2,B_vec_xyz_ISR2,-1);
		EDI_Vec_xyz_ISR2(:,4) = ttt(:,4); 
		clear ttt
		
		[diff_EDIr,E_Vec_xy_ISR2_rEDI] = get_diff_resamp(E_Vec_xy_ISR2, EDI_Vec_xyz_ISR2, t);
	catch
		flag_edi = 0;
		disp('cannot load EDI')
	end
end

%% PEA
pea = my_load(cl_id,'C?_CP_PEA_MOMENTS');
T_PEA_PAR = getmat(pea, ...
	irf_ssub('Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS',cl_id) );
T_PEA_units = getunits(pea, ...
	irf_ssub('Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS',cl_id) );
try
	T_PEA_PERP = getmat(pea, ...
		irf_ssub('Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_',cl_id) );
catch
	disp('trying alternative for Te_perp')
	try
		T_PEA_PERP = getmat(pea, ...
			irf_ssub('Data_Temperature_ComponentPerpendicularToMagField__C?__MOMENTS',cl_id) );
	catch
		disp('no luck here as well')
		T_PEA_PERP = [];
	end
end

if flag_pea
	V_PEA_xyz_gse = getmat(pea, irf_ssub('Data_Velocity_GSE__C?_CP_PEA_MOMENTS',cl_id) );
	V_PEA_xyz_ISR2 = c_gse2dsi(V_PEA_xyz_gse,SAX);
	EVXB_PEA_xyz_ISR2 = irf_tappl(irf_cross(V_PEA_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

	[apar,aperp] = irf_dec_parperp(B_vec_xyz_ISR2,V_PEA_xyz_ISR2); %#ok<ASGLU>
	EVXB_PEA_xyz_ISR2(:,4) = aperp(:,4);
	clear apar aperp
	
	[diff_PEAr,E_Vec_xy_ISR2_rPEA] = get_diff_resamp(E_Vec_xy_ISR2, EVXB_PEA_xyz_ISR2, t);
end

%% CIS-HIA
if flag_hia
	try
		cis_hia = my_load(cl_id,'C?_PP_CIS');
		V_HIA_xyz_gse = getmat(cis_hia, irf_ssub('V_HIA_xyz_gse__C?_PP_CIS',cl_id) );
		V_HIA_xyz_ISR2 = c_gse2dsi(V_HIA_xyz_gse,SAX);
		EVXB_HIA_xyz_ISR2 = irf_tappl(irf_cross(V_HIA_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

		[diff_HIAr,E_Vec_xy_ISR2_rHIA] = get_diff_resamp(E_Vec_xy_ISR2, EVXB_HIA_xyz_ISR2, t);
	catch
		flag_edi = 0;
		disp('cannot load HIA')
	end
end

%% CIS-CODIF
if flag_cod
	try
		cis_codif = my_load(cl_id,'C?_CP_CIS-CODIF_HS_H1_MOMENTS');
		V_COD_xyz_gse = getmat(cis_codif, irf_ssub('velocity__C?_CP_CIS-CODIF_HS_H1_MOMENTS',cl_id) );
		V_COD_xyz_ISR2 = c_gse2dsi(V_COD_xyz_gse,SAX);
		EVXB_COD_xyz_ISR2 = irf_tappl(irf_cross(V_COD_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

		[diff_CODr,E_Vec_xy_ISR2_rCOD] = get_diff_resamp(E_Vec_xy_ISR2, EVXB_COD_xyz_ISR2, t);
	catch
		flag_cod = 0;
		disp('cannot load COD')
	end		
end


%% Plotting

figure(111), clf
if plot_vz
	NPLOTS = 9;
	NCOMP = 3;
else
	NPLOTS = 7;
	NCOMP = 2;
end
DY = .05;
h=3:NPLOTS;

NVARS = 0;
if flag_pea, NVARS = NVARS + 1; end
if flag_hia, NVARS = NVARS + 1; end
if flag_cod, NVARS = NVARS + 1; end
if flag_edi, NVARS = NVARS + 1; end

if NVARS==0, error('No data to compare with'), end

hh = zeros(NVARS,2);
pl = 0;
for comp=1:2
	for cvar=1:NVARS
		pl = pl + 1;
		hh(cvar,comp) = irf_subplot(NPLOTS,NVARS,-pl);
	end
end

pl = 1;
if flag_edi
	for comp=1:2
		plot(hh(pl,comp),...
			E_Vec_xy_ISR2_rEDI(:,comp+1),EDI_Vec_xyz_ISR2(:,comp+1),'.')
	end
	title(hh(pl,1),'EDI')
	pl = pl + 1;
end
if flag_hia
	for comp=1:2
		plot(hh(pl,comp),...
			E_Vec_xy_ISR2_rHIA(:,comp+1),EVXB_HIA_xyz_ISR2(:,comp+1),'.')
	end
	title(hh(pl,1),'HIA')
	pl = pl + 1;
end
if flag_cod
	for comp=1:2
		plot(hh(pl,comp),...
			E_Vec_xy_ISR2_rCOD(:,comp+1),EVXB_COD_xyz_ISR2(:,comp+1),'.')
	end
	title(hh(pl,1),'COD')
	pl = pl + 1;
end
if flag_pea
	for comp=1:2
		plot(hh(pl,comp),...
			E_Vec_xy_ISR2_rPEA(:,comp+1),EVXB_PEA_xyz_ISR2(:,comp+1),'.')
	end
	title(hh(pl,1),'PEA')
end

set(hh(:,1),'XAxisLocation','top');
set(hh(:,2),'XTickLabel',[]);
ylabel(hh(1,1),'Ex [mV/m]')
ylabel(hh(1,2),'Ey [mV/m]')
r = max(range(E_Vec_xy_ISR2(:,2)),range(E_Vec_xy_ISR2(:,3)));
YLim = [min(min(E_Vec_xy_ISR2(:,2:3)))-DY*r max(max(E_Vec_xy_ISR2(:,2:3)))+DY*r];
r = range(E_Vec_xy_ISR2(:,4));
YLimVZ = [min(E_Vec_xy_ISR2(:,4))-DY*r max(E_Vec_xy_ISR2(:,4))+DY*r];
clear r
set(hh,'XLim',YLim,'YLim',YLim,'XGrid','on','YGrid','on')
%set(hh,'DataAspectRatioMode','manual') % Makes axes square

OFF = 2;
ts = [];
for comp=1:NCOMP
	h(OFF+comp) = irf_subplot(NPLOTS,1,-OFF-comp);
	hold on
	
	leg = {};
	if flag_pea
		if plot_range_pea==1
			IDX_ST_PEA = 4;
			irf_plot(get_mm_resamp('min',EVXB_PEA_xyz_ISR2(:,[1 (comp+1)]),...
				t(1:IDX_ST_PEA:fix(length(t)/IDX_ST_PEA*IDX_ST_PEA))),'g')
		elseif plot_range_pea==2
			IDX_ST_PEA = 8;
			plot_area(EVXB_PEA_xyz_ISR2(:,[1 (comp+1)]),...
				t(1:IDX_ST_PEA:fix(length(t)/IDX_ST_PEA*IDX_ST_PEA)),[.788 1 .708])
		else
			irf_plot(EVXB_PEA_xyz_ISR2(:,[1 (comp+1)]),'r')
		end
		leg = {leg{:} 'PEA'}; 
	end
	if flag_hia
		if plot_range_hia==1
			irf_plot(get_mm_resamp('min',EVXB_HIA_xyz_ISR2(:,[1 (comp+1)]),t),'r')
		elseif plot_range_hia==2
			plot_area(EVXB_HIA_xyz_ISR2(:,[1 (comp+1)]),t,[1 .619 .564])
		else
			irf_plot(EVXB_HIA_xyz_ISR2(:,[1 (comp+1)]),'g')
		end
		leg = {leg{:} 'HIA'};
	end
 	if flag_cod
		if plot_range_cod==1
			irf_plot(get_mm_resamp('min',EVXB_COD_xyz_ISR2(:,[1 (comp+1)]),t),'m')
		elseif plot_range_cod==2
			plot_area(EVXB_COD_xyz_ISR2(:,[1 (comp+1)]),t,[1 .913 .947])
		else 
			irf_plot(EVXB_COD_xyz_ISR2(:,[1 (comp+1)]),'m')
		end
 		leg = {leg{:} 'COD'};
 	end
	if flag_edi
		if plot_range_edi==1
			irf_plot(get_mm_resamp('min',EDI_Vec_xyz_ISR2(:,[1 (comp+1)]),t),'k')
		elseif plot_range_edi==2
			plot_area(EDI_Vec_xyz_ISR2(:,[1 (comp+1)]),t,[.8 .8 .8])
		else
			irf_plot(EDI_Vec_xyz_ISR2(:,[1 (comp+1)]),'k.')
		end
		leg = {leg{:} 'EDI'};
	end
	
	irf_plot(E_Vec_xy_ISR2(:,[1 (comp+1)]))
	leg = {leg{:} 'EFW'};
	
	if flag_pea && plot_range_pea==1
		irf_plot(get_mm_resamp('max',EVXB_PEA_xyz_ISR2(:,[1 (comp+1)]),...
			t(1:IDX_ST_PEA:fix(length(t)/IDX_ST_PEA*IDX_ST_PEA))),'r')
	end
	if flag_hia && plot_range_hia==1
		irf_plot(get_mm_resamp('max',EVXB_HIA_xyz_ISR2(:,[1 (comp+1)]),t),'g')
	end
	if flag_cod && plot_range_cod==1
		irf_plot(get_mm_resamp('max',EVXB_COD_xyz_ISR2(:,[1 (comp+1)]),t),'m')
	end
	if flag_edi && plot_range_edi==1
		irf_plot(get_mm_resamp('max',EDI_Vec_xyz_ISR2(:,[1 (comp+1)]),t),'k')
	end
	
	hold off
	
	if comp<=2, set(h(OFF+comp),'YLim',YLim)
	else set(h(OFF+comp),'YLim',YLimVZ)
	end
	if isempty(ts), ts = t_start_epoch(tint(1)); end
end

ylabel(h(OFF+1),'Ex [mV/m]')
ylabel(h(OFF+2),'Ey [mV/m]')
if plot_vz, ylabel(h(OFF+3),'Vz_{\perp} [km/s]'), end
legend(h(OFF+1),leg)
legend(h(OFF+1),'boxoff')

if plot_vz, OFF=5; else OFF = 4; end

comp_s='xyz';
comp_v='EEV';
for comp=1:NCOMP
	h(OFF+comp) = irf_subplot(NPLOTS,1,-OFF-comp);
	hold on
	leg = {};
	if flag_hia, irf_plot(diff_HIAr(:,[1 comp+1])), leg = {leg{:} 'HIA'}; end
	if flag_cod, irf_plot(diff_CODr(:,[1 comp+1]),'m'), leg = {leg{:} 'COD'}; end
	if flag_pea, irf_plot(diff_PEAr(:,[1 comp+1]),'g*'), leg = {leg{:} 'PEA'}; end
	if flag_edi, irf_plot(diff_EDIr(:,[1 comp+1]),'k.'), leg = {leg{:} 'EDI'}; end
	hold off
	ylabel(h(OFF+comp),['log (\Delta' comp_v(comp) comp_s(comp) ')'])
	if comp<=2, set(h(OFF+comp),'YLim',[-1 1.99],'Box','on')
	else set(h(OFF+comp),'YLim',[-1 2.99],'Box','on')
	end
	legend(h(OFF+comp),leg)
	legend(h(OFF+comp),'boxoff')
end

h(NPLOTS) = irf_subplot(NPLOTS,1,-NPLOTS);
irf_plot(ScPot)
set(h(NPLOTS),'YColor','b')
ylabel('ScPot [-V]')

for pl=3:NPLOTS
	set(h(pl),'XLim',tint - ts);
end
ts_s = epoch2iso(tint(1));
add_text(h(NPLOTS),sprintf('Cluster %d %s (position GSE)',cl_id,ts_s(1:10)))

if ~isempty(R), add_position(h(NPLOTS),R), end

ax2 = axes('Position',get(h(NPLOTS),'Position'),...
	'XAxisLocation','top',...
	'YAxisLocation','right',...
	'Color','none',...
	'XColor','k','YColor','r',...
	'XTickLabel',[]);
axes(ax2)
line(T_PEA_PAR(:,1)-ts,T_PEA_PAR(:,2),'Color','k','Marker','d','Parent',ax2);
if ~isempty(T_PEA_PERP)
	line(T_PEA_PERP(:,1)-ts,T_PEA_PERP(:,2),...
		'Color','r','Marker','+','Parent',ax2);
end
set(ax2,'XLim',tint - ts);
ylabel(['Te [' T_PEA_units ']'])

orient tall
  

%% Help function my_load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dobj = my_load(cl_id,prod)

old_pwd = pwd;
d_s = irf_ssub(prod,cl_id);
if ~exist(d_s,'dir')
	disp(['error loading ' d_s ' : no such directory']);
	dobj = [];
else
	disp(['loading ' d_s]);
	try
		cd(d_s)
		dobj = dataobj('*.cdf');
	catch
		disp(['error loading ' d_s ' : no CDF files']);
		dobj = [];
	end
end
cd(old_pwd)

%% Help function t_start_epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_st_e = t_start_epoch(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gives back the value of t_start_epoch of the figure
% if not  set, sets t_start_epoch of the figure
ud = get(gcf,'userdata');
ii = find(~isnan(t));
if ~isempty(ii), valid_time_stamp = t(ii(1)); else valid_time_stamp = []; end

if isfield(ud,'t_start_epoch')
	t_st_e = double(ud.t_start_epoch);
elseif ~isempty(valid_time_stamp)
	if valid_time_stamp > 1e8
		% Set start_epoch if time is in isdat epoch
		% Warn about changing t_start_epoch
		t_st_e = double(valid_time_stamp);
		ud.t_start_epoch = t_st_e;
		set(gcf,'userdata',ud);
		irf_log('proc',['user_data.t_start_epoch is set to ' ...
			epoch2iso(t_st_e,1)]);
	else
		t_st_e = double(0);
	end
else
	t_st_e = double(0);
end

%% Help function get_diff_resamp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res,E_EFW_rOTH] = get_diff_resamp(E_EFW,E_OTH,TREF)

MC = 4;

E_EFW_rOTH = irf_resamp(E_EFW, E_OTH(:,1));

diff_OTH = E_OTH(:,1:MC);
diff_OTH(:,2:MC) = E_OTH(:,2:MC) - E_EFW_rOTH(:,2:MC);
diff_OTH(:,2:MC) = diff_OTH(:,2:MC)*10;

res = irf_resamp(diff_OTH(~isnan(diff_OTH(:,2)),:),TREF);
for comp=1:(MC-1), res(abs(res(:,comp+1))<1,comp+1) = 1; end
res(:,2:MC) = log10(abs(res(:,2:MC))) - 1;

%% Help function get_mm_resamp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = get_mm_resamp(op,E,TREF)

res = zeros(length(TREF),size(E,2));
res(:,1) = TREF;
res(:,2:end) = NaN;

STEP2 = (TREF(2) - TREF(1))/2;

for i=1:length(TREF)
 	switch lower(op)
 		case 'min'
 			tmp = min(E( E(:,1)>=TREF(i)-STEP2 & E(:,1)<TREF(i)+STEP2 , 2:end));
		case 'max'
			tmp = max(E( E(:,1)>=TREF(i)-STEP2 & E(:,1)<TREF(i)+STEP2 , 2:end));
		otherwise
			error('unknown operation')
	end
	if ~isempty(tmp), res(i,2:end) = tmp; end
end
res = res(~isnan(res(:,2)),:);

%% Help function plot_area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_area(E,TREF,COL)

top = get_mm_resamp('max',E,TREF);
bot = get_mm_resamp('min',E,TREF);
data = [top; flipud(bot)];

ts = t_start_epoch(TREF(1));
data(:,1) = data(:,1) - ts;
hh = fill(data(:,1),data(:,2),COL);
set(hh,'EdgeColor',COL)

%% Help function add_text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function add_text(h,txt)

ylim = get(h,'YLim');
xlim = get(h,'XLim');
text(xlim(1) + range(xlim)*.5, ylim(2) - range(ylim)*.15, [' ' txt],...
        'FontWeight','bold','HorizontalAlignment','center')
	
%% Help function add_text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = get_plot_range(plot_range,i)

if ischar(plot_range)
	if length(plot_range)>1
		res = str2double(plot_range(i));
	else
		res = str2double(plot_range);
	end
elseif iscell(plot_range)
	res = plot_range{i};
elseif isnumeric(plot_range)
	if length(plot_range)>1
		res = plot_range(i);
	else
		res = plot_range;
	end
else
	error('bad format for PLOT_RANGE')
end

if (res~=0) && (res~=1) && (res~=2)
	error('bad value for PLOT_RANGE at position %d: %d',i,res)
end

