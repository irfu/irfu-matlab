function caa_pl_efw_pea_hia(cl_id)
%CAA_PL_EFW_PEA_HIA  compare E from EFW, PEACE and CIS_HIA
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%% EFW
efw = my_load(cl_id,'C?_CP_EFW_L3_E');
E_Vec_xy_ISR2 = getmat(efw, irf_ssub('E_Vec_xy_ISR2__C?_CP_EFW_L3_E',cl_id) );
tint = [E_Vec_xy_ISR2(1,1) E_Vec_xy_ISR2(end,1)];
efwp = my_load(cl_id,'C?_CP_EFW_L3_P');
ScPot = getmat(efwp, irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',cl_id) );

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

%% FGM
fgm = my_load(cl_id,'C?_CP_FGM_SPIN');
B_vec_xyz_gse = getmat(fgm, irf_ssub('B_vec_xyz_gse__C?_CP_FGM_SPIN',cl_id) );
B_vec_xyz_ISR2 = c_gse2dsi(B_vec_xyz_gse,SAX);

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

V_PEA_xyz_gse = getmat(pea, irf_ssub('Data_Velocity_GSE__C?_CP_PEA_MOMENTS',cl_id) );
V_PEA_xyz_ISR2 = c_gse2dsi(V_PEA_xyz_gse,SAX);
EVXB_PEA_xyz_ISR2 = irf_tappl(irf_cross(V_PEA_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

%% CIS-HIA
cis_hia = my_load(cl_id,'C?_PP_CIS');
V_HIA_xyz_gse = getmat(cis_hia, irf_ssub('V_HIA_xyz_gse__C?_PP_CIS',cl_id) );
V_HIA_xyz_ISR2 = c_gse2dsi(V_HIA_xyz_gse,SAX);
EVXB_HIA_xyz_ISR2 = irf_tappl(irf_cross(V_HIA_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

%% Computation

STEP = 60;
dt = ceil(range(tint)/STEP)*STEP;
t = tint(1):STEP:tint(1)+dt;
t = t';

E_Vec_xy_ISR2_rPEA = irf_resamp(E_Vec_xy_ISR2, EVXB_PEA_xyz_ISR2(:,1));
E_Vec_xy_ISR2_rHIA = irf_resamp(E_Vec_xy_ISR2, EVXB_HIA_xyz_ISR2(:,1));

diff_HIA = EVXB_HIA_xyz_ISR2(:,1:3);
diff_HIA(:,2:3) = EVXB_HIA_xyz_ISR2(:,2:3) - E_Vec_xy_ISR2_rHIA(:,2:3);
diff_HIA(:,2:3) = diff_HIA(:,2:3)*10;

diff_HIAr = irf_resamp(diff_HIA(~isnan(diff_HIA(:,2)),:),t);
for comp=1:2, diff_HIAr(abs(diff_HIAr(:,comp+1))<1,comp+1) = 1; end
diff_HIAr(:,2:3) = log10(abs(diff_HIAr(:,2:3))) - 1;

diff_PEA = EVXB_PEA_xyz_ISR2(:,1:3);
diff_PEA(:,2:3) = EVXB_PEA_xyz_ISR2(:,2:3) - E_Vec_xy_ISR2_rPEA(:,2:3);
diff_PEA(:,2:3) = diff_PEA(:,2:3)*10;

diff_PEAr = irf_resamp(diff_PEA(~isnan(diff_PEA(:,2)),:),t);
for comp=1:2, diff_PEAr(abs(diff_PEAr(:,comp+1))<1,comp+1) = 1; end
diff_PEAr(:,2:3) = log10(abs(diff_PEAr(:,2:3))) - 1;

%% Plotting

figure(111), clf

NPLOTS = 5;
DY = .05;
h=1:NPLOTS;

for comp=1:2
	h(comp) = irf_subplot(NPLOTS,1,-comp);
	irf_plot({E_Vec_xy_ISR2(:,[1 (comp+1)]),EVXB_PEA_xyz_ISR2(:,[1 (comp+1)]),...
		EVXB_HIA_xyz_ISR2(:,[1 (comp+1)])},'comp')
	r = range(E_Vec_xy_ISR2(:,comp+1));
	set(h(comp),'YLim',...
		[min(E_Vec_xy_ISR2(:,comp+1))-DY*r max(E_Vec_xy_ISR2(:,comp+1))+DY*r])
end

ylabel(h(1),'Ex [mV/m]')
ylabel(h(2),'Ey [mV/m]')
legend(h(1),'EFW','PEA','HIA')
legend(h(1),'boxoff')
title(h(1),irf_ssub('Cluster ?',cl_id))

comp_s='xy';
for comp=1:2
	h(comp+2) = irf_subplot(NPLOTS,1,-2-comp);
	irf_plot(diff_HIAr(:,[1 comp+1]))
	hold on
	irf_plot(diff_PEAr(:,[1 comp+1]),'g*')
	hold off
	ylabel(h(comp+2),['log (\Delta E' comp_s(comp) ') [mV/m]'])
	set(h(comp+2),'YLim',[-1 2.9])
end

h(NPLOTS) = irf_subplot(NPLOTS,1,-NPLOTS);
irf_plot(ScPot)
set(h(NPLOTS),'YColor','b')
ylabel('ScPot [-V]')

ts = t_start_epoch(tint(1));
for pl=1:NPLOTS
	set(h(pl),'XLim',tint - ts);
end

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
function dobj = my_load(cl_id,prod)

old_pwd = pwd;
d_s = irf_ssub(prod,cl_id);
if ~exist(d_s,'dir'), error([d_s ' : no such directory']), end

disp(['loading ' d_s]);
try
	cd(d_s)
	dobj = dataobj('*.cdf');
catch
	disp(['error loading ' d_s]);
	dobj = [];
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
