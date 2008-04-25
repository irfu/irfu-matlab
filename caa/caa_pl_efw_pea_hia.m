function caa_pl_efw_pea_hia(cl_id,vars)
%CAA_PL_EFW_PEA_HIA  compare E from EFW, PEACE and CIS_HIA
%
% CAA_PL_EFW_PEA_HIA(CL_ID,[VARS])
%
% VARS: istrument names separated by '|'
%
% Example:
%     caa_pl_efw_pea_hia(1,'edi|pea|cod')
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

flag_edi = 0;
flag_pea = 0;
flag_hia = 0;
flag_cod = 0;

if nargin==1
	vars = 'edi|hia|cod|pea';
end

toks = tokenize(vars,'|');
if isempty(toks), error('No variables specified to compare with'), end
for i=1:length(toks)
	switch lower(toks{i})
		case {'cod','codif'}
			flag_cod = 1;
		case 'edi'
			flag_edi = 1;
		case 'hia'
			flag_hia = 1;
		case {'pea','peace'}
			flag_pea = 1;
		otherwise
			disp(['unknown variable : ' toks{i}])
	end
end


%% EFW
efw = my_load(cl_id,'C?_CP_EFW_L3_E');
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
B_vec_xyz_gse = getmat(fgm, irf_ssub('B_vec_xyz_gse__C?_CP_FGM_SPIN',cl_id) );
B_vec_xyz_ISR2 = c_gse2dsi(B_vec_xyz_gse,SAX);

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

		diff_EDIr = get_diff_resamp(E_Vec_xy_ISR2, EDI_Vec_xyz_ISR2, t);
	catch
		flag_edi = 0;
		disp('cannot load EDI')
	end
end

%% PEA
if flag_pea
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

	diff_PEAr = get_diff_resamp(E_Vec_xy_ISR2, EVXB_PEA_xyz_ISR2, t);
end

%% CIS-HIA
if flag_hia
	try
		cis_hia = my_load(cl_id,'C?_PP_CIS');
		V_HIA_xyz_gse = getmat(cis_hia, irf_ssub('V_HIA_xyz_gse__C?_PP_CIS',cl_id) );
		V_HIA_xyz_ISR2 = c_gse2dsi(V_HIA_xyz_gse,SAX);
		EVXB_HIA_xyz_ISR2 = irf_tappl(irf_cross(V_HIA_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

		diff_HIAr = get_diff_resamp(E_Vec_xy_ISR2, EVXB_HIA_xyz_ISR2, t);
	catch
		flag_edi = 0;
		disp('cannot load HIA')
	end
end

%% CIS-HIA
if flag_cod
	cis_codif = my_load(cl_id,'C?_CP_CIS-CODIF_HS_H1_MOMENTS');
	V_COD_xyz_gse = getmat(cis_codif, irf_ssub('velocity__C?_CP_CIS-CODIF_HS_H1_MOMENTS',cl_id) );
	V_COD_xyz_ISR2 = c_gse2dsi(V_COD_xyz_gse,SAX);
	EVXB_COD_xyz_ISR2 = irf_tappl(irf_cross(V_COD_xyz_ISR2,B_vec_xyz_ISR2),'*(-1e-3)');

	diff_CODr = get_diff_resamp(E_Vec_xy_ISR2, EVXB_COD_xyz_ISR2, t);
end


%% Plotting

figure(111), clf

NPLOTS = 5;
DY = .05;
h=1:NPLOTS;

for comp=1:2
	h(comp) = irf_subplot(NPLOTS,1,-comp);
	irf_plot(E_Vec_xy_ISR2(:,[1 (comp+1)]))
	r = range(E_Vec_xy_ISR2(:,comp+1));
	set(h(comp),'YLim',...
		[min(E_Vec_xy_ISR2(:,comp+1))-DY*r max(E_Vec_xy_ISR2(:,comp+1))+DY*r])
	leg = {'EFW'};
	hold on
	if flag_pea
		irf_plot(EVXB_PEA_xyz_ISR2(:,[1 (comp+1)]),'r')
		leg = {leg{:} 'PEA'}; 
	end
	if flag_hia
		irf_plot(EVXB_HIA_xyz_ISR2(:,[1 (comp+1)]),'g')
		leg = {leg{:} 'HIA'};
	end
 	if flag_cod
 		irf_plot(EVXB_COD_xyz_ISR2(:,[1 (comp+1)]),'m')
 		leg = {leg{:} 'COD'};
 	end
	if flag_edi
		irf_plot(EDI_Vec_xyz_ISR2(:,[1 (comp+1)]),'k.')
		leg = {leg{:} 'EDI'};
	end
	hold off	
end

ylabel(h(1),'Ex [mV/m]')
ylabel(h(2),'Ey [mV/m]')
legend(h(1),leg)
legend(h(1),'boxoff')
title(h(1),irf_ssub('Cluster ? (position GSE)',cl_id))

comp_s='xy';
for comp=1:2
	h(comp+2) = irf_subplot(NPLOTS,1,-2-comp);
	hold on	
	if flag_hia, irf_plot(diff_HIAr(:,[1 comp+1])), end
	if flag_cod, irf_plot(diff_CODr(:,[1 comp+1]),'m'), end
	if flag_pea, irf_plot(diff_PEAr(:,[1 comp+1]),'g*'), end
	if flag_edi, irf_plot(diff_EDIr(:,[1 comp+1]),'k.'), end
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

function res = get_diff_resamp(E_EFW,E_OTH,TREF)

E_EFW_rOTH = irf_resamp(E_EFW, E_OTH(:,1));

diff_OTH = E_OTH(:,1:3);
diff_OTH(:,2:3) = E_OTH(:,2:3) - E_EFW_rOTH(:,2:3);
diff_OTH(:,2:3) = diff_OTH(:,2:3)*10;

res = irf_resamp(diff_OTH(~isnan(diff_OTH(:,2)),:),TREF);
for comp=1:2, res(abs(res(:,comp+1))<1,comp+1) = 1; end
res(:,2:3) = log10(abs(res(:,2:3))) - 1;
