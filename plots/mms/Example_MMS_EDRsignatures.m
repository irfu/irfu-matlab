% A routine to compute various parameters used to identify electron
% diffusion regions for the four MMS spacecraft. 
% Written by D. B. Graham
%
% Quantities calculated so far are:
% (1) sqrt(Q) - Based on Swisdak, 2015 (Appendix A) for any coordinate
% system. Supposed to be a better measure of agyrotropy than A phi_e. Not
% sure what value is expected in EDR's; lets say ~ 0.1 or higher. 
% (2) A phi_e = 2 abs(1 - alpha)/(1 + alpha), where alpha =
% Pe_perp1/Pe_perp2. This is a measure of electron agyrotropy. Values of 
% O(1) are expected for EDRs. So far for simplicity Pe_perp1 is the 
% component of the pressure tensor most aligned with the X direction and 
% perpendicular to B. (By definition does not include off-diagonal pressure
% terms.)
% (3) A n_e = T_parallel/T_perp. Values much larger than 1 are expected.
% Large T_parallel/T_perp are a feature of the ion diffusion region. For MP
% reconnection ion diffusion regions have A n_e ~ 3 based on MMS observations.
% Scudder says A n_e ~ 7 at IDR-EDR boundary, but this is extremely large
% for MP reconnection.
% (4) Mperp e - electron Mach number: bulk velocity divided by the electron
% thermal speed perpendicular to B. Values of O(1) are expected in EDRs. 
% (5) J.E - J.E > 0 is expected in the electron diffusion region,
% corresponding to dissipation of field energy. J is calculated on each
% spacecraft using the particle moments. 
% (6) epsilon_e - Energy gain per cyclotron period. Values of O(1) are
% expected in EDRs. 
% (7) delta_e - Relative strength of the electric and magnetic force in the
% bulk electron rest frame. N. B. Very sensitive to electron moments and
% electric field. Check version of these quantities. 
%
% Notes: kappa_e (not yet included) is taken to be the largest value of
% epsilon_e and delta_e at any given point. 
% Requires electron distributions with version number v1.0.0 or higher. 
% Calculations of agyrotropy measures (1) and (2) become unreliable at low
% densities n_e <~ 1 cm^-3, when the raw particle counts are low. 


%% Time interval selection
tint = irf.tint('2015-12-02T01:14:15.00Z/2015-12-02T01:15:03.50Z');

%% Load data
ic = 1:4;
tic;
c_eval('Bxyz?=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint);',ic);
c_eval('SCpot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_psp'',tint);',ic);
c_eval('SCpot?.data = -SCpot?.data*1.2+2;',ic);
c_eval('Exyzf? = mms.db_get_ts(''mms?_edp_fast_ql_dce'',''mms?_edp_dce_ql_dsl'',tint);',ic);
toc;

%% Load electron and ion particle data
% This way is fastest. Change directory to appropriate cdf here. 
c_eval('tmpDataObj? = dataobj(''data/mms?_fpi_brst_l1b_des-dist_20151202011414_v1.1.0.cdf'');',ic);
c_eval('pdiste? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_brstSkyMap_dist''));',ic);
c_eval('energye0? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_energy0'');',ic);
c_eval('energye1? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_energy1'');',ic);
c_eval('phie? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_brstSkyMap_phi''));',ic);
c_eval('thetae? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_theta'');',ic);
c_eval('stepTablee? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_stepTable_parity''));',ic);

c_eval('tmpDataObj? = dataobj(''data/mms?_fpi_brst_l1b_dis-dist_20151202011414_v1.1.0.cdf'');',ic);
c_eval('pdisti? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_dist''));',ic);
c_eval('energyi0? = get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_energy0'');',ic);
c_eval('energyi1? = get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_energy1'');',ic);
c_eval('phii? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_phi''));',ic);
c_eval('thetai? = get_variable(tmpDataObj?,''mms?_dis_brstSkyMap_theta'');',ic);
c_eval('stepTablei? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_dis_stepTable_parity''));',ic);

% This way is too slow
%for ii=1:4;
%tic;
%    c_eval('pdiste?=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint);',ii);
%    c_eval('energye0?=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_energy0'',tint);',ii);
%    c_eval('energye1?=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_energy1'',tint);',ii);
%    c_eval('phie?=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_phi'',tint);',ii);
%    c_eval('thetae?=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_theta'',tint);',ii);
%    c_eval('stepTablee?=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_stepTable_parity'',tint);',ii);
%toc;
%end
%for ii=1:4;
%tic;
%    c_eval('pdisti?=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint);',ii);
%    c_eval('energyi0?=mms.db_get_variable(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_energy0'',tint);',ii);
%    c_eval('energyi1?=mms.db_get_variable(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_energy1'',tint);',ii);
%    c_eval('phii?=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_phi'',tint);',ii);
%    c_eval('thetai?=mms.db_get_variable(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_theta'',tint);',ii);
%    c_eval('stepTablei?=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_stepTable_parity'',tint);',ii);
%toc;
%end

%% Compute particle moments and rotate pressure and temperature tensors
c_eval('emoments? = mms.psd_moments(pdiste?,phie?,thetae?,stepTablee?,energye0?,energye1?,SCpot?,''electron'');',ic);
c_eval('imoments? = mms.psd_moments(pdisti?,phii?,thetai?,stepTablei?,energyi0?,energyi1?,SCpot?,''ion'');',ic);

c_eval('Pet? = emoments?.P_psd.tlim(tint);',ic);
c_eval('[PeXXp?,PeXYp?,PeXZp?,PeYYp?,PeYZp?,PeZZp?] = mms.rotate_tensor_fac(Pet?,Bxyz?);',ic);
c_eval('Tet? = emoments?.T_psd.tlim(tint);',ic);
c_eval('[TeXXp?,TeXYp?,TeXZp?,TeYYp?,TeYZp?,TeZZp?] = mms.rotate_tensor_fac(Tet?,Bxyz?);',ic);

%% Compute tests for EDR
% Compute Q from full pressure tensor
c_eval('I1? = Pet?.data(:,1)+Pet?.data(:,4)+Pet?.data(:,6);',ic);
c_eval('I2? = Pet?.data(:,1).*Pet?.data(:,4)+Pet?.data(:,1).*Pet?.data(:,6)+Pet?.data(:,4).*Pet?.data(:,6)-((Pet?.data(:,2)).^2+(Pet?.data(:,3)).^2+(Pet?.data(:,5)).^2);',ic);
c_eval('Q? = 1-4*I2?./((I1?-PeZZp?.data).*(I1?+3*PeZZp?.data));',ic);
c_eval('Q? = irf.ts_scalar(Pet?.time,sqrt(Q?));',ic);

% Compute electron Mach number
Units = irf_units; 
qe = Units.e;
me = Units.me;
c_eval('Ue? = irf.ts_scalar(emoments?.V_psd.time,emoments?.V_psd.abs.data);',ic);
c_eval('Ue? = Ue?.tlim(tint);',ic);
c_eval('Veperp? = sqrt((TeXXp?.data+TeYYp?.data)*qe/me);',ic);
c_eval('Me? = Ue?.data*1000./Veperp?;',ic);
c_eval('Me? = irf.ts_scalar(PeXXp?.time,Me?);',ic);

% Compute agyrotropy Aphi and temperature ratio An
c_eval('alpha? = PeXXp?.data./PeYYp?.data;',ic);
c_eval('agyro? = 2*abs(1-alpha?)./(1+alpha?);',ic);
c_eval('agyro? = TSeries(PeXXp?.time,agyro?);',ic);
c_eval('Temprat? = 2*PeZZp?.data./(PeXXp?.data+PeYYp?.data);',ic);
c_eval('Temprat? = TSeries(PeXXp?.time,Temprat?);',ic);

% Compute current density and J.E
c_eval('ne? = emoments?.n_psd.tlim(tint);',ic);
c_eval('Uevec? = emoments?.V_psd.tlim(tint);',ic);
c_eval('Uivec? = imoments?.V_psd.tlim(tint);',ic);
c_eval('Exyzf? = Exyzf?.resample(Uevec?);',ic);
c_eval('Jmoms? = irf.ts_vec_xyz(Uevec?.time,1e18*qe*[ne?.data ne?.data ne?.data].*(Uivec?.data-Uevec?.data));',ic); % Current density in nA m^-2
c_eval('EdotJ = dot(Exyzf?.data,Jmoms?.data,2)/1000;',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ = irf.ts_scalar(ne?.time.time,EdotJ);',ic);

% Calculate epsilon and delta parameters
c_eval('Bxyzf? = Bxyz?.resample(Uevec?);',ic);
c_eval('Bmagf? = Bxyzf?.abs.data;',ic);
c_eval('omegace? = qe*Bmagf?/me*1e-9;',ic);
c_eval('EdotUe? = dot(Exyzf?.data,Uevec?.data,2)/1000;',ic);
c_eval('epsilone? = abs(6*pi*EdotUe?./(omegace?.*(TeXXp?.data+TeYYp?.data+TeZZp?.data)));',ic);
c_eval('epsilone? = irf.ts_scalar(Uevec?.time,epsilone?);',ic);

c_eval('UexB? = cross(Uevec?,Bxyzf?);',ic);
c_eval('UexB?.data = UexB?.data*1e-3;',ic);
c_eval('UexB?.data = Exyzf?.data+UexB?.data*1e-3;',ic);
c_eval('UexBabs? = UexB?.abs.data;',ic);

c_eval('deltae? = 1e-3*UexBabs?./(Veperp?.*Bmagf?*1e-9);',ic);
c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

%% Plot figure

mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
nplots = 8;
h=irf_plot(nplots,'newfigure');
xSize=750; ySize=750;
set(gcf,'Position',[10 10 xSize ySize]);

h(1)=irf_panel('Bz'); set(h(1),'ColorOrder',mmsColors)
irf_pl_tx(h(1),'Bxyz?',3);
ylabel(h(1),'B_{z} (nT)','Interpreter','tex');
irf_legend(h(1),{'MMS1','MMS2','MMS3','MMS4'},[0.9 0.98],'color','cluster')
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('sqrtQ'); set(h(2),'ColorOrder',mmsColors)
irf_pl_tx(h(2),'Q?',1);
ylabel(h(2),'$$\sqrt{Q}$$','Interpreter','latex');
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('agyro'); set(h(3),'ColorOrder',mmsColors)
irf_pl_tx(h(3),'agyro?',1);
ylabel(h(3),'A\phi_{e}','Interpreter','tex');
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('temprat'); set(h(4),'ColorOrder',mmsColors)
irf_pl_tx(h(4),'Temprat?',1);
ylabel(h(4),'An_{e}','Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

h(5)=irf_panel('Me'); set(h(5),'ColorOrder',mmsColors)
irf_pl_tx(h(5),'Me?',1);
ylabel(h(5),'M_{e \perp}','Interpreter','tex');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k')

h(6)=irf_panel('JdotE'); set(h(6),'ColorOrder',mmsColors)
irf_pl_tx(h(6),'Me?',1);
ylabel(h(6),'E.J (nW m^{-3})','Interpreter','tex');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k')

h(7)=irf_panel('epsilone'); set(h(7),'ColorOrder',mmsColors)
irf_pl_tx(h(7),'epsilone?',1);
ylabel(h(7),'\epsilon_{e}','Interpreter','tex');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k')

h(8)=irf_panel('deltae'); set(h(8),'ColorOrder',mmsColors)
irf_pl_tx(h(8),'deltae?',1);
ylabel(h(8),'\delta_{e}','Interpreter','tex');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k')

irf_plot_axis_align(h(1:nplots));
irf_zoom(h(1:nplots),'x',tint);

%% save figure

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','EDRsignatures2.png');