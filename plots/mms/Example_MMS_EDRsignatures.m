% A routine to compute various parameters used to identify electron
% diffusion regions for the four MMS spacecraft. 
% Written by D. B. Graham
%
% Quantities calculated so far are:
% (1) sqrt(Q) - Based on Swisdak, 2015 (arxiv). Not
% sure what value is expected in EDR's; lets say something above ~ 0.1. 
% (2) Dng - Based on Aunia et al., 2013; Computed based on the off-diagonal
% terms in the pressure tensor for Pe_perp1 = Pe_perp2. 
% (3) A phi_e = 2 abs(Perp1-Perp2)/(Perp1+Perp2). 
% This is a measure of electron agyrotropy. Values of 
% O(1) are expected for EDRs. We transform the pressure tensor into 
% field-aligned coordinates such that the difference in Pe_perp1 and Pe_perp2
% is maximal. This corresponds to P23 being zero. (Note that this definition 
% of agyrotropy neglects the off-diagonal pressure terms P12 and P13, 
% therefore it doesn't capture all agyrotropies.)
% (4) A n_e = T_parallel/T_perp. Values much larger than 1 are expected.
% Large T_parallel/T_perp are a feature of the ion diffusion region. For MP
% reconnection ion diffusion regions have A n_e ~ 3 based on MMS observations.
% Scudder says A n_e ~ 7 at IDR-EDR boundary, but this is extremely large
% for MP reconnection.
% (5) Mperp e - electron Mach number: bulk velocity divided by the electron
% thermal speed perpendicular to B. Values of O(1) are expected in EDRs. 
% (6) J.E - J.E > 0 is expected in the electron diffusion region,
% corresponding to dissipation of field energy. J is calculated on each
% spacecraft using the particle moments. 
% (7) epsilon_e - Energy gain per cyclotron period. Values of O(1) are
% expected in EDRs. 
% (8) delta_e - Relative strength of the electric and magnetic force in the
% bulk electron rest frame. N. B. Very sensitive to electron moments and
% electric field. Check version of these quantities. 
%
% Notes: kappa_e (not yet included) is taken to be the largest value of
% epsilon_e and delta_e at any given point. 
% Requires electron distributions with version number v1.0.0 or higher. 
% Calculations of agyrotropy measures (1)--(3) become unreliable at low
% densities n_e <~ 5 cm^-3, when the raw particle counts are low. 

%% Time interval selection
Tint = irf.tint('2015-12-02T01:14:15.00Z/2015-12-02T01:15:03.50Z');

%% Load data
ic = 1:4;
tic;
c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_dmpa_srvy_l2'',Tint);',ic);
c_eval('E? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',Tint);',ic);
toc;

tic;
for ii=1:4;
   c_eval('ne?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',Tint);',ii);
   c_eval('VeX?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkx_dbcs_brst'',Tint);',ii);
   c_eval('VeY?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulky_dbcs_brst'',Tint);',ii);
   c_eval('VeZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkz_dbcs_brst'',Tint);',ii);
   c_eval('TeXX?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempxx_dbcs_brst'',Tint);',ii);
   c_eval('TeXY?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempxy_dbcs_brst'',Tint);',ii);
   c_eval('TeXZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempxz_dbcs_brst'',Tint);',ii);
   c_eval('TeYY?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempyy_dbcs_brst'',Tint);',ii);
   c_eval('TeYZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempyz_dbcs_brst'',Tint);',ii);
   c_eval('TeZZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_tempzz_dbcs_brst'',Tint);',ii);
   c_eval('PeXX?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_presxx_dbcs_brst'',Tint);',ii);
   c_eval('PeXY?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_presxy_dbcs_brst'',Tint);',ii);
   c_eval('PeXZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_presxz_dbcs_brst'',Tint);',ii);
   c_eval('PeYY?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_presyy_dbcs_brst'',Tint);',ii);
   c_eval('PeYZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_presyz_dbcs_brst'',Tint);',ii);
   c_eval('PeZZ?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_preszz_dbcs_brst'',Tint);',ii);
end
toc;

c_eval('Uevec? = irf.ts_vec_xyz(VeX?.time,[VeX?.data VeY?.data VeZ?.data]);',ic);

for ii=1:4;
   c_eval('ViX?=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkx_dbcs_brst'',Tint);',ii);
   c_eval('ViY?=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulky_dbcs_brst'',Tint);',ii);
   c_eval('ViZ?=mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_bulkz_dbcs_brst'',Tint);',ii);
end

c_eval('Uivec? = irf.ts_vec_xyz(ViX?.time,[ViX?.data ViY?.data ViZ?.data]);',ic);
c_eval('E? = E?.resample(ne?);',ic);
c_eval('Bxyz? = Bxyz?.resample(ne?);',ic);
c_eval('Uivec? = Uivec?.resample(ne?);',ic);

%% Rotate pressure and temperature tensors
c_eval('Pepp? = mms.rotate_tensor(PeXX?,PeXY?,PeXZ?,PeYY?,PeYZ?,PeZZ?,''fac'',Bxyz?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('Peqq? = mms.rotate_tensor(PeXX?,PeXY?,PeXZ?,PeYY?,PeYZ?,PeZZ?,''fac'',Bxyz?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('Tefac? = mms.rotate_tensor(TeXX?,TeXY?,TeXZ?,TeYY?,TeYZ?,TeZZ?,''fac'',Bxyz?);',ic);

%% Compute tests for EDR
% Compute Q and Dng from Pepp
%c_eval('I1? = Pet?.data(:,1)+Pet?.data(:,4)+Pet?.data(:,6);',ic);
%c_eval('I2? = Pet?.data(:,1).*Pet?.data(:,4)+Pet?.data(:,1).*Pet?.data(:,6)+Pet?.data(:,4).*Pet?.data(:,6)-((Pet?.data(:,2)).^2+(Pet?.data(:,3)).^2+(Pet?.data(:,5)).^2);',ic);
%c_eval('Q? = 1-4*I2?./((I1?-PeZZp?.data).*(I1?+3*PeZZp?.data));',ic);
c_eval('Q? = (Pepp?.data(:,1,2).^2+Pepp?.data(:,1,3).^2+Pepp?.data(:,2,3).^2)./(Pepp?.data(:,2,2).^2+2*Pepp?.data(:,2,2).*Pepp?.data(:,1,1));',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(Pepp?.data(:,1,2).^2+Pepp?.data(:,1,3).^2+Pepp?.data(:,2,3).^2))./(Pepp?.data(:,1,1)+2*Pepp?.data(:,2,2));',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from Peqq
c_eval('agyro? = 2*abs(Peqq?.data(:,2,2)-Peqq?.data(:,3,3))./(Peqq?.data(:,2,2)+Peqq?.data(:,3,3));',ic);
c_eval('agyro? = irf.ts_scalar(ne?.time,agyro?);',ic);

% Compute temperature ratio An
c_eval('Temprat? = Pepp?.data(:,1,1)./(Pepp?.data(:,2,2));',ic);
c_eval('Temprat? = irf.ts_scalar(ne?.time,Temprat?);',ic);

% Compute electron Mach number
Units = irf_units; 
qe = Units.e;
me = Units.me;
c_eval('Ue? = Uevec?.abs.data;',ic);
c_eval('Veperp? = sqrt((Tefac?.data(:,2,2)+Tefac?.data(:,3,3))*qe/me);',ic);
c_eval('Me? = Ue?*1000./Veperp?;',ic);
c_eval('Me? = irf.ts_scalar(ne?.time,Me?);',ic);

% Compute current density and J.E
c_eval('Jmoms? = irf.ts_vec_xyz(Uevec?.time,1e18*qe*[ne?.data ne?.data ne?.data].*(Uivec?.data-Uevec?.data));',ic); % Current density in nA m^-2
c_eval('EdotJ? = dot(E?.data,Jmoms?.data,2)/1000;',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ? = irf.ts_scalar(ne?.time,EdotJ?);',ic);

% Calculate epsilon and delta parameters
c_eval('Bmag? = Bxyz?.abs.data;',ic);
c_eval('omegace? = qe*Bmag?/me*1e-9;',ic);
c_eval('EdotUe? = dot(E?.data,Uevec?.data,2);',ic);
c_eval('epsilone? = abs(6*pi*EdotUe?./(omegace?.*(Tefac?.data(:,1,1)+Tefac?.data(:,2,2)+Tefac?.data(:,3,3))));',ic);
c_eval('epsilone? = irf.ts_scalar(Uevec?.time,epsilone?);',ic);

c_eval('UexB? = cross(Uevec?,Bxyz?);',ic);
c_eval('UexB?.data = UexB?.data*1e-3;',ic);
c_eval('UexB?.data = E?.data+UexB?.data;',ic);
c_eval('UexBabs? = UexB?.abs.data;',ic);

c_eval('deltae? = 1e-3*UexBabs?./(Veperp?.*Bmag?*1e-9);',ic);
c_eval('deltae? = irf.ts_scalar(Uevec?.time,deltae?);',ic);

%% Plot figure

mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
nplots = 9;
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

h(3)=irf_panel('Dng'); set(h(2),'ColorOrder',mmsColors)
irf_pl_tx(h(3),'Dng?',1);
ylabel(h(3),'D_{ng}','Interpreter','tex');
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('agyro'); set(h(3),'ColorOrder',mmsColors)
irf_pl_tx(h(4),'agyro?',1);
ylabel(h(4),'A\phi_{e}','Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

h(5)=irf_panel('temprat'); set(h(4),'ColorOrder',mmsColors)
irf_pl_tx(h(5),'Temprat?',1);
ylabel(h(5),'An_{e}','Interpreter','tex');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k')

h(6)=irf_panel('Me'); set(h(5),'ColorOrder',mmsColors)
irf_pl_tx(h(6),'Me?',1);
ylabel(h(6),'M_{e \perp}','Interpreter','tex');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k')

h(7)=irf_panel('JdotE'); set(h(6),'ColorOrder',mmsColors)
irf_pl_tx(h(7),'EdotJ?',1);
ylabel(h(7),'E.J (nW m^{-3})','Interpreter','tex');
irf_legend(h(7),'(g)',[0.99 0.98],'color','k')

h(8)=irf_panel('epsilone'); set(h(7),'ColorOrder',mmsColors)
irf_pl_tx(h(8),'epsilone?',1);
ylabel(h(8),'\epsilon_{e}','Interpreter','tex');
irf_legend(h(8),'(h)',[0.99 0.98],'color','k')

h(9)=irf_panel('deltae'); set(h(8),'ColorOrder',mmsColors)
irf_pl_tx(h(9),'deltae?',1);
ylabel(h(9),'\delta_{e}','Interpreter','tex');
irf_legend(h(9),'(i)',[0.99 0.98],'color','k')

irf_plot_axis_align(h(1:nplots));
irf_zoom(h(1:nplots),'x',Tint);

%% save figure

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','EDRsignatures2.png');
