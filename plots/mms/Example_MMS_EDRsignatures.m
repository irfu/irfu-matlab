% A routine to compute various parameters used to identify electron
% diffusion regions for the four MMS spacecraft.
% Written by D. B. Graham
%
% Quantities calculated so far are:
% (1) sqrt(Q) - Based on Swisdak, GRL ,2016. Values around 0.1 indicate
% electron agyrotropies. Computed based on the off-diagonal
% terms in the pressure tensor for Pe_perp1 = Pe_perp2.
% (1a) Dng - Based on Aunia et al., 2013; Computed based on the off-diagonal
% terms in the pressure tensor for Pe_perp1 = Pe_perp2. Similar to sqrt(Q)
% but with different normalization. Calculated but not plotted.
% (2) AG^(1/3) - Based on Che et al., POP, 2018. Constructed from determinant of
% field-aligned rotation of the electron pressure tensor (Pe_perp1 = Pe_perp2).
% (3) A phi_e/2 = abs(Perp1-Perp2)/(Perp1+Perp2).
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
% thermal speed perpendicular to B. Values of O(1) are expected in EDRs (Scudder et al., 2012, 2015).
% (6) J.E' - J.E > 0 is expected in the electron diffusion region,
% corresponding to dissipation of field energy. J is calculated on each
% spacecraft using the particle moments (Zenitani et al., PRL, 2011).
% (7) epsilon_e - Energy gain per cyclotron period. Values of O(1) are
% expected in EDRs (Scudder et al., 2012, 2015).
% (8) delta_e - Relative strength of the electric and magnetic force in the
% bulk electron rest frame. N. B. Very sensitive to electron moments and
% electric field. Check version of these quantities (Scudder et al., 2012, 2015).
%
% Notes: kappa_e (not yet included) is taken to be the largest value of
% epsilon_e and delta_e at any given point.
% Requires electron distributions with version number v2.0.0 or higher.
% Calculations of agyrotropy measures (1)--(3) become unreliable at low
% densities n_e <~ 2 cm^-3, when the raw particle counts are low.
% Agyrotropies are removed for n_e < 1 cm^-3

%% Time interval selection
Tint = irf.tint('2017-07-06T00:54:00.00Z/2017-07-06T00:54:45.00Z');
%Tint = irf.tint('2015-12-14T01:17:38.00Z/2015-12-14T01:17:41.00Z');

%% Load data
ic = 1:4;
%c_eval('Bxyz?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',Tint);',ic);
%c_eval('E? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',Tint);',ic);
c_eval('Bxyz?=mms.get_data(''B_dmpa_fgm_srvy_l2'',Tint,?);',ic);
c_eval('E? =mms.get_data(''E_dsl_edp_brst_l2'',Tint, ?);',ic);

for ii=1:4
  c_eval('ne?=mms.get_data(''Ne_fpi_brst_l2'',Tint,?);',ii);
  c_eval('Uevec? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',Tint,?);',ii);
  c_eval('Te? = mms.get_data(''Te_dbcs_fpi_brst_l2'',Tint,?);',ii);
  c_eval('Pe? = mms.get_data(''Pe_dbcs_fpi_brst_l2'',Tint,?);',ii);
  c_eval('Uivec? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',Tint,?);',ii);
end

c_eval('E? = E?.resample(ne?);',ic);
c_eval('Bxyz? = Bxyz?.resample(ne?);',ic);
c_eval('Uivec? = Uivec?.resample(ne?);',ic);

%% Rotate pressure and temperature tensors
c_eval('Pepp? = mms.rotate_tensor(Pe?,''fac'',Bxyz?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('Peqq? = mms.rotate_tensor(Pe?,''fac'',Bxyz?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('Tefac? = mms.rotate_tensor(Te?,''fac'',Bxyz?);',ic);

%% Compute tests for EDR
% Compute Q and Dng from Pepp
c_eval('Q? = (Pepp?.xy.data.^2+Pepp?.xz.data.^2+Pepp?.yz.data.^2)./(Pepp?.yy.data.^2+2*Pepp?.yy.data.*Pepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(Pepp?.xy.data.^2+Pepp?.xz.data.^2+Pepp?.yz.data.^2))./(Pepp?.xx.data+2*Pepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy measure AG1/3
c_eval('detP? = Pepp?.xx.data.*Pepp?.yy.data.^2 - Pepp?.xx.data.*Pepp?.yz.data.^2 - Pepp?.yy.data.*(Pepp?.xy.data.^2 + Pepp?.xz.data.^2) + 2*Pepp?.xy.data.*Pepp?.xz.data.*Pepp?.yz.data;',ic);
c_eval('detG? = Pepp?.xx.data.*Pepp?.yy.data.^2;',ic);
c_eval('AG? = abs(detP? - detG?)./(detP? + detG?);',ic);
c_eval('AGcr? = irf.ts_scalar(ne?.time,AG?.^(1/3));',ic);
c_eval('AG? = irf.ts_scalar(ne?.time,AG?);',ic);

% Compute agyrotropy Aphi from Peqq
c_eval('agyro? = abs(Peqq?.yy.data-Peqq?.zz.data)./(Peqq?.yy.data+Peqq?.zz.data);',ic);
c_eval('agyro? = irf.ts_scalar(ne?.time,agyro?);',ic);

% Simple fix to remove spurious points
for xx=1:4
  c_eval('Qdata? = Q?.data;',xx);
  c_eval('for ii=2:1:length(Q?.data)-1;if (Q?.data(ii) > 2*Q?.data(ii-1)) && (Q?.data(ii) > 2*Q?.data(ii+1)); Qdata?(ii) = NaN; end; end;',xx);
  c_eval('Q?.data = Qdata?;',xx);
  c_eval('Dngdata? = Dng?.data;',xx);
  c_eval('for ii=2:1:length(Dng?.data)-1;if (Dng?.data(ii) > 2*Dng?.data(ii-1)) && (Dng?.data(ii) > 2*Dng?.data(ii+1)); Dngdata?(ii) = NaN; end; end;',xx);
  c_eval('Dng?.data = Dngdata?;',xx);
  c_eval('agyrodata? = agyro?.data;',xx);
  c_eval('for ii=2:1:length(agyro?.data)-1;if (agyro?.data(ii) > 2*agyro?.data(ii-1)) && (agyro?.data(ii) > 2*agyro?.data(ii+1)); agyrodata?(ii) = NaN; end; end;',xx);
  c_eval('agyro?.data = agyrodata?;',xx);
  c_eval('AGcrdata? = AGcr?.data;',xx);
  c_eval('for ii=2:1:length(AGcr?.data)-1;if (AGcr?.data(ii) > 2*AGcr?.data(ii-1)) && (AGcr?.data(ii) > 2*AGcr?.data(ii+1)); AGcrdata?(ii) = NaN; end; end;',xx);
  c_eval('AGcr?.data = AGcrdata?;',xx);
end

% remove all points corresponding to densities below 1cm^-3
c_eval('rmpnts? = ones(size(ne?.data));',ic);
c_eval('rmpnts?(ne?.data < 1) = NaN;',ic);
c_eval('Q?.data = Q?.data.*rmpnts?;',ic);
c_eval('agyro?.data = agyro?.data.*rmpnts?;',ic);
c_eval('Dng?.data = Dng?.data.*rmpnts?;',ic);
c_eval('AGcr?.data = AGcr?.data.*rmpnts?;',ic);

% Compute temperature ratio An
c_eval('Temprat? = Pepp?.xx/(Pepp?.yy);',ic);

% Compute electron Mach number
Units = irf_units;
qe = Units.e;
me = Units.me;
c_eval('Ue? = Uevec?.abs.data;',ic);
c_eval('Veperp? = sqrt((Tefac?.yy.data+Tefac?.zz.data)*qe/me);',ic);
c_eval('Me? = Ue?*1000./Veperp?;',ic);
c_eval('Me? = irf.ts_scalar(ne?.time,Me?);',ic);

% Compute current density and J.E
c_eval('Jmoms? = irf.ts_vec_xyz(Uevec?.time,1e18*qe*[ne?.data ne?.data ne?.data].*(Uivec?.data-Uevec?.data));',ic); % Current density in nA m^-2
c_eval('UexB? = cross(Uevec?,Bxyz?);',ic);
c_eval('UexB?.data = UexB?.data*1e-3;',ic);
c_eval('UexB?.data = E?.data+UexB?.data;',ic);
c_eval('EdotJ? = dot(UexB?.data,Jmoms?.data,2)/1000;',ic); %J (nA/m^2), E (mV/m), E.J (nW/m^3)
c_eval('EdotJ? = irf.ts_scalar(ne?.time,EdotJ?);',ic);

% Calculate epsilon and delta parameters
c_eval('Bmag? = Bxyz?.abs.data;',ic);
c_eval('omegace? = qe*Bmag?/me*1e-9;',ic);
c_eval('EdotUe? = dot(E?.data,Uevec?.data,2);',ic);
c_eval('epsilone? = abs(6*pi*EdotUe?./(omegace?.*(Tefac?.trace.data)));',ic);
c_eval('epsilone? = irf.ts_scalar(Uevec?.time,epsilone?);',ic);

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

h(3)=irf_panel('AGcr'); set(h(2),'ColorOrder',mmsColors)
irf_pl_tx(h(3),'AGcr?',1);
ylabel(h(3),'AG^{1/3}','Interpreter','tex');
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

h(4)=irf_panel('agyro'); set(h(3),'ColorOrder',mmsColors)
irf_pl_tx(h(4),'agyro?',1);
ylabel(h(4),'A\phi_{e}/2','Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

h(5)=irf_panel('temprat'); set(h(4),'ColorOrder',mmsColors)
irf_pl_tx(h(5),'Temprat?',1);
ylabel(h(5),'T_{e||}/T_{e \perp}','Interpreter','tex');
irf_legend(h(5),'(e)',[0.99 0.98],'color','k')

h(6)=irf_panel('Me'); set(h(5),'ColorOrder',mmsColors)
irf_pl_tx(h(6),'Me?',1);
ylabel(h(6),'M_{e \perp}','Interpreter','tex');
irf_legend(h(6),'(f)',[0.99 0.98],'color','k')

h(7)=irf_panel('JdotE'); set(h(6),'ColorOrder',mmsColors)
irf_pl_tx(h(7),'EdotJ?',1);
ylabel(h(7),'E''.J (nW m^{-3})','Interpreter','tex');
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
%print('-dpng','-painters','-r300','EDRsignatures.png');
