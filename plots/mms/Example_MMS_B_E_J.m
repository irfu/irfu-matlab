% Plots of B, J, E, JxB electric field, and J.E. Calculates J using
% Curlometer method. 
% Written by D. B. Graham

Tint = irf.tint('2016-06-07T01:53:44.00Z/2016-06-07T19:30:34.00Z');

%%
ic = 1:4; 

c_eval('Bxyz?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',Tint);');
c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',2:4);
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);

c_eval('Exyz?=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',Tint);');
c_eval('Exyz? = Exyz?.resample(Exyz1);',2:4);
%c_eval('[Exyz?,~]=irf_edb(Exyz?,Bxyz?,15,''E.B=0'');',[1:4]); % Removes some wake fields
Exyzav = irf.ts_vec_xyz(Exyz1.time,(Exyz1.data+Exyz2.data+Exyz3.data+Exyz4.data)/4);

%c_eval('ni?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',Tint);'); % For MP phases
c_eval('ni?=mms.db_get_ts(''mms?_hpca_srvy_sitl_moments'',''mms?_hpca_hplus_number_density'',Tint);'); % For 1x phase
c_eval('ni? = ni?.resample(ni1);',2:4);
ni = irf.ts_scalar(ni1.time,(ni1.data+ni2.data+ni3.data+ni4.data)/4);
ni = ni.resample(Bxyz1);

R  = mms.get_data('R_gse',Tint);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?);');
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);');

% Assuming GSE and DMPA are the same coordinate system.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');

divovercurl = divB;
divovercurl.data = abs(divovercurl.data)./j.abs.data;

% Transform current density into field-aligned coordinates
jfac = irf_convert_fac(j,Bxyzav,[1 0 0]);
jfac.coordinateSystem = 'FAC';

%%
h = irf_plot(8,'newfigure');

hca = irf_panel('BMMS');
irf_plot(hca,Bxyzav);
ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
irf_zoom(hca,'y',[-70 70]);
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
hca = irf_panel('niMMS4'); set(hca,'ColorOrder',mmsColors)
irf_pl_tx(hca,'ni?',1);
ylabel(hca,{'n_i','(cm^{-3})'},'Interpreter','tex');
set(hca,'yscale','log');
irf_zoom(hca,'y',[1e-4 10]);
irf_legend(hca,'(b)',[0.99 0.98],'color','k')
irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.99 0.1],'color','cluster')

hca = irf_panel('J');
j.data = j.data*1e9;
irf_plot(hca,j);
ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca = irf_panel('Jfac');
jfac.data = jfac.data*1e9;
irf_plot(hca,jfac);
ylabel(hca,{'J_{FAC}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{\perp 1}','J_{\perp 2}','J_{||}'},[0.88 0.10])
irf_legend(hca,'(d)',[0.99 0.98],'color','k')

hca = irf_panel('divovercurl');
irf_plot(hca,divovercurl);
ylabel(hca,{'|\nabla . B|','|\nabla \times B|'},'Interpreter','tex');
irf_legend(hca,'(e)',[0.99 0.98],'color','k')

hca = irf_panel('EMMS1');
irf_plot(hca,Exyzav);
ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('jxB');
jxB.data = jxB.data./[ni.data ni.data ni.data]; 
jxB.data = jxB.data/1.6e-19/1000; %Convert to (mV/m)
jxB.data(abs(jxB.data) > 100) = NaN; % Remove some questionable fields
irf_plot(hca,jxB);
ylabel(hca,{'J \times B/n_{e} q_{e}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,'(f)',[0.99 0.98],'color','k')

j = j.resample(Exyzav);
EdotJ = dot(Exyzav.data,j.data,2)/1000; %J (nA/m^2), E (mV/m), E.J (nW/m^3)
EdotJ = TSeries(Exyzav.time,EdotJ);

hca = irf_panel('jdotE');
irf_plot(hca,EdotJ);
ylabel(hca,{'E . J','(nW m^{-3})'},'Interpreter','tex');
irf_legend(hca,'(g)',[0.99 0.98],'color','k')

title(h(1),'MMS - Current density and fields');

irf_plot_axis_align(1,h(1:8))
irf_zoom(h(1:8),'x',Tint);