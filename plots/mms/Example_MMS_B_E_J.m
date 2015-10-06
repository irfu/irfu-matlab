% Plots of B, J, E, JxB electric field, and J.E. Calculates J using
% Curlometer method. 
% Written by D. B. Graham

mms.db_init('local_file_db','/data/mms');

tint = irf.tint('2015-10-01T06:30:00.00Z/2015-10-01T07:20:00.00Z');

ic = 1; %get density from this spacecraft

c_eval('Bxyz?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',[1:4]);
c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',[2:4]);
Bxyzav = (Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4;
Bxyzav = TSeries(Bxyz1.time,Bxyzav,'to',1);

c_eval('Exyz?=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',[1:4]);
c_eval('Exyz? = Exyz?.resample(Exyz1);',[2:4]);
c_eval('Exyz?.data(find(abs(Exyz?.data) > 100)) = NaN;',[1:4]); %Remove some questionable fields
Exyzav = (Exyz1.data+Exyz2.data+Exyz3.data+Exyz4.data)/4;
Exyzav = TSeries(Exyz1.time,Exyzav,'to',1);

c_eval('ne=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);
ne = TSeries(ne.time,ne.data,'to',1);
ne = ne.resample(Bxyz1);

%old data files 
%load('/data/mms/irfu/mmsR.mat');
%time = EpochTT(R.time);
%c_eval('Rxyz? = TSeries(time,R.gseR?,''to'',1);',[1:4]);
%c_eval('Rxyz? = Rxyz?.tlim(tint);',[1:4]);
%c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',[1:4]);
%clear R;

R  = mms.get_data('R_gse',tint);
c_eval('Rxyz? = TSeries(R.time,R.gseR?,''to'',1);',[1:4]);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',[1:4]);

% Assuming GSE and DMPA are the same coordinate system.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');

divovercurl = divB;
divovercurl.data = abs(divovercurl.data)./j.abs.data;

% Transform current density into field-aligned coordinates
SCpos = [0 1 0];

Bmag = Bxyz1.abs.data;
Rpar = Bxyz1.data./[Bmag Bmag Bmag];
Rperpy = irf_cross(Rpar,SCpos);
Rmag   = irf_abs(Rperpy,1);
Rperpy = Rperpy./[Rmag Rmag Rmag];
Rperpx = irf_cross(Rperpy, Rpar);
Rmag   = irf_abs(Rperpx,1);
Rperpx = Rperpx./[Rmag Rmag Rmag];

jpar = dot(Rpar,j.data,2);
jperp = dot(Rperpx,j.data,2);
jperp2 = dot(Rperpy,j.data,2);

jfac = TSeries(j.time,[jperp jperp2 jpar],'to',1);


h = irf_plot(7,'newfigure');

hca = irf_panel('BMMS1');
irf_plot(hca,Bxyzav);
ylabel(hca,{'B_{DMPA}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

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
set(hca,'yscale','log');

hca = irf_panel('EMMS1');
irf_plot(hca,Exyzav);
ylabel(hca,{'E_{DSL}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(hca,{'E_{x}','E_{y}','E_{z}'},[0.88 0.10])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('jxB');
jxB.data = jxB.data./[ne.data ne.data ne.data]; 
jxB.data = jxB.data/1.6e-19/1000; %Convert to (mV/m)
jxB.data(find(abs(jxB.data) > 100)) = NaN; % Remove some questionable fields
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

title(h(1),strcat('MMS',num2str(ic),'- Current density and fields'));

irf_plot_axis_align(1,h(1:7))
irf_zoom(h(1:7),'x',tint);