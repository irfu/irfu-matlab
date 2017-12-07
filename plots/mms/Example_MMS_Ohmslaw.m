% Compute the terms in the generalized Ohm's law equation: Ion convection,
% Hall, and electron pressure divergence terms. Hall and pressure terms are
% computed using four-spacecraft methods. The observed electric fields and
% convection terms are averaged over the four spacecraft. Terms computed in GSE
% coordinates.

Tint = irf.tint('2015-10-30T05:15:40.00Z/2015-10-30T05:15:55.00Z');
Tintlong = Tint+[-60 60]; % Take longer time interval for loading position data
%% Load all data and constants
ic = 1:4; % Use all MMS spacecraft

% Define constants
Units = irf_units; 
qe = Units.e;
me = Units.me;

% Load FPI data
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',Tint,?);',ic);
c_eval('Ve? = mms.get_data(''Ve_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('Pe? = mms.get_data(''Pe_gse_fpi_brst_l2'',Tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',Tint,?);',ic);
c_eval('Vi? = mms.get_data(''Vi_gse_fpi_brst_l2'',Tint,?);',ic);

% Load FGM data
c_eval('Bxyz? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',Tint);',ic);

% Load spacecraft position
R  = mms.get_data('R_gse',Tintlong);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

% Load Electric field
c_eval('Exyz? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',Tint);',ic);

%% Resample and compute averages
c_eval('ne? = ne?.resample(ne1);',2:4);
c_eval('Ve? = Ve?.resample(ne1);',ic);
c_eval('Pe? = Pe?.resample(ne1);',ic);
c_eval('ni? = ni?.resample(ne1);',ic);
c_eval('Vi? = Vi?.resample(ne1);',ic);
c_eval('Exyz? = Exyz?.resample(ne1);',ic);
c_eval('Bxyz? = Bxyz?.resample(ne1);',ic);
c_eval('Rxyz? = Rxyz?.resample(ne1);',ic);

neav = irf.ts_scalar(ne1.time,(ne1.data+ne2.data+ne3.data+ne4.data)/4);
Bxyzav = irf.ts_vec_xyz(Bxyz1.time,(Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4);
Exyzav = irf.ts_vec_xyz(Exyz1.time,(Exyz1.data+Exyz2.data+Exyz3.data+Exyz4.data)/4);

%% Compute convection terms
% Compute ion convection term (MMS average)
c_eval('EvxBi? = cross(Vi?,Bxyz?);',ic);
c_eval('EvxBi?.data = -EvxBi?.data*1e-3;',ic);
EvxBi = irf.ts_vec_xyz(EvxBi1.time,(EvxBi1.data+EvxBi2.data+EvxBi3.data+EvxBi4.data)/4);

% Compute electron convection term (MMS average)
c_eval('EvxBe? = cross(Ve?,Bxyz?);',ic);
c_eval('EvxBe?.data = -EvxBe?.data*1e-3;',ic);
EvxBe = irf.ts_vec_xyz(EvxBe1.time,(EvxBe1.data+EvxBe2.data+EvxBe3.data+EvxBe4.data)/4);

%% Compute pressure divergence term
c_eval('Pe?.data = Pe?.data/1e9;',ic); % Unit conversion
EPeXX = c_4_grad('Rxyz?','Pe?.xx','grad');
EPeXY = c_4_grad('Rxyz?','Pe?.xy','grad');
EPeXZ = c_4_grad('Rxyz?','Pe?.xz','grad');
EPeYY = c_4_grad('Rxyz?','Pe?.yy','grad');
EPeYZ = c_4_grad('Rxyz?','Pe?.yz','grad');
EPeZZ = c_4_grad('Rxyz?','Pe?.zz','grad');
EPeX = -(EPeXX.data(:,1)+EPeXY.data(:,2)+EPeXZ.data(:,3))./(neav.data*1e6*qe);
EPeY = -(EPeXY.data(:,1)+EPeYY.data(:,2)+EPeYZ.data(:,3))./(neav.data*1e6*qe);
EPeZ = -(EPeXZ.data(:,1)+EPeYZ.data(:,2)+EPeZZ.data(:,3))./(neav.data*1e6*qe);
EPe = irf.ts_vec_xyz(EPeXX.time,[EPeX EPeY EPeZ]);

%% Compute Hall term and current density using curlometer
[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');
jxB = cross(j,Bxyzav);
jxB.data = jxB.data*1e-9;
jxB.data = jxB.data./[neav.data neav.data neav.data]; 
jxB.data = jxB.data/qe/1000; %Convert to (mV/m)
j.data = j.data*1e9;

%% Plot figure
h = irf_plot(5,'newfigure');
xSize=750; ySize=650;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.85;
ywidth = 0.18;
set(h(1),'position',[0.11 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.11 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.11 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.11 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.11 0.97-5*ywidth xwidth ywidth]);

h(1) = irf_panel('BMMSall');
irf_plot(h(1),Bxyzav);
ylabel(h(1),{'B','(nT)'},'Interpreter','tex');
irf_legend(h(1),{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
irf_legend(h(1),'(a)',[0.99 0.90],'color','k')

h(2) = irf_panel('Jxyz');
irf_plot(h(2),j);
ylabel(h(2),{'J','(nA m^{-2})'},'Interpreter','tex');
irf_legend(h(2),{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3) = irf_panel('EMMS1');
irf_plot(h(3),Exyzav);
ylabel(h(3),{'E','(mV m^{-1})'},'Interpreter','tex');
irf_legend(h(3),{'E_{x}','E_{y}','E_{z}'},[0.88 0.95])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')

Exterms = TSeries(Exyzav.time,[Exyzav.data(:,1) jxB.data(:,1) EvxBi.data(:,1) EPe.data(:,1) EvxBe.data(:,1)]);

h(4) = irf_panel('Exterms');
irf_plot(h(4),Exterms);
ylabel(h(4),{'E_{x}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(h(4),{'E','J \times B/q_{e}n','-V_{i} \times B','-\nabla \cdot P_{e}/q_{e}n','-V_{e} \times B'},[0.9 0.98])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

ELHS = irf.ts_vec_xyz(Exyzav.time,Exyzav.data-EvxBi.data);
ERHS = irf.ts_vec_xyz(Exyzav.time,jxB.data+EPe.data);

h(5) = irf_panel('ELHSvsRHS');
irf_plot(h(5),irf.ts_scalar(ELHS.time,ELHS.data(:,1)),'color','k');
hold(h(5),'on');
irf_plot(h(5),irf.ts_scalar(ELHS.time,ERHS.data(:,1)),'color','r');
hold(h(5),'off');
ylabel(h(5),{'E_{x}','(mV m^{-1})'},'Interpreter','tex');
irf_legend(h(5),'E+V_{i} \times B',[0.01 0.98])
irf_legend(h(5),'J \times B/q_{e}n - \nabla \cdot P_{e}/q_{e}n',[0.13 0.98],'color','r')
irf_legend(h(5),'(e)',[0.99 0.98],'color','k')

title(h(1),'MMS - 4 Spacecraft average')

irf_plot_axis_align(1,h(1:5))
irf_zoom(h(1:5),'x',Tint);

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r400','Ohmslaw.png');