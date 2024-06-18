% Example of dispersion relation properties using the four MMS spacecraft.
% Default example is of whistler waves.
%
% Notes: Should only be used when the half wavelength of the wave is larger
% than the spacecraft separations, otherwise spatial aliasing will occur.

%% Load data
ic = 1:4;
Tint = irf.tint('2015-10-16T13:05:24.000Z/2015-10-16T13:05:50.000Z');
% Take longer time than used for the dispersion relation, so edge effects are not included
c_eval('Bxyz?=mms.get_data(''B_gse_brst_l2'',Tint,?);',ic);
c_eval('Bscm?=mms.get_data(''B_gse_scm_brst_l2'',Tint,?);',ic);
Tintlong = Tint+[-60 60];
R  = mms.get_data('R_gse',Tintlong);
c_eval('Rxyz? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

c_eval('Bscmfac? = irf_convert_fac(Bscm?, Bxyz?, [1, 0, 0]);',ic);

%% Compute dispersion relation
Tints = irf.tint('2015-10-16T13:05:26.500Z/2015-10-16T13:05:27.000Z');

[xvecs,yvecs,Power] = mms.fk_powerspec4SC('Bscmfac?.x','Rxyz?','Bxyz?',Tints,'linear',10,'numk',500,'cav',4,'wwidth',2);

%% Plot figure all quantities

fn=figure;
xwidth = 0.20;
ywidth = 0.38;

set(fn,'Position',[10 10 1400 600])
h(1)=axes('position',[0.035 0.58 xwidth ywidth]); % [x y dx dy]
h(2)=axes('position',[0.285 0.58 xwidth ywidth]);
h(3)=axes('position',[0.535 0.58 xwidth ywidth]);
h(4)=axes('position',[0.785 0.58 xwidth ywidth]);
h(5)=axes('position',[0.035 0.09 xwidth ywidth]);
h(6)=axes('position',[0.285 0.09 xwidth ywidth]);
h(7)=axes('position',[0.535 0.09 xwidth ywidth]);
h(8)=axes('position',[0.785 0.09 xwidth ywidth]);
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2);
set(fn,'defaultAxesFontSize',14)

pcolor(h(1),xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf));
shading(h(1),'flat');
xlabel(h(1),'|k| (m^{-1})');
ylabel(h(1),'f (kHz)');
c=colorbar('peer',h(1),'ver');
ylabel(c,'log_{10} P(f,k)/P_{max}');
colormap('jet');
irf_legend(h(1),'(a)',[0.99 0.99],'color','k','fontsize',14)
vph = mms.estimate_phase_speed(Power.Powerkmagf,yvecs.fkmag,xvecs.kmag,100);
kfit = 2*pi*yvecs.fkmag/vph;
hold(h(1),'on')
plot(h(1),kfit,yvecs.fkmag*1e-3,'k','linewidth',1)
hold(h(1),'off')
axis(h(1),[0 3e-4 0 2])
caxis(h(1),[-6 0]);
irf_legend(h(1),strcat('v_{ph} = ',num2str(round(vph/1e3)),'km s^{-1}'),[0.9 0.05],'color','k','fontsize',14)
set(gcf,'color','w')

pcolor(h(2),xvecs.kperp,yvecs.kpar,log10(Power.Powerkperpkpar));
shading(h(2),'flat');
xlabel(h(2),'|k_{\perp}| (m^{-1})');
ylabel(h(2),'k_{||} (m^{-1})');
c=colorbar('peer',h(2),'ver');
ylabel(c,'log_{10} P(k_{\perp},k_{||})/P_{max}');
colormap('jet');
caxis(h(2),[-6 0]);
irf_legend(h(2),'(b)',[0.99 0.99],'color','k','fontsize',14)

pcolor(h(3),xvecs.kxkxky,yvecs.kykxky,log10(Power.Powerkxky));
shading(h(3),'flat');
hold(h(3),'on')
plot(h(3),[0 0],max(yvecs.kykxky)*[-1 1],'color','k','Linewidth',1)
plot(h(3),max(xvecs.kxkxky)*[-1 1],[0 0],'color','k','Linewidth',1)
hold(h(3),'off')
xlabel(h(3),'k_{x} (m^{-1})');
ylabel(h(3),'k_{y} (m^{-1})');
c=colorbar('peer',h(3),'ver');
ylabel(c,'log_{10} P(k_x,k_y)/P_{max}');
colormap('jet');
caxis(h(3),[-6 0]);
irf_legend(h(3),'(c)',[0.99 0.99],'color','k','fontsize',14)

pcolor(h(4),xvecs.kxkxkz,yvecs.kzkxkz,log10(Power.Powerkxkz));
shading(h(4),'flat');
hold(h(4),'on')
plot(h(4),[0 0],max(yvecs.kzkxkz)*[-1 1],'color','k','Linewidth',1)
plot(h(4),max(xvecs.kxkxkz)*[-1 1],[0 0],'color','k','Linewidth',1)
hold(h(4),'off')
xlabel(h(4),'k_{x} (m^{-1})');
ylabel(h(4),'k_{z} (m^{-1})');
c=colorbar('peer',h(4),'ver');
ylabel(c,'log_{10} P(k_x,k_z)/P_{max}');
colormap('jet');
caxis(h(4),[-6 0]);
irf_legend(h(4),'(d)',[0.99 0.99],'color','k','fontsize',14)

pcolor(h(5),xvecs.kykykz,yvecs.kzkykz,log10(Power.Powerkykz));
shading(h(5),'flat');
hold(h(5),'on')
plot(h(5),[0 0],max(yvecs.kzkykz)*[-1 1],'color','k','Linewidth',1)
plot(h(5),max(xvecs.kykykz)*[-1 1],[0 0],'color','k','Linewidth',1)
hold(h(5),'off')
xlabel(h(5),'k_{y} (m^{-1})');
ylabel(h(5),'k_{z} (m^{-1})');
c=colorbar('peer',h(5),'ver');
ylabel(c,'log_{10} P(k_y,k_z)/P_{max}');
colormap('jet');
caxis(h(5),[-6 0]);
irf_legend(h(5),'(e)',[0.99 0.99],'color','k','fontsize',14)

pcolor(h(6),xvecs.kxf,yvecs.fkxf*1e-3,log10(Power.Powerkxf));
shading(h(6),'flat');
xlabel(h(6),'k_{x} (m^{-1})');
ylabel(h(6),'f (kHz)');
c=colorbar('peer',h(6),'ver');
ylabel(c,'log_{10} P(k_x,f)/P_{max}');
colormap('jet');
caxis(h(6),[-6 0]);
irf_legend(h(6),'(f)',[0.99 0.99],'color','k','fontsize',14)

pcolor(h(7),xvecs.kyf,yvecs.fkyf*1e-3,log10(Power.Powerkyf));
shading(h(7),'flat');
xlabel(h(7),'k_{y} (m^{-1})');
ylabel(h(7),'f (kHz)');
c=colorbar('peer',h(7),'ver');
ylabel(c,'log_{10} P(k_y,f)/P_{max}');
colormap('jet');
caxis(h(7),[-6 0]);
irf_legend(h(7),'(g)',[0.99 0.99],'color','k','fontsize',14)

pcolor(h(8),xvecs.kzf,yvecs.fkzf*1e-3,log10(Power.Powerkzf));
shading(h(8),'flat');
xlabel(h(8),'k_{z} (m^{-1})');
ylabel(h(8),'f (kHz)');
c=colorbar('peer',h(8),'ver');
ylabel(c,'log_{10} P(k_z,f)/P_{max}');
colormap('jet');
caxis(h(8),[-6 0]);
irf_legend(h(8),'(h)',[0.99 0.99],'color','k','fontsize',14)
