% Script to plot electron PSD around pitch angles 0, 90, and 180 deg 
% and PSD versus pitch angle L1b brst data 
%
% Written by D. B. Graham

ic = 1; % Spacecraft number

Tintr = irf.tint('2015-12-30T00:30:00.00Z/2015-12-30T00:30:33.00Z');

%% Load data

tic;
c_eval('diste = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',Tintr);',ic);
c_eval('energy0=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_energy0'',Tintr);',ic);
c_eval('energy1=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_energy1'',Tintr);',ic);
c_eval('phi=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_phi'',Tintr);',ic);
c_eval('theta=mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_theta'',Tintr);',ic);
c_eval('stepTable=mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_stepTable_parity'',Tintr);',ic);
toc;

c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);

%% Produce a single PAD at a selected time

tint = irf_time('2015-12-30T00:30:11.010000Z','utc>epochTT');
[paddist,thetapad,energypad,tintpad] = mms.get_pitchangledist(diste,phi,theta,stepTable,energy0,energy1,Bxyz,tint); 
paddist = paddist*1e30; %convert to commonly used s^3 km^-6

%% Plot PAD

fn=figure;
set(fn,'Position',[10 10 600 250])
    h(1)=axes('position',[0.08 0.12 0.4 0.82]); 
    h(2)=axes('position',[0.58 0.12 0.4 0.82]); 
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 

ymin = 10^-4;
ymax = ceil(max(max(log10(paddist))));
yrange = [ymin 10^ymax];
plot(h(1),energypad,paddist(:,1),'k',energypad,mean(paddist(:,[6 7]),2),'r',energypad,paddist(:,12),'b');
ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)')
set(h(1),'yscale','log');
set(h(1),'xscale','log');
irf_zoom(h(1),'y',yrange);
irf_zoom(h(1),'x',[10 3e4]);
irf_legend(h(1),{'0 deg'},[0.91 0.92],'color','k')
irf_legend(h(1),{'90 deg'},[0.91 0.84],'color','r')
irf_legend(h(1),{'180 deg'},[0.91 0.76],'color','b')    

jetcolor = colormap('jet');
lll = length(jetcolor(:,1));
vcolors = floor(lll/length(energypad));
vcolors = [1:length(energypad)]*vcolors;
    
c_eval('plot(h(2),thetapad,squeeze(paddist(?,:)),''color'',jetcolor(vcolors(?),:));',1);
hold(h(2),'on');
c_eval('plot(h(2),thetapad,squeeze(paddist(?,:)),''color'',jetcolor(vcolors(?),:));',[2:32]);
hold(h(2),'off')
ylabel(h(2),'f_e (s^3 km^{-6})');
xlabel(h(2),'\theta (deg.)')
set(h(2),'yscale','log');
irf_zoom(h(2),'x',[0 180]);
irf_zoom(h(2),'y',yrange);
tintutc = tintpad.utc;

title(h(1),strcat(tintutc(12:23),'UT'));
title(h(2),strcat(tintutc(12:23),'UT'));
set(gcf,'color','w');