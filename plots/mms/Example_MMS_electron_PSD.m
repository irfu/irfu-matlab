% Script to plot electron PSD around pitch angles 0, 90, and 180 deg 
% and PSD versus pitch angle L1b brst data 
%
% Written by D. B. Graham

ic = 1; % Spacecraft number

Tintr = irf.tint('2015-10-30T05:15:20.00Z/2015-10-30T05:16:20.00Z');

%% Load data

tic;
c_eval('ePDist = mms.get_data(''PDe_fpi_brst_l2'',Tintr,?);',ic)
c_eval('Bxyz=mms.get_data(''B_dmpa_brst_l2'',Tintr,?);',ic);
c_eval('SCpot=mms.get_data(''V_edp_brst_l2'',Tintr,?);',ic);
toc;
ePDist = ePDist.convertto('s^3/km^6');
SCpot = SCpot.resample(ePDist);

%% Produce a single PAD at a selected time

tint = irf_time('2015-10-30T05:15:45.740000Z','utc>epochTT');
[paddist,thetapad,energypad,tintpad] = mms.get_pitchangledist(ePDist,Bxyz,tint,'angles',13); 
[~,idx] = min(abs(SCpot.time-tint));
energypad = energypad-SCpot.data(idx);

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
plot(h(1),energypad,paddist(:,1),'k',energypad,paddist(:,7),'r',energypad,paddist(:,13),'b');
ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)')
set(h(1),'yscale','log');
set(h(1),'xscale','log');
axis(h(1),[5 3e4 yrange])
set(h(1),'xtick',[1e0 1e1 1e2 1e3 1e4 1e5])
irf_legend(h(1),{'0 deg'},[0.91 0.92],'color','k')
irf_legend(h(1),{'90 deg'},[0.91 0.84],'color','r')
irf_legend(h(1),{'180 deg'},[0.91 0.76],'color','b')    

jetcolor = colormap('jet');
lll = length(jetcolor(:,1));
vcolors = floor(lll/length(energypad));
vcolors = (1:length(energypad))*vcolors;
    
c_eval('plot(h(2),thetapad,squeeze(paddist(?,:)),''color'',jetcolor(vcolors(?),:));',1);
hold(h(2),'on');
c_eval('plot(h(2),thetapad,squeeze(paddist(?,:)),''color'',jetcolor(vcolors(?),:));', 2:32);
hold(h(2),'off')
ylabel(h(2),'f_e (s^3 km^{-6})');
xlabel(h(2),'\theta (deg.)')
set(h(2),'yscale','log');
axis(h(2),[0 180 yrange])
set(h(2),'xtick',[0 45 90 135 180])
irf_zoom(h(2),'y',yrange);
tintutc = tintpad.utc;

title(h(1),strcat(tintutc(12:23),'UT'));
title(h(2),strcat(tintutc(12:23),'UT'));
set(gcf,'color','w');