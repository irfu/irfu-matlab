% Script to plot electron PSD around pitch angles 0, 90, and 180 deg using
% L1b brst data. 
% Written by D. B. Graham

ic = 1; % Spacecraft number

mms.db_init('local_file_db','/data/mms');

tint = irf.tint('2015-09-19T10:05:32.900000Z/2015-09-19T10:05:33.000000Z'); % Time interval to average over
tintl = tint+[-1 1];

%%

c_eval('disttemp = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint);',ic);
% This is rather slow so the data can be read locally using these commands:
% tmpDataObj = dataobj('data/mms1_fpi_brst_l1b_des-dist_20150919100500_v0.2.0.cdf');
% dist = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_brstSkyMap_dist'));
% disttemp = dist.tlim(tint);

dist1 = squeeze(irf.nanmean(disttemp.data,1));


dangle = 180/16;
phi = dangle*[0:31]+dangle/2;
theta = dangle*[0:15]+dangle/2;
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

x = -cosd(phi')*sind(theta);
y = -sind(phi')*sind(theta);
z = -ones(length(phi),1)*cosd(theta);

tid = tint(1);
c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tintl);',ic);
Bmag = Bxyz.abs.data;
[~,tB] = min(abs(Bxyz.time-tid));
Bvec = Bxyz.data(tB,:)/Bmag(tB)

%Calculate angles between B and theta and phi
thetab = acosd(x*Bvec(1)+y*Bvec(2)+z*Bvec(3));
pospar = ones(length(phi),length(theta)); 
pospar(find(thetab > 15)) = NaN;
posperp = ones(length(phi),length(theta)); 
posperp(find(thetab < 82.5)) = NaN;
posperp(find(thetab > 97.5)) = NaN;
posapar = ones(length(phi),length(theta)); 
posapar(find(thetab < 165)) = NaN;

distpar = dist1;
distperp = dist1;
distapar = dist1;

for ii = 1:length(energy);
    distpar(ii,:,:)  = squeeze(dist1(ii,:,:)).*pospar;
    distperp(ii,:,:) = squeeze(dist1(ii,:,:)).*posperp;
    distapar(ii,:,:) = squeeze(dist1(ii,:,:)).*posapar;
end

distpar =  squeeze(irf.nanmean(irf.nanmean(distpar,3),2))*1e30; %Convert units (s^3 cm^-6 -> s^3 km^-6)
distperp = squeeze(irf.nanmean(irf.nanmean(distperp,3),2))*1e30;
distapar = squeeze(irf.nanmean(irf.nanmean(distapar,3),2))*1e30;

fn=figure;
set(fn,'Position',[10 10 400 300])
    h(1)=axes('position',[0.15 0.2 0.80 0.75]); 
    ud=get(fn,'userdata');
    ud.subplot_handles=h;
    set(fn,'userdata',ud);
    set(fn,'defaultLineLineWidth',2); 

h(1)=irf_panel('EpsdMMS');
plot(energy,distpar,'k',energy,distperp,'r',energy,distapar,'b');
ylabel(h(1),'f_e (s^3 km^{-6})');
xlabel(h(1),'E (eV)')
set(h(1),'yscale','log');
set(h(1),'xscale','log');
irf_legend(h(1),{'0 deg'},[0.91 0.92],'color','k')
irf_legend(h(1),{'90 deg'},[0.91 0.84],'color','r')
irf_legend(h(1),{'180 deg'},[0.91 0.76],'color','b')