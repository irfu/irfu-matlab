function [starttime1,endtime1,starttime3,endtime3] = mms_getprobefields(Exyz,Bxyz,SCpot,zphase) 
%
% [starttime1,endtime1,starttime3,endtime3] = mms_getprobefields(Exyz,Bxyz,SCpot,zphase)
%
% plot probe angles and fields data on MMS and print times when probes satisfy 
% alignment conditions used in Graham et al., JGR, 2015. Used to identify
% times when field-aligned waves can be characterized using interferomtry
% techniques. Can view time delays between electric fields aligned with B.
% Currently p5-p6 are not used in this routine; the analysis is the same as
% the one using Cluster.
% Written by D. B. Graham.
%
% Panels are: (a) B in DMPA Coordinates, (b) Magnitude of B in and out of
% the spin plane, (c) Angles between B and probes 1 and 3 in the spin
% plane (angle between 0 and 90 degrees), (d) Spacecraft potential from
% probes perpendicular to B, (e) E fields from p1-p4 and SC for probes
% closely aligned with B, (f) E in field-aligned coordinates, (g) E from
% probes p1-p2 and p3-p4.
% 
% Input: (only works with TSeries format data.)
%       Exyz -   Electric field in DSL coordinates, brst mode.
%       Bxyz -   Magnetic field in DMPA coordinates.
%       SCpot -  L2 SCpot data. 
%       zphase - Spacecraft phase (zphase). Obtained from ancillary_defatt.
%       e.g. zphase = mms.db_get_variable('mms3_ancillary_defatt','zphase',tint);
%            zphase = irf.ts_scalar(zphase.time, zphase.zphase);
%
% Outputs: Start and end times of intervals which satisfy the probe
% alignment conditions, for probe combinates p1-p2 and p3-p4. 
% 
%
% Temporary fix for axial field direction is used before running this function: 
% Exyz.data(:,3) = -Exyz.data(:,3);

tlimit = irf.tint(SCpot.time.start.utc,SCpot.time.stop.utc);
tlimitl = tlimit+[-10 10];

Bxyz = Bxyz.tlim(tlimitl);
Bxyz = Bxyz.resample(SCpot);
Exyz = Exyz.resample(SCpot);

zphase = zphase.tlim(tlimitl);

nph = length(zphase.data);
for ii=[2:nph]
    if(zphase.data(ii) < zphase.data(ii-1));
        zphase.data(ii:end) = zphase.data(ii:end)+double(360.0);
    end
end

zphase = TSeries(zphase.time,zphase.data,'to',1);
zphase = zphase.resample(SCpot);

%Perform rotation on Exyz into field-aligned coordinates 
SCpos = [0 1 0];

Bmag = Bxyz.abs.data;
Rpar = Bxyz.data./[Bmag Bmag Bmag];
Rperpy = irf_cross(Rpar,SCpos);
Rmag   = irf_abs(Rperpy,1);
Rperpy = Rperpy./[Rmag Rmag Rmag];
Rperpx = irf_cross(Rperpy, Rpar);
Rmag   = irf_abs(Rperpx,1);
Rperpx = Rperpx./[Rmag Rmag Rmag];

Epar = dot(Rpar,Exyz.data,2);
Eperp = dot(Rperpx,Exyz.data,2);
Eperp2 = dot(Rperpy,Exyz.data,2);

Efac = TSeries(Exyz.time,[Eperp Eperp2 Epar],'to',1);

% Probe angles in DSL or whatever
phase_p1=zphase.data/180*pi + pi/6;
phase_p3=zphase.data/180*pi + 2*pi/3;
phase_p2=zphase.data/180*pi + 7*pi/6;
phase_p4=zphase.data/180*pi + 5*pi/3;
rp1=[60*cos(phase_p1) 60*sin(phase_p1)]; 
rp2=[60*cos(phase_p2) 60*sin(phase_p2)];
rp3=[60*cos(phase_p3) 60*sin(phase_p3)];
rp4=[60*cos(phase_p4) 60*sin(phase_p4)];

% Calculate angles between probes and direction of B in the spin plane.
thetap1b = (rp1(:,1).*Bxyz.data(:,1)+rp1(:,2).*Bxyz.data(:,2))./(sqrt(rp1(:,1).^2+rp1(:,2).^2).*sqrt(Bxyz.data(:,1).^2+Bxyz.data(:,2).^2));
thetap1b = acosd(abs(thetap1b));
thetap3b = (rp3(:,1).*Bxyz.data(:,1)+rp3(:,2).*Bxyz.data(:,2))./(sqrt(rp3(:,1).^2+rp3(:,2).^2).*sqrt(Bxyz.data(:,1).^2+Bxyz.data(:,2).^2));
thetap3b = acosd(abs(thetap3b));

SCV12 = (SCpot.data(:,1)+SCpot.data(:,2))/2;
SCV34 = (SCpot.data(:,3)+SCpot.data(:,4))/2;
E1 = (SCpot.data(:,1)-SCV34)*1e3/60;
E2 = (SCV34-SCpot.data(:,2))*1e3/60;
E3 = (SCpot.data(:,3)-SCV12)*1e3/60;
E4 = (SCV12-SCpot.data(:,4))*1e3/60;

E12 = (SCpot.data(:,1)-SCpot.data(:,2))*1e3/120;
E34 = (SCpot.data(:,3)-SCpot.data(:,4))*1e3/120;

idxB = find(sqrt(Bxyz.data(:,1).^2+Bxyz.data(:,2).^2) < abs(Bxyz.data(:,3)));
thresang = 25.0;

E1(find(thetap1b > thresang)) = NaN;
E2(find(thetap1b > thresang)) = NaN;
E3(find(thetap3b > thresang)) = NaN;
E4(find(thetap3b > thresang)) = NaN;
SCV12(find(thetap3b > thresang)) = NaN;
SCV34(find(thetap1b > thresang)) = NaN;
E1(idxB) = NaN;
E2(idxB) = NaN;
E3(idxB) = NaN;
E4(idxB) = NaN;
SCV12(idxB) = NaN;
SCV34(idxB) = NaN;

Bplane = TSeries(Bxyz.time,[sqrt(Bxyz.data(:,1).^2+Bxyz.data(:,2).^2) abs(Bxyz.data(:,3))]);
thetas = TSeries(zphase.time,[thetap1b thetap3b]);
Scpots1234 = TSeries(SCpot.time,[SCV12 SCV34]);
E1234 = TSeries(SCpot.time,[E1 E2 E3 E4]);
Epp1234 = TSeries(SCpot.time,[E12 E34]);

h=irf_plot(7,'newfigure');

h(1)=irf_panel('Bdmpa');
irf_plot(h(1),Bxyz);
ylabel(h(1),'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(h(1),{'B_x','B_y','B_z'},[0.98 0.1])
irf_legend(h(1),'(a)',[0.99 0.98],'color','k')

h(2)=irf_panel('Bdmpa2');
irf_plot(h(2),Bplane);
ylabel(h(2),'B_{DMPA} (nT)','Interpreter','tex');
irf_zoom(h(2),'y',[0,100])
irf_legend(h(2),{'|B_{plane}|','|B_z|'},[0.98 0.1])
irf_legend(h(2),'(b)',[0.99 0.98],'color','k')

h(3)=irf_panel('theta');
irf_plot(h(3),thetas)
ylabel(h(3),'\theta (deg)','Interpreter','tex');
irf_legend(h(3),{'\theta_{1}','\theta_{3}'},[0.98 0.1])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k')
irf_zoom(h(3),'y',[0 90])

h(4)=irf_panel('Vp');
irf_plot(h(4),Scpots1234)
ylabel(h(4),'Pot (V)');
irf_legend(h(4),{'V_{SC12}','V_{SC34}'},[0.98 0.1])
irf_legend(h(4),'(d)',[0.99 0.98],'color','k')

h(5)=irf_panel('Ep');
irf_plot(h(5),E1234)
ylabel(h(5),'E (mV/m)');
irf_legend(h(5),{'E_{SC-p1}','E_{p2-SC}','E_{SC-p3}','E_{p4-SC}'},[0.98 0.1])
irf_legend(h(5),'(e)',[0.99 0.98],'color','k')

h(6)=irf_panel('Efac');
irf_plot(h(6),Efac)
ylabel(h(6),'E (mV/m)');
irf_legend(h(6),{'E_{\perp 1}','E_{\perp 2}','E_{||}'},[0.98 0.1])
irf_legend(h(6),'(f)',[0.99 0.98],'color','k')

h(7)=irf_panel('Eplane');
irf_plot(h(7),Epp1234)
ylabel(h(7),'E (mV/m)');
irf_legend(h(7),{'E_{p2-p1}','E_{p4-p3}'},[0.98 0.1])
irf_legend(h(7),'(g)',[0.99 0.98],'color','k')

irf_plot_axis_align(7)
irf_zoom(h,'x',tlimit);
irf_timeaxis(h);

%set(gcf,'paperpositionmode','auto')
%print('-dpng','-painters','-r600','ESWprobefields.png');

%Find and print times when the probes satisfy alignment conditions
E1temp = isnan(E1)-[-1; isnan(E1(3:length(E1))); 1];
starttime1 = SCpot.time(find(E1temp == 1));
endtime1 = SCpot.time(find(E1temp == -1));

fprintf('Intervals when Probe 1 and 2 satisfy alignment conditions \n')
for q = 1:length(starttime1(:,1))
    fprintf(strcat('No.',num2str(q),' : ',starttime1(q).toUtc,'/',endtime1(q).toUtc,'\n'))
end

E3temp = isnan(E3)-[-1; isnan(E3(3:length(E3))); 1];
starttime3 = SCpot.time(find(E3temp == 1),1);
endtime3 = SCpot.time(find(E3temp == -1),1);

fprintf('Intervals when Probe 3 and 4 satisfy alignment conditions \n')
for q = 1:length(starttime3(:,1))
    fprintf(strcat('No.',num2str(q),': ',starttime3(q).toUtc,'/',endtime3(q).toUtc,'\n'))
end

end
    