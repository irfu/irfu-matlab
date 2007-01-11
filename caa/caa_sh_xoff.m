function [dE1,dE2,dE3,dE4] = caa_sh_xoff(st,dt)
%CAA_SH_XOFF  sunward offset in the magnetosheath/sw
%
% caa_sh_xoff(st,dt)
%
% Study sunward offset (X GSE) by comparing EFW data with CIS HIA
%
% See also CAA_SH_PLAN
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev

STEP = 600; % Averaging window
DEY = 0.5;  % Good Ey correspondence in mV/m

dE1 = NaN; dE2 = NaN; dE3 = NaN; dE4 = NaN;

dt = ceil(dt/STEP)*STEP;
t = st:STEP:st+dt;
t = t';

Ps1 = caa_get(st,dt,1,'Ps?');
Ps2 = caa_get(st,dt,2,'Ps?');
Ps3 = caa_get(st,dt,3,'Ps?');
Ps4 = caa_get(st,dt,4,'Ps?');

diEs1 = caa_get(st,dt,1,'diEs?p34');
if ~isempty(diEs1)
    diEs1(isnan(diEs1(:,2)),:)=[];
    E1=irf_resamp(diEs1,t);
else E1 = [];
end

diVCEh1 = caa_get(st,dt,1,'diVCEh?');
if ~isempty(diVCEh1)
    diVCEh1(isnan(diVCEh1(:,2)),:)=[];
    CE1=irf_resamp(diVCEh1,t);
else CE1 = [];
end

diEs2 = caa_get(st,dt,2,'diEs?p34');
if ~isempty(diEs2)
    diEs2(isnan(diEs2(:,2)),:)=[];
    E2 = irf_resamp(diEs2,t);
else E2 = [];
end

diEs3 = caa_get(st,dt,3,'diEs?p34');
if ~isempty(diEs3)
    diEs3(isnan(diEs3(:,2)),:)=[];
    E3 = irf_resamp(diEs3,t);
else E3 = [];
end

diVCEh3 = caa_get(st,dt,3,'diVCEh?');
if ~isempty(diVCEh3)
    diVCEh3(isnan(diVCEh3(:,2)),:)=[];
    CE3=irf_resamp(diVCEh3,t);
else CE3 = [];
end

diEs4 = caa_get(st,dt,4,'diEs?p34');
if ~isempty(diEs4)
    diEs4(isnan(diEs4(:,2)),:)=[];
    E4 = irf_resamp(diEs4,t);
else E4 = [];
end

tt = [1.0097e+09 1 ;1.0097e+09+1 1 ;];
h = irf_plot({tt,tt,tt,tt,tt});

on = 0;
for co=1:2
    axes(h(co)), cla
    if ~isempty(E1), irf_plot(E1(:, [1 1+co]),'k'), on = 1; end
    if on, hold on, end
    if ~isempty(E2), irf_plot(E2(:, [1 1+co]),'r'), on = 1; end
    if on, hold on, end
    if ~isempty(E3), irf_plot(E3(:, [1 1+co]),'g'), on = 1; end
    if on, hold on, end
    if ~isempty(E4), irf_plot(E4(:, [1 1+co]),'b'), on = 1; end
    if on, hold on, end
    if ~isempty(CE1), irf_plot(CE1(:, [1 1+co]),'k+'), on = 1; end
    if on, hold on, end
    if ~isempty(CE3), irf_plot(CE3(:, [1 1+co]),'g+'), on = 1; end
    hold off
    add_timeaxis
    set(gca,'YLimMode','auto', 'XLimMode','auto','XTickLabel','')
    xlabel('')
    if co==1
        ylabel('Ex [mV/m]')
        title(['SH/SW ' epoch2iso(st,1) ' -- ' epoch2iso(st+dt,1)...
            '  EFW (--), CIS HIA (+)'])
    else ylabel('Ey [mV/m]')
    end
end

Eref = [];
axes(h(3)), cla
on = 0;
if ~isempty(E1) && ~isempty(CE1)
    ii = find( abs(CE1(:,3)-E1(:,3)) < DEY );
    if ~isempty(ii)
        dEx = E1(ii,1:2);
        dEx(:,2) = E1(ii,2) - CE1(ii,2);
        irf_plot(dEx,'k'), on = 1;
        dEx1 = mean(dEx(:,2));
    end
    Eref = E1(:,1:2);
    Eref(:,2) = Eref(:,2) - dEx1;
end
if ~isempty(E3) && ~isempty(CE3)
    ii = find( abs(CE3(:,3)-E3(:,3)) < DEY );
    if ~isempty(ii)
        dEx = E3(ii,1:2);
        dEx(:,2) = E3(ii,2) - CE3(ii,2);
        if on, hold on, end
        irf_plot(dEx,'g')
        dEx3 = mean(dEx(:,2));
        if isempty(Eref)
            Eref = E3(:,1:2);
            Eref(:,2) = Eref(:,2) - dEx3;
        else
            Eref(:,2) = (Eref(:,2) + E3(:,2) - dEx3)/2;
            irf_log('proc','using two signals')
        end
    end
end


set(gca,'YLimMode','auto', 'XLimMode','auto','XTickLabel','')
xlabel('')
ylabel('dEx [mV/m]')

axes(h(4)), cla
c_pl_tx('Ps?')
set(gca,'XTickLabel','')
xlabel('')
ylabel('Sc pot [-V]')

axes(h(5)), cla
if ~isempty(Eref)
    irf_plot(Eref,'o'), hold on
    leg = '';
    if ~isempty(E1)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E1(:,2)) );
        dE1 = mean(E1(ii,2)-Eref(ii,2)); 
        irf_plot([E1(:,1) E1(:,2)-dE1],'k')
        leg = num2str(dE1,'dEx1 = %.2f');
    end
    if ~isempty(E2)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E2(:,2)) );
        dE2 = mean(E2(ii,2)-Eref(ii,2)); 
        irf_plot([E2(:,1) E2(:,2)-dE2],'r')
        l = num2str(dE2,'dEx2 = %.2f');
        if isempty(leg), leg = l; else leg = [leg ', ' l]; end
    end
    if ~isempty(E3)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E3(:,2)) );
        dE3 = mean(E3(ii,2)-Eref(ii,2)); 
        irf_plot([E3(:,1) E3(:,2)-dE3],'g')
        l = num2str(dE3,'dEx3 = %.2f');
        if isempty(leg), leg = l; else leg = [leg ', ' l]; end
    end
    if ~isempty(E4)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E4(:,2)) );
        dE4 = mean(E4(ii,2)-Eref(ii,2)); 
        irf_plot([E4(:,1) E4(:,2)-dE4],'b')
        l = num2str(dE4,'dEx4 = %.2f');
        if isempty(leg), leg = l; else leg = [leg ', ' l]; end
    end
    if ~isempty(CE1), irf_plot(CE1(:, [1 2]),'k+'), end
    if ~isempty(CE3), irf_plot(CE3(:, [1 2]),'g+'), end
    hold off
end
title(leg), ylabel('Ex [mV/m]')
irf_zoom(st+[0 dt],'x',h)

if nargout<=1, dE1 = [dE1 dE2 dE3 dE4]; end