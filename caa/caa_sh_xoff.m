function [dE, dAmp, weight] = caa_sh_xoff(st,dt,flag_amp)
%CAA_SH_XOFF  sunward offset and amplitude correction in the sh/sw
%
% [dE, dAmp, weight] = caa_sh_xoff(st,dt [,flag_amp])
% [dE, dAmp, weight] = caa_sh_xoff(iso_st,iso_et [,flag_amp])
%
% Study sunward offset (X GSE) and amplitude correction
% factor by comparing EFW data with CIS HIA
%
% if FLAG_AMP is zero (default), use default amplitude correction
% factor = 1.1
%
% See also CAA_SH_PLAN, CAA_COROF_DSI
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev

STEP = 600; % Averaging window
DEY = 0.5;  % Good Ey correspondence in mV/m
DAMP_DEF = 1.1; % Default amplitude correction factor

dE1 = NaN; dE2 = NaN; dE3 = NaN; dE4 = NaN;
dAmp1 = NaN; dAmp2 = NaN; dAmp3 = NaN; dAmp4 = NaN;
weight = 0;

if nargin<3, flag_amp = 0; end
if flag_amp~=0, flag_amp = 1; end

[st,dt] = irf_stdt(st,dt);

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
h = irf_plot({tt,tt,tt,tt,tt,tt});

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
		if flag_amp, dAmp = find_damp(E1(ii,3), CE1(ii,3));
		else dAmp = DAMP_DEF;
		end
        dEx(:,2) = dAmp*E1(ii,2) - CE1(ii,2);
        irf_plot(dEx,'k'), on = 1;
        dEx1 = mean(dEx(:,2));
        Eref = E1(:,1:3);
        Eref(:,2) = dAmp*Eref(:,2) - dEx1;
		Eref(:,3) = dAmp*Eref(:,3);
    end
end
if ~isempty(E3) && ~isempty(CE3)
    ii = find( abs(CE3(:,3)-E3(:,3)) < DEY );
    if ~isempty(ii)
        dEx = E3(ii,1:2);
		if flag_amp, dAmp = find_damp(E3(ii,3), CE3(ii,3));
		else dAmp = DAMP_DEF;
		end
        dEx(:,2) = dAmp*E3(ii,2) - CE3(ii,2);
        if on, hold on, end
        irf_plot(dEx,'g')
        dEx3 = mean(dEx(:,2));
        if isempty(Eref)
            Eref = E3(:,1:3);
            Eref(:,2) = dAmp*Eref(:,2) - dEx3;
			Eref(:,3) = dAmp*Eref(:,3);
        else
            Eref(:,2) = ( Eref(:,2) + dAmp*E3(:,2) - dEx3 )/2;
			Eref(:,3) = ( Eref(:,3) + dAmp*E3(:,3) )/2;
            irf_log('proc','using two signals')
        end
    end
end


set(gca,'YLimMode','auto', 'XLimMode','auto','XTickLabel','')
xlabel('')
ylabel('dEx [mV/m]')

if ~isempty(Eref)
	weight = length(find(~isnan(Eref(:,2))))/(dt/STEP+1);
    axes(h(5)), cla, irf_plot(Eref(:,[1 2]),'o'), hold on
	axes(h(6)), cla, irf_plot(Eref(:,[1 3]),'o'), hold on
    legx = ''; legy = '';
    if ~isempty(E1)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E1(:,2)) );
		if ~isempty(ii)
			if flag_amp, dAmp1 = find_damp(E1(ii,3), Eref(ii,3));
			else dAmp1 = DAMP_DEF;
			end
			dE1 = mean(dAmp1*E1(ii,2)-Eref(ii,2))/dAmp1;
			axes(h(5)), irf_plot([E1(:,1) dAmp1*E1(:,2)-dE1],'k')
			axes(h(6)), irf_plot([E1(:,1) dAmp1*E1(:,3)],'k')
			legx = num2str(dE1,'dEx1 = %.2f');
			legy = num2str(dAmp1,'dAm1 = %.2f');
		end
    end
    if ~isempty(E2)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E2(:,2)) );
		if ~isempty(ii)
			if flag_amp, dAmp2 = find_damp(E2(ii,3), Eref(ii,3));
			else dAmp2 = DAMP_DEF;
			end
			dE2 = mean(dAmp2*E2(ii,2)-Eref(ii,2))/dAmp2;
			axes(h(5)), irf_plot([E2(:,1) dAmp2*E2(:,2)-dE2],'r')
			axes(h(6)), irf_plot([E2(:,1) dAmp2*E2(:,3)],'r')
			l = num2str(dE2,'dEx2 = %.2f');
			if isempty(legx), legx = l; else legx = [legx ', ' l]; end
			l = num2str(dAmp2,'dAm2 = %.2f');
			if isempty(legy), legy = l; else legy = [legy ', ' l]; end
		end
    end
    if ~isempty(E3)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E3(:,2)) );
		if ~isempty(ii)
			if flag_amp, dAmp3 = find_damp(E3(ii,3), Eref(ii,3));
			else dAmp3 = DAMP_DEF;
			end
			dE3 = mean(dAmp3*E3(ii,2)-Eref(ii,2))/dAmp3; 
			axes(h(5)), irf_plot([E3(:,1) dAmp3*E3(:,2)-dE3],'g')
			axes(h(6)), irf_plot([E3(:,1) dAmp3*E3(:,3)],'g')
			l = num2str(dE3,'dEx3 = %.2f');
			if isempty(legx), legx = l; else legx = [legx ', ' l]; end
			l = num2str(dAmp3,'dAm3 = %.2f');
			if isempty(legy), legy = l; else legy = [legy ', ' l]; end
		end
    end
    if ~isempty(E4)
        ii = find( ~isnan(Eref(:,2)) & ~isnan(E4(:,2)) );
		if ~isempty(ii)
			if flag_amp, dAmp4 = find_damp(E4(ii,3), Eref(ii,3));
			else dAmp4 = DAMP_DEF;
			end
			dE4 = mean(E4(ii,2)-Eref(ii,2));
			axes(h(5)), irf_plot([E4(:,1) dAmp4*E4(:,2)-dE4],'b')
			axes(h(6)), irf_plot([E4(:,1) dAmp4*E4(:,3)],'b')
			l = num2str(dE4,'dEx4 = %.2f');
			if isempty(legx), legx = l; else legx = [legx ', ' l]; end
			l = num2str(dAmp4,'dAm4 = %.2f');
			if isempty(legy), legy = l; else legy = [legy ', ' l]; end
		end
    end
	if ~isempty(CE1)
		axes(h(5)), irf_plot(CE1(:, [1 2]),'k+')
		axes(h(6)), irf_plot(CE1(:, [1 3]),'k+')
	end
	if ~isempty(CE3)
		axes(h(5)), irf_plot(CE3(:, [1 2]),'g+')
		axes(h(6)), irf_plot(CE3(:, [1 3]),'g+')
	end
    hold off
	title(h(5),legx), ylabel(h(5),'Ex [mV/m]')
	if flag_amp, title(h(6),legy), end, ylabel(h(6),'Ey [mV/m]')
end

axes(h(4)), cla
c_pl_tx('Ps?')
ylabel('Sc pot [-V]')
if ~isempty(Eref)
	title(sprintf('%d reference points (%d%% data coverage)',...
		length(find(~isnan(Eref(:,2)))), round(weight*100)))
end

irf_zoom(st+[0 dt],'x',h)
axes(h(6)), add_timeaxis

dE = [dE1 dE2 dE3 dE4];
dAmp = [dAmp1 dAmp2 dAmp3 dAmp4];

function res = find_damp(Ey,CEy)
% find amlitude correction factor by searching for 
% minimum( std( ECISy -dAMP*Ey ) )

damp = 1:0.025:1.4;
dstd = damp;
for i=1:length(damp), dstd(i)=std(CEy - damp(i)*Ey); end
res = damp(dstd==min(dstd));