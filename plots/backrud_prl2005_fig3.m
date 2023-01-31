% read all data
%c_get_batch(toepoch([2003 02 14 14 30 0]),5*60)
%c_get_batch(toepoch([2003 02 14 14 30 0]),5*60,'vars','pburst')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot with wave length, interferometry and particle flux coefficient
%%%%%%%%%%%%%%% p23->p41,p42->p13 %%%%%%%%%%%%%%
%
% In interferometry correlate the maxima, mininma and steepest gradients of
% both signals. This means taking first and second time derivatives
% and checking for zero crossings.

mode = 0; % correlating p23->p41
dist=.062; % distance in km between p23 and p41
disp('p23->p41, p42->p13 interferometry');

%flag_corr_deriv = 0; % Use max gradient for correlation
flag_corr_deriv = 1; % Use zero crossings for correlation

ic=2; % sc number
ic_str=['s/c ' num2str(ic)];

start_time = [2003 02 14 14  32 59.9];
Dt=0.2;
tint = toepoch(start_time) + [0 Dt];

ff = [40 200];% frequency interval in which filter data
fstr = ['filter [' num2str(ff(1)) ' ' num2str(ff(2)) '] Hz'];
ff_ref = (ff(1)+ff(2))/2;
disp(['Using f=' num2str(ff_ref) ' Hz to estimate phase velocities']);

Fs = 9000; % sampling frequency in Hz
kHz = '4'; % 4kHz filter for internal burst

c_eval(['load mEFWburst P' kHz 'kHz?p1 P' kHz 'kHz?p2 P' kHz 'kHz?p3 P' kHz 'kHz?p4;'],ic)
c_eval(['p1=P' kHz 'kHz?p1; p2=P' kHz 'kHz?p2; p3=P' kHz 'kHz?p3; p4=P' kHz 'kHz?p4;'],ic)
c_eval(['clear P' kHz 'kHz?p1 P' kHz 'kHz?p2 P' kHz 'kHz?p3 P' kHz 'kHz?p4'],ic)
c_load diB2;
A = c_load('A?',ic,'var');

p13=[p1(:,1) (p1(:,2)-p3(:,2))/.044/sqrt(2)];
p42=[p4(:,1) (p4(:,2)-p2(:,2))/.044/sqrt(2)];
p23=[p2(:,1) (p2(:,2)-p3(:,2))/.044/sqrt(2)];
p41=[p4(:,1) (p4(:,2)-p1(:,2))/.044/sqrt(2)];

p1234=[p4(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
n = c_efw_scp2ne(p1234);
n(:,2)=n(:,2)./5; % callibrate density dividing by a factor of 5
nf = irf_filt(n,ff(1),ff(2),Fs,5);
nf = irf_tlim(nf,tint);
n = irf_tlim(n,tint);

e = c_load('dibE?p1234',ic,'var');
ef=irf_filt(e,ff(1),ff(2),Fs,5);
ef = irf_tlim(ef,tint);
e = irf_tlim(e,tint);

disp('...data loaded');

bn=irf_norm([diB2(:,1:3) diB2(:,4)*0]);
bp=[0 1 0]; % make bp as close as possible to given vector but still perp to bn
bp=irf_norm(irf_cross(irf_cross(bn,bp),bn));
%bp=irf_cross([0 0 1],bn);
b=irf_abs(irf_tlim(diB2,tint));

efn=irf_dot(ef,bn);
efp=irf_dot(ef,bp);

nfef=irf_vec_x_scal(ef,nf,1); % dE*dn
nfef_bn=irf_dot(nfef,bn); % dE*dn component in bn direction
nfef_bp=irf_dot(nfef,bp); % dE*dn component in bp direction


psignal={'p13','p42','p23','p41'};
psignal_f={'p13f','p42f','p23f','p41f'};

legend_corr=[psignal{1} '->' psignal{2} ' ' psignal{3} '->' psignal{4}];
legend_corr12=[psignal{1} '->' psignal{2}];
legend_corr34=[psignal{3} '->' psignal{4}];

tref=toepoch(start_time);

% Filter the data and Crop the data
ff_str=['f_{filter}=[' num2str(ff(1),3) ' ' num2str(ff(2),3) '] Hz'];
for j=1:4
  eval([psignal_f{j} '=irf_filt(' psignal{j} ',ff(1),ff(2),[],5);'])
  eval([psignal_f{j} '=irf_tlim(' psignal_f{j} ',tint);'])
  eval([psignal{j} '=irf_tlim(' psignal{j} ',tint);'])
end

[t23_d,t41_d,t23_dd,t41_dd]=irf_corr_deriv(p23f,p41f,flag_corr_deriv);
vi_d_23_41=[t23_d (t41_d-t23_d)/dist];
vi_dd_23_41=[t23_dd (t41_dd-t23_dd)/dist];
vi_23_41=sortrows([vi_d_23_41;vi_dd_23_41]);

[t42_d,t13_d,t42_dd,t13_dd]=irf_corr_deriv(p42f,p13f,flag_corr_deriv);
vi_d_42_13=[t42_d (t13_d-t42_d)/dist];
vi_dd_42_13=[t42_dd (t13_dd-t42_dd)/dist];
vi_42_13=sortrows([vi_d_42_13;vi_dd_42_13]);

ii = find(diff(vi_23_41(:,1))==0); %Remove repeating points
vi_23_41(ii,:) = [];
ii = find(diff(vi_42_13(:,1))==0); %Remove repeating points
vi_42_13(ii,:) = [];
[t1,t2,i1,i2]=irf_find_closest(vi_23_41(:,1),vi_42_13(:,1));
vi23=[vi_23_41(i1,1) vi_23_41(i1,2) vi_42_13(i2,2) 0*t1];

k=[vi23(:,1) ff_ref*irf_abs(vi23,1)];

vphiDS=c_despin_new(vi23,A,'23');
vphiDS(:,3) = -vphiDS(:,3); % convert to DSI

vphi_bn=irf_dot(vphiDS,bn);
vphi_bp=irf_dot(vphiDS,bp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;tref=0;
h=irf_plot({irf_tlim(p1234,tint),n,n,n,nf,[nfef_bn nfef_bp(:,2)]});
ipanel=1;npanel=6;

%%%%%% subplot 1 %%%%%%
axes(h(ipanel));ipanel=ipanel+1;
grid on
irf_zoom([-6.9999 -3.0001],'y');
ylabel('V_{ps} [V]');
title_text=['s/c' num2str(ic) '  ' legend_corr ff_str '.'];
%irf_pl_info([mfilename '  ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss")) '. ' title_text]); % add information to the plot

%%%%%% subplot 2 %%%%%%
axes(h(ipanel));ipanel=ipanel+1;
irf_plot(efn,'k');
hold on
irf_plot(efp,'r');
grid on;
ylabel('E_{f} [mV/m]');
irf_zoom([-149.999 149.999],'y');
irf_timeaxis(gca,'nolabels');
legend('|| B','\perp B');


%%%%%% subplot 3 %%%%%%
axes(h(ipanel));ipanel=ipanel+1;
irf_plot(vphi_bn(:,1:2),'k*');
hold on;
irf_plot(vphi_bp(:,1:2),'r*');
grid on;
irf_zoom([-0.0199 .0199],'y');
ylabel('k/\omega [s/km]');
irf_timeaxis(gca,'nolabels');
legend('|| B','\perp B');

%%%%%% subplot 4 %%%%%%
axes(h(ipanel));ipanel=ipanel+1;
irf_plot(k,'k*');
irf_zoom([0 1.999],'y');
grid on;
irf_timeaxis(gca,'nolabels');
ylabel('\lambda^{-1} [1/km]')

%%%%%% subplot 5 %%%%%%
axes(h(ipanel));ipanel=ipanel+1;
irf_plot(nf)
irf_zoom([-1.499 1.499],'y');
grid on
irf_timeaxis(gca,'nolabels');
ylabel('dn [cm^{-3}]');

%%%%%% subplot 6 %%%%%%
axes(h(ipanel));ipanel=ipanel+1;
irf_plot(nfef_bn,'k');
hold on;
irf_plot(nfef_bp,'r');
grid on
irf_zoom([-99.9 99.9],'y');
ylabel('dn dE [cm^{-3} mV/m]')
legend('|| B','\perp B');

irf_figmenu

%irf_zoom(tint,'x',h)
numb={'A','B','C','D','E','F','G','H','I'};
for ip=1:ipanel-1
  axes(h(ip));
  ht=irf_pl_info(numb{ip},gca,[0.02,.8]);
  set(ht,'fontsize',12);
end

