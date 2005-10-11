%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot with wave length, interferometry and particle flux coefficient
%%%%%%%%%%%%%%% p23->p41,p42->p13 %%%%%%%%%%%%%%
%
% In interferometry correlate the maxima, mininma and steepest gradients of 
% both signals. This means taking first and second time derivatives 
% and checking for zero crossings. 
%

mode = irf_ask('Mode p1_34->p_34_2(1) or p23->p41(0)? [%]>','mode',0);

if mode, disp('p1_34->p_34_2, p1_34->p_34_2 interferometry');
else, disp('p23->p41, p42->p13 interferometry');
end

flag_corr_deriv = irf_ask('Use max gradient (0) or zero crossings(1)? [%]>',...
	'flag_corr_deriv',0);

ic=irf_ask('Which s/c? [%]>','ic',2);
ic_str=['s/c ' num2str(ic)];
 
start_time = irf_ask('Start time [y m d h m s] [%]>','start_time',[2004 01 04 12 47 9.4]);

Dt = irf_ask('Time interval in seconds? [%]>','Dt',6);
tint = toepoch(start_time) + [0 Dt];

ff = irf_ask('Frequency interval to filter? [%]>','ff',[10 25]);
fstr = ['filter [' num2str(ff(1)) ' ' num2str(ff(2)) '] Hz'];
 
ff_ref = (ff(1)+ff(2))/2;
disp(['Using f=' num2str(ff_ref) ' Hz to estimate phase velocities']);

% Distance between the signals
if mode, dist=.044;
else, dist=.062;
end

Fs = irf_ask('sampling frequnecy? [%]>','Fs',9000);

kHz = irf_ask('Use P4kHz?px or P32kHz?px? 4/32 [%]>','kHz','4');
if ~strcmp(kHz,'4') & ~strcmp(kHz,'32')
	kHz = '4';
	disp('using 4kHz')
end
c_eval(['load mEFWburst P' kHz 'kHz?p1 P' kHz 'kHz?p2 P' kHz 'kHz?p3 P' kHz 'kHz?p4;'],ic)
c_eval(['p1=P' kHz 'kHz?p1; p2=P' kHz 'kHz?p2; p3=P' kHz 'kHz?p3; p4=P' kHz 'kHz?p4;'],ic)
c_eval(['clear P' kHz 'kHz?p1 P' kHz 'kHz?p2 P' kHz 'kHz?p3 P' kHz 'kHz?p4'],ic)
A = c_load('A?',ic,'var');

if mode
	p1s=[p1(:,1) (-p1(:,2) +(p3(:,2)+p4(:,2))/2)/.044];
	ps2=[p2(:,1) (+p2(:,2) -(p3(:,2)+p4(:,2))/2)/.044];
	ps3=[p3(:,1) (-p3(:,2) +(p1(:,2)+p2(:,2))/2)/.044];
	p4s=[p4(:,1) ( p4(:,2) -(p1(:,2)+p2(:,2))/2)/.044];
else
	p13=[p1(:,1) (p1(:,2)-p3(:,2))/.044/sqrt(2)];
	p42=[p4(:,1) (p4(:,2)-p2(:,2))/.044/sqrt(2)];
	p23=[p2(:,1) (p2(:,2)-p3(:,2))/.044/sqrt(2)];
	p41=[p4(:,1) (p4(:,2)-p1(:,2))/.044/sqrt(2)];
end
p1234=[p4(:,1) (p1(:,2)+p2(:,2)+p3(:,2)+p4(:,2))/4];
n = c_efw_scp2ne(p1234);
nf = irf_filt(n,ff(1),ff(2),Fs,5);
nf = irf_tlim(nf,tint);
n = irf_tlim(n,tint);

e = c_load('dibE?p1234',ic,'var');
ef=irf_filt(e,ff(1),ff(2),Fs,5);
ef = irf_tlim(ef,tint);
e = irf_tlim(e,tint);

disp('...data loaded');

% two directions
bn = irf_ask('N direction DSI? [%]>', 'bn', [1 0 0]);
bn=irf_norm(bn);
bp = irf_ask('M direction DSI? [%]>', 'bp', [0 1 0]);
bp=irf_norm(irf_cross(irf_cross(bn,bp),bn));
leg = {'N','M'};

nfef=irf_vec_x_scal(ef,nf,1); % dE*dn 
nfef_bn=irf_dot(nfef,bn);
nfef_bp=irf_dot(nfef,bp);

if mode
	psignal={'p1s','ps2','ps3','p4s'};
	psignal_f={'p1sf','ps2f','ps3f','p4sf'};
	%psignal_fh={'p1sfh','ps2fh','ps3fh','p4sfh'};
else
	psignal={'p13','p42','p23','p41'};
	psignal_f={'p13f','p42f','p23f','p41f'};
	%psignal_fh={'p13fh','p42fh','p23fh','p41fh'};
end
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

if mode
	[ts3_d,t4s_d,ts3_dd,t4s_dd] = irf_corr_deriv(ps3f,p4sf,flag_corr_deriv);
	%{ 
	vi_d_s3_4s = [ts3_d (t4s_d-ts3_d)/dist];
	vi_dd_s3_4s = [ts3_dd (t4s_dd-ts3_dd)/dist];
	vi_s3_4s = sortrows([vi_d_s3_4s;vi_dd_s3_4s]);
	ii = find(diff(vi_s3_4s(:,1))==0); %Remove repeating points
	vi_s3_4s(ii,:) = [];
	 %}
	vi_d_4s_s3 = [ts3_d (ts3_d-t4s_d)/dist];
	vi_dd_4s_s3 = [ts3_dd (ts3_dd-t4s_dd)/dist];
	vi_4s_s3 = sortrows([vi_d_4s_s3;vi_dd_4s_s3]);
	ii = find(diff(vi_4s_s3(:,1))==0); %Remove repeating points
	vi_4s_s3(ii,:) = [];
	
	[ts2_d,t1s_d,ts2_dd,t1s_dd] = irf_corr_deriv(ps2f,p1sf,flag_corr_deriv);
	vi_d_s2_1s = [ts2_d (t1s_d-ts2_d)/dist];
	vi_dd_s2_1s = [ts2_dd (t1s_dd-ts2_dd)/dist];
	vi_s2_1s = sortrows([vi_d_s2_1s;vi_dd_s2_1s]);
	ii = find(diff(vi_s2_1s(:,1))==0); %Remove repeating points
	vi_s2_1s(ii,:) = [];
	
	%{ 
	[t1,t2,i1,i2] = irf_find_closest(vi_s3_4s(:,1),vi_s2_1s(:,1));
	vi23 = [vi_s3_4s(i1,1) vi_s3_4s(i1,2) vi_s2_1s(i2,2) 0*t1];
	 %}
	
	[t1,t2,i1,i2] = irf_find_closest(vi_4s_s3(:,1),vi_s2_1s(:,1));
	vi23 = [vi_4s_s3(i1,1) vi_s2_1s(i2,2) vi_4s_s3(i1,2)];
	
	k = [vi23(:,1) ff_ref*irf_abs([vi23 vi_4s_s3(i1,1)*0],1)];
	
	vphiDS = c_efw_despin(vi23,A);
	vphiDS(:,3) = -vphiDS(:,3); % convert to DSI
else
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
end

vphi_bn=irf_dot(vphiDS,bn);
vphi_bp=irf_dot(vphiDS,bp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;tref=0;
h=irf_plot({n,n,n,n,nf,[nfef_bn nfef_bp(:,2)]});

%%%%%% subplot 1 %%%%%%
axes(h(1));ylabel('N_{Vps} [cm^{-3}]');
title_text=['s/c' num2str(ic) '  ' legend_corr ff_str '.'];
irf_pl_info([mfilename '  ' datestr(now) '. ' title_text]); % add information to the plot
%irf_zoom([0 2.5],'y',gca);

ud=get(gcf,'userdata');
if isfield(ud,'t_start_epoch'), ts = ud.t_start_epoch;
else, ts = 0;
end

%%%%%% subplot 2 %%%%%%
axes(h(2));
%eval(['plot(' psignal{1} '(:,1)-tref, ' psignal_f{1} '(:,2),' psignal{2} '(:,1)-tref, ' psignal_f{2} '(:,2),' psignal{3} '(:,1)-tref, ' psignal_f{3} '(:,2),' psignal{4} '(:,1)-tref, ' psignal_f{4} '(:,2))' ]);
eval(['plot(' psignal{1} '(:,1)-ts, ' psignal_f{1} '(:,2),' psignal{3} '(:,1)-ts, ' psignal_f{3} '(:,2))']);
ylabel('E_{f} [mV/m]');
irf_zoom([-49.9 49.9],'y');
grid on;

%%%%%% subplot 3 %%%%%%
axes(h(3));
plot(vphi_bn(:,1)-ts,vphi_bn(:,2),'.',vphi_bp(:,1)-ts,vphi_bp(:,2),'.');grid on;
%hold on;plot(vphi_n_sc(:,1)-tref,vphi_n_sc(:,2),vphi_m_sc(:,1)-tref,vphi_m_sc(:,2));grid on;
ylabel('k/\omega [s/km]');
irf_zoom([-0.0199 .0199],'y');

%%%%%% subplot 4 %%%%%%
axes(h(4));
%plot(k(:,1)-tref, k(:,2),'b.',k(ind,1)-tref,k(ind,2),'c.');
plot(k(:,1)-ts, k(:,2),'k.');
ylabel('\lambda^{-1} [1/km]')
irf_zoom([0 4.99],'y');
grid on;

axes(h(5));ylabel('dn [cm^-3]');
%irf_zoom([-0.1999 0.1999],'y');
%axis tight;

axes(h(6));ylabel('dn dE [cc mV/m]');axis tight;

numb={'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
for ip=1:6,
  axes(h(ip));
  ht=irf_pl_info(numb{ip},gca,[0.05,.7]);
  set(ht,'fontsize',10);
end

for j=1:6,irf_zoom(tint,'x',h(j));end
add_timeaxis(h)
legend
irf_figmenu

axes(h(2)); legend(psignal{1},psignal{3});
axes(h(3)); legend(leg{1},leg{2});
axes(h(6)); legend(leg{1},leg{2})
