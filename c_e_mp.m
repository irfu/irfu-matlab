function [elmn,h]=c_e_mp(vngse,tint,e,b,sc_list);
% [elmn,h]=c_e_mp(vngse,tint,e,b,sc_list);
% [elmn,h]=c_e_mp(vngse,tint);
% [elmn,h]=c_e_mp(vngse,tint,sc_list);
%
% plots electric field in the magnetopause reference frame
% returns handles to plots and electric field
%
% vngse - magnetopause velocity in GSE
% tint=[t1 t2] - time interval to plot in isdat epoch time
% ic - spacecraft number
% e - electric field in DS ref frame, if not given loaded from mE.mat
% b - b field in DS ref frame, if not given loaded from mB.mat
% elmn = [t El Em En] field in NML reference frame

q_flag=av_q('LMN frame defined by \n 0) L || B and N closest to vn \n 1) N || vn (stationary frame), L along mean B \n 2) N||vn, L closest to the specified direction [%]>','q_flag',2);
if q_flag == 0, title_lmn='L-along B, N-closest to vn, M=NxL';flag=0;end
if q_flag == 1, title_lmn='N-along vn, L-mean direction of B, M=NxL';flag=1;end
if q_flag == 2,
  q_L_direction=av_q('L direction closest to specified vector [%]>','q_L_direction',[0 0 1]);
  L_dir=av_norm(q_L_direction);flag=L_dir;
  title_lmn=['N-along vn, L-closest to [' num2str(L_dir,'%7.2f')  '], M=NxL'];
end

if nargin<1, disp('Not enough arguments. See usage:');help c_e_mp;return; end

if length(vngse) == 3, vngse=[tint(1) vngse];end

if nargin==1 | isempty(tint),   load mE dE1;tint=[E(1,1) E(end,1)]; end
if nargin == 3;  sc_list=e; end
if nargin<3, sc_list=1:4; end
if nargin<=3,
  for ic=sc_list,eval(av_ssub('load mEDSI diE?p1234;die?=av_t_lim(diE?p1234,tint);clear diEp?1234;disp(''..diE?p1234'');',ic)),end
  for ic=sc_list,eval(av_ssub('load mB diB?;dib?=av_t_lim(diB?,tint);clear diB?;disp(''..diB?'')',ic)),end
end
if nargin == 4,  sc_list=1;die1=e;dib1=b; end
if nargin == 5,  eval(av_ssub('die?=av_t_lim(e,tint);dib?=av_t_lim(b,tint);',sc_list));  end

for ic=sc_list, % which satellite
c_eval('vn=c_gse2dsc(vngse,[tint(1) ic],2);b=dib?;e=die?;',ic);
bpol=av_car2pol(b);b_angle=[bpol(:,1) bpol(:,3)];
be=av_interp(b,e);
% make assumption that E.B=0
[eb,deg]=av_ed(e,be,5);
% estimate E in boundary system Ev=E+(v x B)
evxb=av_t_appl(av_cross(be,vn),'*1e-3*(-1)');
ebv=av_add(1,eb,1,evxb);
ev=av_add(1,e,1,evxb);

ev_lmn=av_c_mp(ev,be,vn,flag);
eb_lmn=av_c_mp(eb,be,vn,flag);
ebv_lmn=av_c_mp(ebv,be,vn,flag);
evxb_lmn=av_c_mp(evxb,be,vn,flag);
b_lmn=av_c_mp(b,be,vn,flag);
enml=ebv_lmn;

figure(ic);clf
npl=7;ipl=1;
h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
if exist('mP.mat'), eval(av_ssub('load mP P?;p=av_t_lim(P?,tint);',ic));av_tplot(p);end
title(['sc ' num2str(ic) ' vn_{GSE}=' num2str(av_abs(vngse,1),3) ' [' num2str(av_norm(vngse(1,2:4)),2) '] km/s. ' title_lmn]);
ylabel('Vps [V]');
av_pl_info(['c\_e\_mp() ' datestr(now)]); % add information to the plot

av_zoom([-35 -2],'y');

h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
av_tplot(b_lmn);grid on;hold on;
ylabel('B [nT] LMN');

h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
av_tplot(eb);grid on;av_zoom([-15 15],'y'); ylabel('E [mV/m] DSI');

h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
av_tplot(evxb_lmn);av_zoom([-10 10],'y');grid on; ylabel('Vn x B [mV/m]');

h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
av_tplot(ev_lmn);grid on; av_zoom([-15 15],'y');ylabel('E+vxB [mV/m]');
%legend('E_l','E_m','E_n');

h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
av_tplot(ebv_lmn);grid on; av_zoom([-15 15],'y');ylabel('E_{E.B=0}+vxB [mV/m]');
%legend('E_l','E_m','E_n');

h(ic,ipl)=av_subplot(npl,1,ipl);ipl=ipl+1;
av_tplot(b_angle);grid on;av_zoom([-90 90],'y');title('B elevation angle. +-5 deg limits marked.');ylabel('[degrees]');
ll=line(b_angle([1 end end 1],1),[5 5 -5 -5],'Color',[.8 .8 .8]);

end

av_zoom(tint,'x',h(sc_list,:));
av_zoom([-25 -2],'y',h(sc_list,1)); % V_sc
av_zoom([-30 30],'y',h(sc_list,2)); % B field
av_zoom([-10 10],'y',h(sc_list,3:6)); % measured E field
for ic=sc_list,
add_timeaxis(h(ic,:));
legend(h(ic,2),'B_L','B_M','B_N')
legend(h(ic,3),'Ex','Ey')
legend(h(ic,4),'E_L','E_M','E_N');
legend(h(ic,5),'E_L','E_M','E_N');
legend(h(ic,6),'E_L','E_M','E_N');
end

