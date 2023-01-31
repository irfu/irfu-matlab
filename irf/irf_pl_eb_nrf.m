function [elmn,h]=irf_pl_eb_nrf(vngse,tint,e,b,sc_list)
%IRF_PL_EB_NRF   Plot E and B in the new reference frame
%
% OLD ROUTINE!!!!!
%
% [elmn,h]=irf_pl_eb_nrf(vngse,tint,e,b,sc_list);
% [elmn,h]=irf_pl_eb_nrf(vngse,tint);
% [elmn,h]=irf_pl_eb_nrf(vngse,tint,sc_list);
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
%

q_flag=irf_ask('LMN frame defined by \n 0) L || B and N closest to vn \n 1) N || vn (stationary frame), L along mean B \n 2) N||vn, L closest to the specified direction [%]>','q_flag',2);
if q_flag == 0, title_lmn='L-along B, N-closest to vn, M=NxL';flag=0;end
if q_flag == 1, title_lmn='N-along vn, L-mean direction of B, M=NxL';flag=1;end
if q_flag == 2
  q_L_direction=irf_ask('L direction closest to specified vector [%]>','q_L_direction',[0 0 1]);
  L_dir=irf_norm(q_L_direction);flag=L_dir;
  title_lmn=['N-along vn, L-closest to [' num2str(L_dir,'%7.2f')  '], M=NxL'];
end

if nargin<1, disp('Not enough arguments. See usage:');help irf_pl_eb_nrf;return; end

if length(vngse) == 3, vngse=[tint(1) vngse];end

if nargin==1 || isempty(tint),   load mE dE1;tint=[E(1,1) E(end,1)]; end
if nargin == 3;  sc_list=e; end
if nargin<3, sc_list=1:4; end
if nargin<=3
  for ic=sc_list,eval(irf_ssub('load mEDSI diE?p1234;die?=irf_tlim(diE?p1234,tint);clear diEp?1234;disp(''..diE?p1234'');',ic)),end
  for ic=sc_list,eval(irf_ssub('load mB diB?;dib?=irf_tlim(diB?,tint);clear diB?;disp(''..diB?'')',ic)),end
end
if nargin == 4,  sc_list=1;die1=e;dib1=b; end %#ok<NASGU>
if nargin == 5,  eval(irf_ssub('die?=irf_tlim(e,tint);dib?=irf_tlim(b,tint);',sc_list));  end

for ic=sc_list % which satellite
  c_eval('vn=c_gse2dsc(vngse,[tint(1) ic],2);b=dib?;e=die?;',ic);
  bpol=av_car2pol(b);b_angle=[bpol(:,1) bpol(:,3)];
  be=irf_resamp(b,e);
  % make assumption that E.B=0
  [eb,deg]=irf_edb(e,be,5);
  % estimate E in boundary system Ev=E+(v x B)
  evxb=irf_tappl(irf_cross(be,vn),'*1e-3*(-1)');
  ebv=irf_add(1,eb,1,evxb);
  ev=irf_add(1,e,1,evxb);
  
  ev_lmn=irf_eb_nrf(ev,be,vn,flag);
  eb_lmn=irf_eb_nrf(eb,be,vn,flag);
  ebv_lmn=irf_eb_nrf(ebv,be,vn,flag);
  evxb_lmn=irf_eb_nrf(evxb,be,vn,flag);
  b_lmn=irf_eb_nrf(b,be,vn,flag);
  enml=ebv_lmn;
  
  figure(ic);clf
  npl=7;ipl=1;
  h(ic,ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  if exist('mP.mat','file'), eval(irf_ssub('load mP P?;p=irf_tlim(P?,tint);',ic));irf_plot(p);end
  title(['sc ' num2str(ic) ' vn_{GSE}=' num2str(irf_abs(vngse,1),3) ' [' num2str(irf_norm(vngse(1,2:4)),2) '] km/s. ' title_lmn]);
  ylabel('Vps [V]');
  irf_pl_info(['c\_e\_mp() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))]); % add information to the plot
  
  irf_zoom([-35 -2],'y');
  
  h(ic,ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(b_lmn);grid on;hold on;
  ylabel('B [nT] LMN');
  
  h(ic,ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(eb);grid on;irf_zoom([-15 15],'y'); ylabel('E [mV/m] DSI');
  
  h(ic,ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(evxb_lmn);irf_zoom([-10 10],'y');grid on; ylabel('Vn x B [mV/m]');
  
  h(ic,ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(ev_lmn);grid on; irf_zoom([-15 15],'y');ylabel('E+vxB [mV/m]');
  %legend('E_l','E_m','E_n');
  
  h(ic,ipl)=irf_subplot(npl,1,-ipl);ipl=ipl+1;
  irf_plot(ebv_lmn);grid on; irf_zoom([-15 15],'y');ylabel('E_{E.B=0}+vxB [mV/m]');
  %legend('E_l','E_m','E_n');
  
  h(ic,ipl)=irf_subplot(npl,1,ipl);ipl=ipl+1;
  irf_plot(b_angle);grid on;irf_zoom([-90 90],'y');title('B elevation angle. +-5 deg limits marked.');ylabel('[degrees]');
  ll=line(b_angle([1 end end 1],1),[5 5 -5 -5],'Color',[.8 .8 .8]);
  
end

irf_zoom(tint,'x',h(sc_list,:));
irf_zoom([-25 -2],'y',h(sc_list,1)); % V_sc
irf_zoom([-30 30],'y',h(sc_list,2)); % B field
irf_zoom([-10 10],'y',h(sc_list,3:6)); % measured E field
for ic=sc_list
  irf_timeaxis(h(ic,:));
  legend(h(ic,2),'B_L','B_M','B_N')
  legend(h(ic,3),'Ex','Ey')
  legend(h(ic,4),'E_L','E_M','E_N');
  legend(h(ic,5),'E_L','E_M','E_N');
  legend(h(ic,6),'E_L','E_M','E_N');
end

