function [elmn,h]=c_ev_mp(vngse,tint,e,b,sc_list);
% [elmn,h]=c_ev_mp(vngse,tint,e,b,sc_list);
% [elmn,h]=c_ev_mp(vngse,tint);
% [elmn,h]=c_ev_mp(vngse,tint,sc_list);
%
% plots electric field  and convection velocities in the magnetopause reference frame
% returns handles to plots and electric field
%
% vngse - magnetopause velocity in GSE
% tint=[t1 t2] - time interval to plot in isdat epoch time
% ic - spacecraft number
% e - electric field in DS ref frame, if not given loaded from mE.mat
% b - b field in DS ref frame, if not given loaded from mB.mat
% elmn = [t El Em En] field in NML reference frame

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'irf_pl_ebv_nrf')

persistent q_flag  q_L_direction

q_flag=av_q('LMN frame defined by \n 0) L || B and N closest to vn \n 1) N || vn (stationary frame), L along mean B \n 2) N||vn, L closest to the specified direction [%]>','q_flag',2);
if q_flag == 0, title_lmn='L-along B, N-closest to vn, M=NxL';flag=0;end
if q_flag == 1, title_lmn='N-along vn, L-mean direction of B, M=NxL';flag=1;end
if q_flag == 2,
  q_L_direction=av_q('L direction (GSE) closest to specified vector [%]>','q_L_direction',[0 0 1]);
  L_dir=av_norm(q_L_direction);
  title_lmn=['N-along vn, L-closest to [' num2str(L_dir,'%7.2f')  '], M=NxL'];
end

if nargin<1, disp('Not enough arguments. See usage:');help c_ev_mp;return; end

if length(vngse) == 3, vngse=[tint(1) vngse];end

if nargin==1 | isempty(tint),   load mE dE1;tint=[E(1,1) E(end,1)]; end
if nargin == 3;  sc_list=e; end
if nargin<3, sc_list=1:4; end
if nargin<=3,
  for ic=sc_list,eval(av_ssub('load mE dE?;de?=av_t_lim(dE?,tint);clear dE?;disp(''..dE?'');',ic)),end
  for ic=sc_list,eval(av_ssub('load mB dB?;db?=av_t_lim(dB?,tint);clear dB?;disp(''..dB?'')',ic)),end
end
if nargin == 4,  sc_list=1;de1=e;db1=b; end
if nargin == 5,  eval(av_ssub('de?=av_t_lim(e,tint);db?=av_t_lim(b,tint);',sc_list));  end

for ic=sc_list, % which satellite

  %%%%%%%%%%%%%%%%%%%%% Convert to MP boundary reference frame %%%%%%%%%%%%
  eval(av_ssub('vn=c_gse2dsc(vngse,[tint(1) ic]);b=db?;e=de?;',ic));
  if q_flag == 2, flag=c_gse2dsc(L_dir,[tint(1) ic]);end
  bpol=av_car2pol(b);b_angle=[bpol(:,1) bpol(:,3)];
  be=av_interp(b,e);
  % make assumption that E.B=0
  [eb,deg]=av_ed(e,be,10);
  % estimate E in boundary system Ev=E+(v x B)
  evxb=av_t_appl(av_cross(be,vn),'*1e-3*(-1)');
  ebv=av_add(1,eb,1,evxb);
  ev=av_add(1,e,1,evxb);
  veb=av_e_vxb(ebv,be,-1);

  %%%%%%%%%%%%%%%%%%%%%%%%   Convert to LMN %%%%%%%%%%%%%%%%%%%%%%%%
  [ev_lmn,L,M,N]=av_c_mp(ev,be,vn,flag);
  eb_lmn=av_c_mp(eb,be,vn,flag);
  ebv_lmn=av_c_mp(ebv,be,vn,flag);
  evxb_lmn=av_c_mp(evxb,be,vn,flag);
  elmn=ebv_lmn;
  b_lmn=av_c_mp(be,be,vn,flag);
  v_lmn=av_c_mp(veb,be,vn,flag);

  %%%%%%%%%%%%%%%%  Get the directions of L,M,N in GSE coordinates %%%%%%%%%%%%%%%%
  L=c_gse2dsc(L,[tint(1) ic],-1);
  M=c_gse2dsc(M,[tint(1) ic],-1);
  N=c_gse2dsc(N,[tint(1) ic],-1);
  if size(L,1)>1,
    L_str=['L is timevarying, L_1=[' num2str(L(1,2:4),'%7.2f') ']'];
  else
    L_str=['L=[' num2str(L(1,2:4),'%7.2f') ']'];  
  end
  if size(M,1)>1,
    M_str=['M is timevarying, M_1=[' num2str(M(1,2:4),'%7.2f') ']'];  
  else
    M_str=['M=[' num2str(M(1,2:4),'%7.2f') ']'];  
  end
  if size(N,1)>1,
    N_str=['N is timevarying, N_1=[' num2str(N(1,2:4),'%7.2f') ']'];  
  else
    N_str=['N=[' num2str(N(1,2:4),'%7.2f') ']'];  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%   PLOT %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(ic);clf
  npl=7;ipl=1;
  h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  if exist('mP.mat'), 
    eval(av_ssub('load mP P? NVps?;p=av_t_lim(P?,tint);n=av_t_lim(NVps?,tint);',ic));
    av_tplot(n);
    ylabel('n_{sc} [cm^{-3}]');
    ax=get(gca,'ylim'); ylim=[0 ax(2)]; set(gca,'ylim',ylim);
    ytick=get(gca,'YTick');
    %set(gca,'YTick',[-35 -30 -25 -20]);
    ax2 = axes('Position',get(gca,'Position'),'YAxisLocation','right','XTick',[],'Color','none');
    ntick=c_n_Vps([0 ytick],-1);ntick=ntick(2:end);
    set(ax2,'YLim',ylim,'YTick',ytick);
    set(ax2,'YTickLabel',num2str(ntick',2));ylabel('V_{ps} [cm^{-3}]');
    title(['sc ' num2str(ic) ' vn_{GSE}=' num2str(av_abs(vngse,1),3) ' [' num2str(av_norm(vngse(1,2:4)),2) '] km/s. ' title_lmn '\newline' L_str ', ' M_str ', ' N_str '.']);
  else
    disp('No satellite potential data');
  end
  av_pl_info([mfilename '  '  datestr(now)]); % add information to the plot
  
  h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(b_lmn);grid on;axis tight;ylabel('B_{LMN} [nT]');
  
  h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot([eb(:,1:2) -eb(:,3:4)]);grid on; ylabel('E_{DSI} [mV/m]');
  
  h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(ev_lmn);grid on; ylabel('E+vxB [mV/m]');
  
  h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(ebv_lmn);grid on; ylabel('E_{E.B=0}+vxB [mV/m]');
  
  h(ic,ipl)=av_subplot(npl,1,-ipl);ipl=ipl+1;
  av_tplot(v_lmn);grid on; ylabel('v_{ExB}-Vmp [km/s]');
  
  h(ic,ipl)=av_subplot(npl,1,ipl);ipl=ipl+1;
  av_tplot(b_angle);grid on;av_zoom([-90 90],'y');title('B elevation angle. +-10 deg limits marked.');ylabel('[degrees]');
  ll=line(b_angle([1 end end 1],1),[10 10 -10 -10],'Color',[.8 .8 .8]);
  
end

av_zoom(tint,'x',h(sc_list,:));
av_zoom([-15 15],'y',h(sc_list,3:5)); % measured E field
av_zoom([-400 400],'y',h(sc_list,6)); % ExB velocity
for ic=sc_list,
  add_timeaxis(h(ic,:));
  legend(h(ic,2),'L','M','N')
  legend(h(ic,3),'X','Y','Z')
  legend(h(ic,4),'L','M','N');
  legend(h(ic,5),'L','M','N');
  legend(h(ic,6),'L','M','N');
end

% 1. Vps
% 2. Blmn
% 3. E DSI
% 4. E lmn
% 5. Elmn_E.B=0
% 6. Vlmn_E.B=0
% 7. angle
