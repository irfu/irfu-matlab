%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
function status = polarplot2_manager_v(arg)
%
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
global ud hf
persistent flag_calculate_v

gcbf=hf;
if gcbf,
  ud = get(gcbf, 'userdata');
else
  par=get(gca,'parent');ud = get(par, 'userdata');
end

switch arg
  case 'ok'
    uiresume;
  case 'recalculate'
    flag_calculate_v=1;
  case 'fpl'
    fplmin=str2num(get(ud.fplminh,'string'));
    fplmax=str2num(get(ud.fplmaxh,'string'));
    av_zoom([fplmin fplmax],'x',ud.h2([1:3]));set(h2(1),'ylim','auto');


end

if flag_calculate_v==1,
  flag_fitzero=get(ud.vfitzero,'Value');
  sampling_distance=str2num(get(ud.sampling_distance_h,'string'));
  fmin=str2num(get(ud.fminh,'string'));
  fmax=str2num(get(ud.fmaxh,'string'));
  freqgood1=ud.freq(ud.ind_coh_good1);
  phgood1=ud.ph1(ud.ind_coh_good1);
  freqgood1=freqgood1(:);  phgood1=phgood1(:);
  ind1=find(freqgood1 < fmin | freqgood1 > fmax);
  freqgood1(ind1)=[];phgood1(ind1)=[];
  freqgood2=ud.freq(ud.ind_coh_good2);
  phgood2=ud.ph2(ud.ind_coh_good2);
  freqgood2=freqgood2(:);  phgood2=phgood2(:);
  ind2=find(freqgood2 < fmin | freqgood2 > fmax);
  freqgood2(ind2)=[];phgood2(ind2)=[];

  if flag_fitzero == 1,    % least square fit
    [vphase1, flag_lsqr_1]=lsqr(freqgood1,phgood1,[],20);
    [vphase2, flag_lsqr_2]=lsqr(freqgood2,phgood2,[],20);
    line1y=[fmin*vphase1 fmax*vphase1];
    line2y=[fmin*vphase2 fmax*vphase2];
  elseif flag_fitzero == 0, % polynomial fit
    p=polyfit(freqgood1,phgood1,1);
    vphase1=p(1);
    line1y=[polyval(p,fmin) polyval(p,fmax)];
    p=polyfit(freqgood2,phgood2,1);
    vphase2=p(1);
    line2y=[polyval(p,fmin) polyval(p,fmax)];
  end
  set(ud.fitline1,'xdata',[fmin fmax],'ydata',line1y);
  set(ud.fitline2,'xdata',[fmin fmax],'ydata',line2y);
  text_vphase1=['v_phase x =' num2str(vphase1,2) ' deg/s, or ' num2str(sampling_distance*360/vphase1,4) ' phys.un./s'];
  text_vphase2=['v_phase y =' num2str(vphase2,2) ' deg/s, or ' num2str(sampling_distance*360/vphase2,4) ' phys.un./s'];
  set(ud.vphase1,'string',text_vphase1);
  set(ud.vphase2,'string',text_vphase2);
  flag_recalculate_v=0;
end
% update the minvar plot
flag_calculate_v=0;
return
X=av_t_lim(ud.X,tlim);
clear ud.Xminvar;
[ud.Xminvar, l, v]=av_minvar(X);
ud.l=l;ud.v=v;ud.v1=v(1,:);ud.v2=v(2,:);ud.v3=v(3,:);
axes(ud.h(2));
av_tplot([ud.Xminvar X(:,5)]);
axis tight;add_timeaxis(ud.h(2),'date');
legend('max','mean','min','abs');
axes(ud.h(3));
plot(ud.Xminvar(:,4),ud.Xminvar(:,2));xlabel('min');ylabel('max');
axis tight;axis equal; ax=axis;grid on;
axes(ud.h(4))
plot(ud.Xminvar(:,3),ud.Xminvar(:,2));xlabel('mean');ylabel('max');
axis equal; grid on;
l_str=['L1=' num2str(l(1),3) ' L2=' num2str(l(2),3) ' L3=' num2str(l(3),3) '\newline'];
v1_str=['v1=[' num2str(v(1,:),'%6.2f') '] \newline'];
v2_str=['v2=[' num2str(v(2,:),'%6.2f') '] \newline'];
v3_str=['v3=[' num2str(v(3,:),'%6.2f') '] \newline'];
v_str=[v1_str v2_str v3_str];
set(ud.result_text,'string',[l_str v_str],'verticalalignment','top');

