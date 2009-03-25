function [vht,corr_coef]=irf_vht_plot(e,b,tint,vht_flag,vht)
%IRF_VHT_PLOT computer and plot de Hoffman-Teller velocity
%
% [VHT,CORR_COEF] = IRF_VHT_PLOT(E,B,[TINT,VHT_FLAG,VHT])
%
% Make the standard plot for estimate of the goodness of HT-frame
% E vs V_HT x B
%
% e - E field, first column time
% b - B field, first column time, inside code b is interpolated to e 
% c - returns the value of correlation coefficient
% flag - see irf_VHT, default 2 (use two E field components for VHT estimate)
%                    if flag ==3 use all 3 components
% vht - VHT velocity vector 
%           if given do not calculate VHT but plot with a specified vht
%
% See also IRF_VHT.

if nargin < 3 || isempty(tint)
   tint=[min([e(1,1),b(1,1)]) max([e(end,1),b(end,1)])];
end
if nargin < 4, vht_flag=2;end
if nargin == 5, vht_is='given'; else vht_is='calculated';end

strint = irf_disp_iso_range(tint,1);
disp(strint);

e=irf_tlim(e,tint);
b=irf_resamp(b,e);
if strcmp(vht_is,'given'),
  vht=vht(1,end-2:end); % assume that vht=[[t] vx vy vz]; use GSE notation
elseif strcmp(vht_is,'calculated');
  if vht_flag == 2,
    vht=irf_vht(e,b,2);
  else
    vht=irf_vht(e,b,1);
  end
end

strvht=['V_{HT}=' num2str(irf_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] km/s GSE'];
eht=irf_e_vxb([0 vht],b); % evht=irf_add(1,e,-1,eht); evht would be E field in VHT frame

if vht_flag == 2,
  ep=[e(:,2);e(:,3)];xp=[min(ep) max(ep)];
  ehtp=[eht(:,2);eht(:,3)];
else
  ep=[e(:,2);e(:,3);e(:,4)];xp=[min(ep) max(ep)];
  ehtp=[eht(:,2);eht(:,3);eht(:,4)];
end
p=polyfit( ehtp,ep,1);
cc=corrcoef(ep,ehtp);
corr_coef=cc(1,2);

if vht_flag == 2,
   plot(eht(:,2),e(:,2),'b.',eht(:,3),e(:,3),'r.');
else
   plot(eht(:,2),e(:,2),'b.',eht(:,3),e(:,3),'r.',eht(:,4),e(:,4),'g.');
end
axis equal;grid on;
ht=irf_pl_info([mfilename ' ' datestr(now)]); set(ht,'interpreter','none','FontSize', 5);
title('deHoffmann-Teller frame');
xlabel('E_{HT} [mV/m] DS');ylabel('E [mV/m] DS')
if vht_flag == 2, legend('x','y');
else legend('x','y','z');
end

ax=axis;
ymax=ax(4);ymin=ax(3);dy=(ymax-ymin)/20;
ytext=ymax-dy;
xtext=ax(1)+(ax(2)-ax(1))/10;
text(xtext,ytext,strint);ytext=ytext-dy;
text(xtext,ytext,strvht);ytext=ytext-dy;
hold on
plot(xp,polyval(p,xp),'k-');
text(xtext,ytext,['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)]);ytext=ytext-dy;
text(xtext,ytext,['cc=' num2str(cc(1,2),3)]);ytext=ytext-dy;
if strcmp(vht_is,'given'),
   text(xtext,ytext,'V_{HT} given as input');
end

