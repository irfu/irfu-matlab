function [vht,corr_coef,eht]=irf_vht_plot(e,b,tint,vht_flag,vht)
%IRF_VHT_PLOT compute and plot de Hoffmann-Teller velocity
%
% [VHT,CORR_COEF] = IRF_VHT_PLOT(E,B,[TINT,VHT_FLAG,VHT])
%
% Make the standard plot for estimate of the goodness of HT-frame
% E vs V_HT x B
%
% e - E field, first column time
% b - B field, first column time, inside code b is interpolated to e
% c - returns the value of correlation coefficient
% flag - see irf_VHT
%        flag == 2 (default), use two E field components for VHT estimate
%        flag == 3, use all 3 components
% vht - VHT velocity vector
%           if given do not calculate VHT but plot with a specified vht
%
% See also IRF_VHT.

if isa(e,'TSeries') && isa(b,'TSeries')
	e=[e.time.epochUnix double(e.data)];
	b=[b.time.epochUnix double(b.data)];
end

if nargin < 3 || isempty(tint)   
  tint=[min([e(1,1),b(1,1)]) max([e(end,1),b(end,1)])];
end
if nargin < 4, vht_flag=2;end
if nargin == 5, vht_is='given'; else, vht_is='calculated';end

if isa(tint,'EpochTT')
  tint = tint.epochUnix';
end
strint = irf_disp_iso_range(tint,1);
%strint = tint.utc;
disp(strint);

e=irf_tlim(e,tint);
b=irf_resamp(b,e);
if strcmp(vht_is,'given')  % VHT is given
    vht=vht(1,end-2:end); % assume that vht=[[t] vx vy vz]; use GSE notation
    eht=irf_e_vxb([0 vht],b); % evht=irf_add(1,e,-1,eht); evht would be E field in VHT frame
    ind_nan=find(isnan(e(:,2) + b(:,2))); % remove NaN in B or E from calculation
    if vht_flag == 2
        ep=[e(ind_nan,2);e(ind_nan,3)];
        ehtp=[eht(ind_nan,2);eht(ind_nan,3)];
    else
        ep=[e(ind_nan,2);e(ind_nan,3);e(ind_nan,4)];
        ehtp=[eht(ind_nan,2);eht(ind_nan,3);eht(ind_nan,4)];
    end
    p=polyfit( ehtp,ep,1);
    cc=corrcoef(ep,ehtp);
    corr_coef=cc(1,2);
elseif strcmp(vht_is,'calculated')
    if vht_flag == 2
        [vht,eht,dvht,p,cc]=irf_vht(e,b,2);
    else
        [vht,eht,dvht,p,cc]=irf_vht(e,b,1);
    end
    corr_coef=cc(1,2);
end

strvht=['V_{HT}=' num2str(irf_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] km/s GSE'];
strdvht=['\Delta V_{HT}=' num2str(irf_abs(dvht,1),3) ' [ ' num2str(irf_norm(dvht),' %5.2f') '] km/s GSE'];

if vht_flag == 2
    plot(eht(:,2),e(:,2),'b.',eht(:,3),e(:,3),'r.');
else
    plot(eht(:,2),e(:,2),'b.',eht(:,3),e(:,3),'r.',eht(:,4),e(:,4),'g.');
end
axis equal;grid on;
title('deHoffmann-Teller frame');
xlabel('E_{HT} [mV/m] DS');ylabel('E [mV/m] DS')
if vht_flag == 2, legend({'x','y'},'location','southeast');
else, legend({'x','y','z'});
end

ax=axis;
ymax=ax(4);ymin=ax(3);dy=(ymax-ymin)/20;
ytext=ymax-dy;
xtext=ax(1)+(ax(2)-ax(1))/20;
text(xtext,ytext,strint);ytext=ytext-dy;
text(xtext,ytext,strvht);ytext=ytext-dy;
text(xtext,ytext,strdvht);ytext=ytext-dy;
hold on
xp=[min(min(eht(:,2:end))) max(max(eht(:,2:end)))];
plot(xp,polyval(p,xp),'k-');
axis(ax);
text(xtext,ytext,['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)]);ytext=ytext-dy;
text(xtext,ytext,['cc=' num2str(cc(1,2),3)]);ytext=ytext-dy;
if strcmp(vht_is,'given')
    text(xtext,ytext,'V_{HT} given as input');
end

