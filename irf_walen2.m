function [slope,cc]=irf_walen2(v,b,n,tpar,tperp,vht,tint);
%IRF_WALEN2   Walen test and estimate its goodness
%
% [slope,cc]=irf_walen2(v,b,n,tpar,tperp,vht,tint)
%
% INPUT (all vectors in the same coordinate system e.g GSE or GSM)
% v - ion velocity [time vx vy vz] km/s
% b - B field [time bx by bz] nT,inside code b is interpolated to v 
% n - density [time n] cm-3,inside code n is interpolated to v
% tpar - parallel temperature [time tpar] MK,inside code temp is interpolated to v
% tpar - perpendicular temperature [time tperp] MK,inside code temp is
%interpolated to v
% vht - HT velocity vector [vhtx vhty vhtz] km/s
% tint - toepoch([yyyy mm dd hh mm ss]) + [0 ss]
%
% OUTPUT
% [slope, cc] - slope and correlation coefficient of the Walen plot 
%
% See also IRF_WALEN
%
% $Id$

% Alessandro Retino`
% 2004/11/26


if time interval is not given define it from the size of the input vectors
if nargin < 7, % if tint is not given
   tint=[min([v(1,1),b(1,1),n(1,1)]) max([v(end,1),b(end,1),n(end,1)])];
elseif isempty(tint), % if tint=[], define tint from v and b time axis
   tint=[min([v(1,1),b(1,1),n(1,1)]) max([v(end,1),b(end,1),n(end,1)])];  
end

%display time interval
strint=[epoch2iso(tint(1)) ' -- ' epoch2iso(tint(2)) ];
disp(strint);
%display HT velocity
%strvht=['V_{HT}=' num2str(av_abs(vht,1),3) ' [ ' num2str(irf_norm(vht),' %5.2f') '] km/s GSE'];
%disp(strvht);


%define common time interval for input vectors
n=irf_tlim(n,tint);

irf_plot(n)

%interpolate b,n,tpat,tperp to v

v=av_interp(v,n);
%tpar=av_interp(tpar,n);
%tperp=av_interp(tperp,n);
b=av_interp(b,n);

%
n= [n repmat(n(:,2),1,2)];
%tpar= [tpar repmat(tpar(:,2),1,2)];
%tperp= [tperp repmat(tperp(:,2),1,2)];




%calculate velocity in HT frame, Alfven velocity and pressure anisotropy
vtransf(:,1)=n(:,1);
vtransf(:,2:4)=v(:,2:4)-repmat(vht,size(n,1),1);

alpha(:,1)=n(:,1);
%alpha(:,2)=17.33*n(:,2).*( tpar(:,2)-tperp(:,2) )./( b(:,2).^2+b(:,3).^2+b(:,4).^2);
alpha(:,2)=0.4;
alpha= [alpha repmat(alpha(:,2),1,2)];


valfv(:,1)=n(:,1);

%valfv(:,2:4)=22*b(:,2:4)./ sqrt(n(:,2:4));
valfv(:,2:4)=22*b(:,2:4).*sqrt( 1 - alpha(:,2:4) )./ sqrt( n(:,2:4));

irf_plot(valfv)
irf_plot({alpha,n,valfv})


%no anisotropy
%valfv(:,2:4)=22*b(:,2:4)./ sqrt( n(:,2:4));




%calculate the correlation coefficient between all three components of
%vtransf and valfv


vtransf3=vtransf(:,2:4);
valfv3=valfv(:,2:4);


vtransf3=vtransf3';
valfv3=valfv3';


vtransftot=vtransf3(:);
valfvtot=valfv3(:);





corr=corrcoef(vtransftot,valfvtot);
cc=corr(1,2)
[p,s]=polyfit(valfvtot,vtransftot,1);



p(1)
p(2)



plot(valfv(:,2),vtransf(:,2),'b.',valfv(:,3),vtransf(:,3),'g.',valfv(:,4),vtransf(:,4),'r.');

axis equal;grid on;

ht=irf_pl_info([mfilename ' ' datestr(now)]); set(ht,'interpreter','none','FontSize', 5);

title(['Walen test']);
xlabel('V_{A} [km/s] GSE');ylabel('V-V_{HT} [km/s] GSE')
legend('x','y','z');





ax=axis;
ymax=ax(4);ymin=ax(3);dy=(ymax-ymin)/20;
ytext=ymax-dy;
xtext=ax(1)+(ax(2)-ax(1))/40;
text(xtext,ytext,strint);ytext=ytext-dy;
text(xtext,ytext,strvht);ytext=ytext-dy;
hold on



plot(xp,polyval(p,xp),'k-');
text(xtext,ytext,['slope=' num2str(p(1),3) '  offs=' num2str(p(2),2)]);ytext=ytext-dy;
text(xtext,ytext,['cc=' num2str(cc(1,2),3)]);ytext=ytext-dy;
if strcmp(vht_is,'given'),
   text(xtext,ytext,['V_{HT} given as input']);ytext=ytext-dy;
end




