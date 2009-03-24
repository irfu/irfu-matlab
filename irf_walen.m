function [slope,cc]=irf_walen(v,b,n,vht,tint,tint_ex)
%IRF_WALEN   Walen test and estimate its goodness
%
% [slope,cc]=irf_walen(v,b,n,vht,tint,tint_ex);
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
% tint_ex - time intervals to exclude
%
% OUTPUT
% [slope, cc] - slope and correlation coefficient of the Walen plot 
%
% See also IRF_WALEN2, IRF_VHT
%
% $Id$

% Alessandro Retino`
% 2004/11/26


%if time interval is not given define it from the size of the input vectors
if nargin < 5, % if tint is not given
  tint=[min([v(1,1),b(1,1),n(1,1)]) max([v(end,1),b(end,1),n(end,1)])];
elseif isempty(tint), % if tint=[], define tint from v and b time axis
   tint=[min([v(1,1),b(1,1),n(1,1)]) max([v(end,1),b(end,1),n(end,1)])];  
end

%display time interval
strint=[epoch2iso(tint(1)) ' -- ' epoch2iso(tint(2)) ];

%define common time interval for input vectors
n = irf_tlim(n,tint);

%interpolate b,n,tpat,tperp to v
v = irf_resamp(v,n);
b = irf_resamp(b,n);

% exclude subintervals
if nargin > 5
	for i=1:size(tint_ex)
		disp(['excluding ' irf_disp_iso_range(tint_ex(i,:),1)])
		n = irf_tlim(n,tint_ex(i,1),tint_ex(i,2),1);
		v = irf_tlim(v,tint_ex(i,1),tint_ex(i,2),1);
		b = irf_tlim(b,tint_ex(i,1),tint_ex(i,2),1);
	end
end

n = [n repmat(n(:,2),1,2)];

%tpar=irf_interp(tpar,n);
%tperp=irf_interp(tperp,n);

%calculate velocity in HT frame, Alfven velocity and pressure anisotropy
vtransf(:,1)=n(:,1);
vtransf(:,2:4)=v(:,2:4)-repmat(vht,size(n,1),1);

%alpha(:,1)=n(:,1);
%b=irf_abs(b);
%alpha(:,2)=17.33*n(:,2).*( tpar(:,2)-tperp(:,2) )./( b(:,5).^2);
%alpha(:,2)=17.33*n(:,2).*( tpar(:,2)-tperp(:,2) )./( b(:,2).^2+b(:,3).^2+b(:,4).^2);
%alpha(:,2)=0.0;
%alpha= [alpha repmat(alpha(:,2),1,2)];
%irf_plot(alpha)


valfv(:,1)=n(:,1);
%valfv1(:,1)=n(:,1);

valfv(:,2:4)=22*b(:,2:4)./ sqrt(n(:,2:4));
%valfv1(:,2:4)=22*b(:,2:4)./ sqrt( n(:,2:4));

%valfv(:,2:4)=valfv(:,2:4).*sqrt( 1 - alpha(:,2:4))


%diff(:,1)=n(:,1);
%diff(:,2:4)=valfv1(:,2:4)-valfv(:,2:4);
%irf_plot({vtransf,valfv})

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
cc=corr(1,2);
p=polyfit(valfvtot,vtransftot,1);



slope=p(1);
if nargout>0
    disp(strint);
    disp(['Offset: ' num2str(p(2))])
    return
end

plot(valfv(:,2),vtransf(:,2),'b.',valfv(:,3),vtransf(:,3),'g.',valfv(:,4),vtransf(:,4),'r.');
axis equal;grid on;
ht=irf_pl_info([mfilename ' ' datestr(now)]); set(ht,'interpreter','none','FontSize', 5);

title(['Walen test ' strint])
xlabel('V_{A} [km/s] GSE');ylabel('V-V_{HT} [km/s] GSE')
legend('x','y','z');

xx=get(gca,'XLim');
yy=get(gca,'YLim');
dx=(xx(2)-xx(1))/10;
dy=(yy(2)-yy(1))/10;
text(xx(1)+dx,yy(1)+2*dy,  [' slope: ' num2str(slope,'%1.2f')])
text(xx(1)+dx,yy(1)+1.5*dy,['    cc: ' num2str(cc,'%1.2f')])
text(xx(1)+dx,yy(1)+dy,    ['offset: ' num2str(p(2),'%1.2f')])
text(xx(1)+dx,yy(1)+.5*dy, ['   vht: [' num2str(vht) '] km/s GSE'])