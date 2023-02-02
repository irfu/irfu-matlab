function [slope,cc]=irf_walen(v,b,n,vht,tpar,tperp,tint,tint_ex)
%IRF_WALEN   Walen test and estimate its goodness
%
% [slope,cc]=irf_walen(v,b,n,vht,[tperp,tpar,tint,tint_ex]);
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
% See also IRF_VHT, IRF_VHT_PLOT
%

% Alessandro Retino`
% 2004/11/26

if nargin < 5, anys_mode = 0;
else
  if ~isempty(tperp) &&  ~isempty(tpar), anys_mode = 1;
  else, anys_mode = 0;
  end
end

%if time interval is not given define it from the size of the input vectors
if nargin < 7 || isempty(tint)
  tint=[min([v(1,1),b(1,1),n(1,1)]) max([v(end,1),b(end,1),n(end,1)])];
  tint_ex = [];
end

%display time interval
strint = irf_disp_iso_range(tint,1);

%define common time interval for input vectors
n = irf_tlim(n,tint);

%interpolate b,n,tpat,tperp to v
v = irf_resamp(v,n);
b = irf_resamp(b,n);
if anys_mode
  tperp = irf_resamp(tperp,n);
  tpar = irf_resamp(tpar,n);
end

NTHRESH=0.6;
ii = find(n(:,2)>NTHRESH);
n = n(ii,:);
v = v(ii,:);
b = b(ii,:);
if anys_mode
  tperp = tperp(ii,:);
  tpar = tpar(ii,:);
end

% exclude subintervals
if nargin > 5
  for i=1:size(tint_ex)
    disp(['excluding ' irf_disp_iso_range(tint_ex(i,:),1)])
    n = irf_tlim(n,tint_ex(i,1),tint_ex(i,2),1);
    v = irf_tlim(v,tint_ex(i,1),tint_ex(i,2),1);
    b = irf_tlim(b,tint_ex(i,1),tint_ex(i,2),1);
    if anys_mode
      tperp = irf_tlim(tperp,tint_ex(i,1),tint_ex(i,2),1);
      tpar = irf_tlim(tpar,tint_ex(i,1),tint_ex(i,2),1);
    end
  end
end

n = [n repmat(n(:,2),1,2)];

%calculate velocity in HT frame, Alfven velocity and pressure anisotropy

vtransf = n;
vtransf(:,2:4)=v(:,2:4)-repmat(vht,size(n,1),1);

valfv = b;
valfv(:,2:4) = 22*b(:,2:4)./ sqrt(n(:,2:4));

if anys_mode
  alpha = n(:,1:2);
  b = irf_abs(b);
  alpha(:,2)=17.33*n(:,2).*( tpar(:,2)-tperp(:,2) )./( b(:,5).^2);
  %alpha(alpha(:,1)<=iso2epoch('2007-03-27T05:07:58.000Z'),2) = 0;
  %irf_plot(alpha), keyboard
  alpha(:,2) = my_smooth(alpha(:,2),3);
  alpha= [alpha repmat(alpha(:,2),1,2)];
  valfv(:,2:4)=valfv(:,2:4).*sqrt( 1 - alpha(:,2:4));
end

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
  % return         % why RETURN; wyli @ IRFU 2015-11-21;
end
figure(117), clf
plot(valfv(:,2),vtransf(:,2),'b.',valfv(:,3),vtransf(:,3),'g.',valfv(:,4),vtransf(:,4),'r.');
axis equal;grid on;
ht=irf_pl_info([mfilename ' ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))]); set(ht,'interpreter','none','FontSize', 5);

title(['Walen test ' strint])
xlabel('V_{A} [km/s]');ylabel('V-V_{HT} [km/s]')
legend('x','y','z');

xx=get(gca,'XLim');
yy=get(gca,'YLim');
dx=(xx(2)-xx(1))/10;
dy=(yy(2)-yy(1))/10;
text(xx(1)+dx,yy(1)+2*dy,  [' slope: ' num2str(slope,'%1.2f')])
text(xx(1)+dx,yy(1)+1.5*dy,['    cc: ' num2str(cc,'%1.2f')])
text(xx(1)+dx,yy(1)+dy,    ['offset: ' num2str(p(2),'%1.2f')])
text(xx(1)+dx,yy(1)+.5*dy, ['   vht: [' num2str(vht) '] km/s GSE'])

figure(118), clf
irf_plot({valfv,vtransf},'comp','-')            % change '.' to '_'
legend('V_A','V-V_{HT}')
ylabel('V [km/s]');

figure(119), clf
irf_subplot(2,1,-1)
if anys_mode
  irf_plot(alpha,'.')
end
ylabel('\alpha');
irf_subplot(2,1,-2)
irf_plot(b,'.')
ylabel('B [nT]');

function sig = my_smooth(sig, niter)
% Remove spikes and smoothen the signal

%clf, plot(sig,'k.'), hold on

if niter > 0
  for i=1:niter
    mm = mean(sig);
    ii = find(abs(sig - mm) > 1.2*std(sig));
    sig_tmp = sig;
    sig_tmp(ii) = mm;
    sig_tmp = smooth(sig_tmp);
    sig(ii) = sig_tmp(ii);
    %plot(sig,irf_lstyle(i));
  end
end

sig = smooth(sig);
%plot(sig,irf_lstyle(niter+1)); hold off