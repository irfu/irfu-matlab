function [iV, dV] = c_4_v_xcorr(tint,B1,B2,B3,B4,R1,R2,R3,R4)
%C_4_V_XCORR  Automatically estimate boundary velocity
%
%  [V, dV] = c_4_v_xcorr(tint,B1,B2,B3,B4,R1,R2,R3,R4)
%
%  Automatically estimate velocity of a boundary. The time difference
%  between two of the spacecrfat if found from minimum of the sum
%  S = SUM ( [data1 - data2]^2 ). The error is defined as DT at which S is
%  double the minimum value.
%
%  Output:
%     V - velocity in GSE
%    dV - error on velocity


% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin == 0, help c_4_v_xcorr; return; end

if isa(tint,'GenericTimeArray')
  tint = [tint.start.utc '/' tint.stop.utc];
end
if ischar(tint) % assume time interval input in character format
  tint=irf_time(tint,'utc>tint');
end

if isa(B1,'TSeries'), B1 = [B1.time.epochUnix double(B1.data)]; end %#ok<NASGU>
if isa(B2,'TSeries'), B2 = [B2.time.epochUnix double(B2.data)]; end %#ok<NASGU>
if isa(B3,'TSeries'), B3 = [B3.time.epochUnix double(B3.data)]; end %#ok<NASGU>
if isa(B4,'TSeries'), B4 = [B4.time.epochUnix double(B4.data)]; end %#ok<NASGU>

if isa(R1,'TSeries'), R1 = [R1.time.epochUnix double(R1.data)]; end %#ok<NASGU>
if isa(R2,'TSeries'), R2 = [R2.time.epochUnix double(R2.data)]; end %#ok<NASGU>
if isa(R3,'TSeries'), R3 = [R3.time.epochUnix double(R3.data)]; end %#ok<NASGU>
if isa(R4,'TSeries'), R4 = [R4.time.epochUnix double(R4.data)]; end %#ok<NASGU>

b1 = [];
c_eval('b?=irf_tlim(B?(:,1:4),tint);')
c_eval('r?= interp1(R?(:,1),R?(:,[2 3 4]),tint(1),''spline'',''extrap'');')

%% Initial velocity estimate
fs = 1/(b1(2,1)-b1(1,1)); % Sampling frequency

R = 0; dR12 = zeros(6,3); pos = 0; t0 = zeros(6,1); dt = 0;
for cli1 = 1:3
  for cli2 = (cli1+1):4
    pos = pos + 1;
    %sprintf('(%d) - > (%d)',cli1,cli2)
    c_eval('data1=b?(:,2:4); data2=b!(:,2:4); dt=b?(1,1)-b!(1,1);',cli1,cli2);
    
    t0(pos) = mycorr( data1, data2, fs );
    t0(pos) = t0(pos) + dt; % Correct for time shift between the timelines
    
    eval(irf_ssub('dR12(pos,:) = r? - r!;',cli1,cli2))
    R = R + dR12(pos,:)'*dR12(pos,:);
  end
end

R = R/6;
S = sum((t0*ones(1,3)).*dR12,1) / 6;
M = S/R;
V = M/sum(M.^2);

vn = irf_norm(V);
disp([' V=' num2str(irf_abs(V,1),3) '*[' num2str(vn(end-2:end),' %5.2f') '] km/s GSE [initial]'])

%% Improve the velocity and estimate errors
[R,~,dR1,dR2,dR3,dR4]=c_4_r(r1,r2,r3,r4); % R Volumetric tensor

dt = zeros(4,1);
c_eval('dt(?) = sum(dR?.*V)/sum(V.*V);') % Time shifts using the initial V
tint1 = tint + [min(dt) max(dt)];
it = irf_tlim(b1(:,1),tint1);
iB1 = []; iB2 = []; iB3 = []; iB4 = [];
c_eval('iB? = irf_resamp([b?(:,1)-dt(?) b?(:,2:end)],it);')
ii1 = irf_find_comm_idx(b1,it);
B0 = [b1(:,1) ones(size(b1,1),3)*NaN];
B0(ii1,:) = [it (iB1(:,2:4) + iB2(:,2:4) + iB3(:,2:4) + iB4(:,2:4))/4]; %#ok<NASGU> % Average B profile

it0 = zeros(4,1); ddt = zeros(4,1);
for cli = 1:4
  btmp = [];
  c_eval('btmp=b?;',cli);
  [ it0(cli), ddt(cli) ] = mycorr(btmp,B0,fs);
  dt = btmp(1,1) - B0(1,1);
  it0(cli) = it0(cli) + dt;
end

M = ( it0(1)*dR1 + it0(2)*dR2 + it0(3)*dR3 + it0(4)*dR4)/4/R;
dM = sqrt( (ddt(1)*dR1/R).^2 );
c_eval('dV?=0; Mplus=M+0.25*ddt(?)*dR?/R;Mminus=M-0.25*ddt(?)*dR?/R;dV?=(Mplus/sum(Mplus.^2)-Mminus/sum(Mminus.^2))/2;');
dV=sqrt(dV1.^2+dV2.^2+dV3.^2+dV4.^2);

%dM = sqrt( (ddt(1)*dR1/R).^2 + (ddt(2)*dR2/R).^2 + (ddt(3)*dR3/R).^2 + (ddt(4)*dR4/R).^2);

iV = M/sum(M.^2);
%dV = dM/sum((M+dM).^2);
%dV = dM/sum(M.^2);

ivn = irf_norm(iV);
disp([' V=' num2str(irf_abs(iV,1),3) '*[' num2str(ivn(end-2:end),' %5.2f') '] km/s GSE [final]'])
disp([' dV= +/- [ ' num2str(dV,' %5.2f') '] km/s GSE'])

end

%% Help functions
function [ t0, dt, tplus, tminus, corr ] = mycorr( data1, data2, fs )
%MYCORR  Cross correlation "by eye"
%
%   Minimize S = SUM ( [data1 - data2]^2 )
%
%   [ t0, dt, tplus, tminus, corr ] = mycorr( data1, data2, fs )
%
%   If x and y are not the same length, the shorter vector is zero-padded
%   to the length of the longer vector
%
%   t0 - time delay between data1 and data2
%   dt = (tplus + tminus)/2

if nargin <3, fs = 1; end

ndata = max(size(data1,1),size(data2,1));

if ndata < 5, error('Must have at least 5 points!'), end

ncomp = size(data1,2);

% Remove mean - XXX: seems to increase the error
%data1 = data1 - ones(size(data1,1),1)*mean(data1,1);
%data2 = data2 - ones(size(data2,1),1)*mean(data2,1);

if length(data1) < ndata
  data1(end:ndata,:) = NaN;
elseif length(data2) < ndata
  data2(end:ndata,:) = NaN;
end

corr = zeros(ndata*2-1,1);

data2 = [NaN*ones(ndata-1,ncomp); data2; NaN*ones(ndata-1,ncomp)];


for i=(-ndata+1):ndata-1
  S = ( data1 - data2(i+ndata:i+2*ndata-1,:) ).^2;
  if any(~isnan(S))
    corr(i+ndata) = sum(S(~isnan(S)))/(ndata - abs(i))/ncomp;
  else, corr(i+ndata) = NaN;
  end
end

OFF = 1;
ii = find(corr == min(corr(1+OFF:end-OFF))); % take out the last points for security

warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
p = polyfit((ii-OFF:ii+OFF)',corr(ii-OFF:ii+OFF),2);
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
t0 = -p(2)/p(1)/2;
if abs(t0-ii)<1
  corr_min = polyval(p,t0);
else
  disp('Discarding parabolic fit!')
  t0 = ii;
  corr_min = corr(ii(1));
end

% We define the error as DT at which the correlation doubles
corr_plus = corr - 2*corr_min;

ii_plus = find(corr_plus > 0);
if isempty(ii_plus), tplus = NaN; tminus = NaN;
else
  ii_plus = ii_plus - ii;
  ii1 = find(ii_plus>0);
  ii1 = ii_plus(ii1(1)) + ii;
  tplus = find_zero(corr_plus,ii1,ii1-1);
  ii1 = find(ii_plus<0);
  ii1 = ii_plus(ii1(end)) + ii;
  tminus = find_zero(corr_plus,ii1,ii1+1);
end

dt = ( tplus - tminus) /2;
tplus = tplus - ndata;
tminus = tminus - ndata;
t0 = t0 - ndata;

% Diagnostics
if false
  figure; %#ok<UNRCH>
  dt_tmp = (1:(ndata*2-1)) - ndata;
  plot(dt_tmp/fs,corr,'.-')
  set(gca, 'YLim',[0 3*corr_min])
  hold on
  plot(t0/fs,corr_min,'r*')
  text(t0/fs,corr_min,sprintf(' t0 = %.3f s +/- %.3f s (%.2f +/- %.2f)',t0/fs,dt/fs,t0,dt))
  plot(tplus/fs,corr_min*2,'b*')
  text(tplus/fs,corr_min*2,sprintf(' t+ = %.2f s (%.2f)',tplus/fs,tplus))
  plot(tminus/fs,corr_min*2,'b*')
  text(tminus/fs,corr_min*2,sprintf(' t- = %.2f s (%.2f)',tminus/fs,tminus))
  hold off
  grid
end

% Prepare the output
tplus = - tplus/fs;
tminus = - tminus/fs;
t0 = - t0/fs;
dt = dt/fs;

end

function res = find_zero(y,i1,i2)

a = ( y(i1) - y(i2) ) / ( i1 - i2 );
b = ( y(i1) + y(i2) - a*( i1 + i2 ) ) / 2;

res = -b/a;
end
