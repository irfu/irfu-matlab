function [spin_period,err_angle,err_angle_mean, phc_coef] = c_spin_period(phase, flag)
%C_SPIN_PERIOD  Compute cluster spin period
%
% [SPIN_PERIOD,ERR_ANGLE,ERR_ANGLE_MEAN, POL_COEF] = C_SPIN_PERIOD(PHASE_2, flag)
%
% Compute Cluster spin period from PHASE_2 or PHASE (FLAG=1)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
if nargin<2, flag = 0; end
spin_period = [];err_angle  = []; err_angle_mean = [];

ph = phase; tref = phase(1,1); ph(:,1) = ph(:,1) - tref;
if flag
  phc = unwrap(ph(:,2)/180*pi); phc_coef = polyfit(ph(:,1),phc,1);
else
  phc = unwrap(ph(1:2,2)/180*pi); phc_coef = polyfit(ph(:,1),phc,1);
end
if isnan(phc_coef(1))
  irf_log('proc','Cannot determine spin period!');
  return
end
for j=1:floor(log10(length(ph(:,1))))
  ii = 10^j;
  dp = ph(ii,2) - mod(polyval(phc_coef,ph(ii,1))*180/pi,360);
  dpm = [dp dp-360 dp+360];
  dph = dpm(abs(dpm)<180);
  phc_coef(1) = phc_coef(1) + dph*pi/180/ph(ii,1);
end
diffangle = mod(ph(:,2) - polyval(phc_coef,ph(:,1))*180/pi,360);
diffangle = abs(diffangle); diffangle = min([diffangle';360-diffangle']);
err_angle_mean = mean(diffangle); err_angle = std(diffangle);
if err_angle>1 || err_angle_mean>1
  irf_log('proc',['Spin period is changing! Phase errors>1deg. err=' num2str(err_angle) 'deg.']);
end
spin_period = 2*pi/phc_coef(1);
