function caa_ms_pl_xoff_stats(xoff)
%CAA_MS_PL_XOFF_STATS plot statistics for DdsiX1-4
%
% caa_ms_pl_xoff_stats(DdsiXn)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

STEP = 0.1;
d = -1:STEP:5;

yy=2015:2019; % yy length c size
c = 'krgbm';
leg = {};

clf
for i=1:length(yy)
  st = toepoch([yy(i) 1 1 0 0 0]);
  et = toepoch([yy(i)+1 1 1 0 0 0]);
  x = irf_tlim(xoff,st,et);
  if isempty(x)
    irf_log('proc',sprintf('no offsets for %d',yy(i)))
    continue
  end
  x = x(:,2);

  pdf = zeros(size(d));

  for j=1:length(d)
    pdf(j) = length(find( x>=d(j) & x<d(j)+STEP ));
  end

  pdf = pdf/abs(sum(x));

  im = find( pdf==max(pdf), 1, 'last');
  ii = im-3:im+3; ii(ii<=0) = [];
  cf = fit(d(ii)'+STEP/2,(pdf(ii))','gauss1');
  cf = coeffvalues(cf);

  [xx,ii] = sort([d+STEP/2 cf(2)]);
  pdf = [pdf cf(1)]; %#ok<AGROW>
  pdf = pdf(ii);

  plot(xx,pdf,[c(i) '-o'])
  hold on
  text(cf(2),cf(1),sprintf('. %.2f mV/m',cf(2)))
  leg = [leg {sprintf('%d %.2f mV/m',yy(i),cf(2))}]; %#ok<AGROW>
end
legend(leg)
hold off
grid
ylabel('Offset PDF')
xlabel('\Delta E_x [mV/m]')