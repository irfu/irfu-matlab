function caa_sh_pl_xoff(dE, dAmp, dt_int, weight, tint)
%CAA_SH_PL_XOFF  visualize offset study results
%
% caa_sh_pl_xoff(dE, dAmp, dt_int, weight,[tint])
%
% See also CAA_SH_XOFF_BATCH

% Copyright 2007 Yuri Khotyaintsev

narginchk(4,5)

if nargin==5
	if ~all( size([1 2]) == [1 2] )
		error('TINT must be [START_EPO STOP_EPO]')
	end
	if tint(2) <= tint(1)
		error('TINT must be [START_EPO STOP_EPO], STOP_EPO>START_EPO')
	end
	ii = find( dE(:,1)>=tint(1) & dE(:,1)<tint(2) );
	dE = dE(ii,:);
	dAmp = dAmp(ii,:);
	dt_int = dt_int(ii,:);
	weight = weight(ii,:);
end

for cli=1:4, dE(abs(dE(:,cli+1)) > 5, cli+1) = NaN; end

clrs='krgb';
figure(77), clf
subplot(3,1,1)
for cli=1:4
	irf_plot(dE(:,[1 cli+1]),[clrs(cli) '.-'])
	if cli==1, hold on, end
end
hold off
%legend('C1','C2','C3','C4')
ylabel('dE [mV/m]')
xlabel('')
set(gca,'XTickLabel',[], 'YLimMode','auto')
title('Offset summary')

ddE = zeros(1,4);
subplot(3,1,2)
for cli=1:4
	ii = find( ~isnan(dE(:,cli+1)) );
	ii = find( ~isnan(dE(:,cli+1)) & (...
		abs( dE(:,cli+1) - mean( dE(ii,cli+1) ) ) < std( dE(ii,cli+1) ) ));
	ci = dt_int(ii).*weight(ii)/sum(dt_int(ii).*weight(ii));
	ddE(cli) = sum( dE(ii,cli+1).*ci );
	irf_plot([dE(:,1) dE(:,cli+1)-ddE(cli)],[clrs(cli) 'x'])
	if cli==1, hold on, end
	irf_plot([dE(ii,1) dE(ii,cli+1)-ddE(cli)],[clrs(cli) 'O'])
end
hold off
set(gca,'YLimMode','auto')
ylabel('dE-<dE> [mV/m]')

subplot(3,2,5)
for cli=1:4
	di = dist_dE(dE(:,cli+1),ddE(cli));
	plot(di(:,1),di(:,2), [clrs(cli) 'x'])
	if cli==1, hold on, end
end
hold off, grid on
set(gca,'YLim',[0 1.1],'XTick',[-2 -1 -.5 0 .5 1 2])
xlabel('dE-<dE> [mV/m]')
ylabel('distribution')

legend(sprintf('<C1>=%.2f',ddE(1)),sprintf('<C2>=%.2f',ddE(2)),...
	sprintf('<C3>=%.2f',ddE(3)),sprintf('<C4>=%.2f',ddE(4)),...
	'Location','NorthEastOutside')

function res = dist_dE(dE,ddE)
% find distribution of dE

STEP2 = .05;

res = zeros(41,2);
res(:,1) = -2:(STEP2*2):2;
dE = dE(~isnan(dE)) - ddE;

for i=1:length(dE)
	ii = find( res(:,1)-STEP2 < dE(i) & dE(i) <= res(:,1)+STEP2 );
	res(ii,2) = res(ii,2) + 1;
end
res(:,2) = res(:,2)/max(res(:,2));
