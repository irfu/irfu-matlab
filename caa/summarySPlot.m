function h=summarySPlot(sp,flag)
%summarySPlot make a simple summary plot with B (full res), P and E/Es DSI.
%
% summarySPlot(sp)
% 
% Input:
% sp -storage path (data directory)
%
% $Revision$  $Date$
%
% See also ClusterProc/summaryPlot, summaryRPlot

% Copyright 2004 Yuri Khotyaintsev

if nargin<2, flag=0; end

load mP
load mEDSI diEs1p34 diE2p1234 diEs3p34 diE4p1234
load mB B1 B2 B3 B4

tt = diEs1p34(:,1)*ones(1,8);

h = c_pl_tx(tt,tt,tt,tt,2:8);
st=min(B1(1,1),B2(1,1));
se=max(B1(end,1),B2(end,1));

for j=1:4
	axes(h(j))
	c_pl_tx(av_abs(B1),av_abs(B2),av_abs(B3),av_abs(B4),1+j);
	ylabel(['B' n2c(j) ' [nT]'])
	set(gca,'XLim',[st se])
	set(gca,'XTickLabel',[]),xlabel('')
	set(gca,'YLim',.99*get(gca,'YLim'))
end

axes(h(5))
c_pl_tx(NVps1,NVps2,NVps3,NVps4,2);
set(gca,'XLim',[st se])
set(gca,'XTickLabel',[]),xlabel('')
set(gca,'YScale','log')
set(gca,'YLim',.99*get(gca,'YLim'))
ylabel('n [cm^{-3}]')

for j=1:2
	axes(h(5+j))
	c_pl_tx(diEs1p34,diE2p1234,diEs3p34,diE4p1234,1+j);
	ylabel(['E' n2c(j) ' [mV/m]'])
	set(gca,'XLim',[st se])
	if j==1, set(gca,'XTickLabel',[]),xlabel(''), end
	set(gca,'YLim',.99*get(gca,'YLim'))
end

addPlotInfo
for j=1:7, axes(h(j)), end

if nargout < 1, clear h, end

function res=n2c(n)
switch n
case 1,res='x';
case 2,res='y';
case 3,res='z';
case 4,res='t';
otherwise
	res='undef';
end

