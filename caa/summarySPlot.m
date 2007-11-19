function h=summarySPlot(sp,flag)
%summarySPlot make a simple summary plot with B (full res), P and E/Es DSI.
%
% summarySPlot(sp)
% 
% Input:
% sp -storage path (data directory)
%
% $Id$
%
% See also ClusterProc/summaryPlot, summaryRPlot

% Copyright 2004 Yuri Khotyaintsev

if nargin<2, flag=0; end

old_pwd = pwd;
cd(sp)

load mP
c_load('diEs1p34');
c_load('diE2p1234');
c_load('diEs3p34');
c_load('diE4p1234');
c_load('B?');


tt = diEs1p34(:,1)*ones(1,8);

h = c_pl_tx(tt,tt,tt,tt,2:8);
st=min(B1(1,1),B2(1,1));
se=max(B1(end,1),B2(end,1));

for j=1:4
	axes(h(j))
	c_pl_tx(irf_abs(B1),irf_abs(B2),irf_abs(B3),irf_abs(B4),1+j);
	ylabel(['B' n2c(j) ' [nT]'])
	set(gca,'XLim',[st se])
	set(gca,'XTickLabel',[]),xlabel('')
	set(gca,'YLim',.99*get(gca,'YLim'))
end

axes(h(5))
c_pl_tx(NVps1,NVps2,NVps3,NVps4,2);
set(gca,'XLim',[st se])
set(gca,'XTickLabel',[]),xlabel('')
%set(gca,'YScale','log')
set(gca,'YLim',.99*get(gca,'YLim'))
%set(gca,'YTick',[.01 .1 1 10])
ylabel('n [cm^{-3}]')

for j=1:2
	axes(h(5+j))
	c_pl_tx(diEs1p34,diE2p1234,diEs3p34,diE4p1234,1+j);
	ylabel(['E' n2c(j) ' [mV/m]'])
	set(gca,'XLim',[st se])
	if j==1, set(gca,'XTickLabel',[]),xlabel(''), end
	set(gca,'YLim',.99*get(gca,'YLim'))
end

irf_pl_add_info
for j=1:7, axes(h(j)), end

if nargout < 1, clear h, end

	
cd(old_pwd)

function res=n2c(n)
switch n
case 1,res='x';
case 2,res='y';
case 3,res='z';
case 4,res='t';
otherwise
	res='undef';
end

