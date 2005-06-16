function h=c_staff_saspec(time,dat)
%function h=c_staff_saspec(dat,time)
%plots STAFF SA spectrgram
%	dat - data matrix with frequencies
%	time - time axis
%Examples:
%	[timeB,datB]=c_staff_getsa([2002 03 03 22 23 00],60,1);
%	c_staff_saspec(timeB,datB);
%
%See also: c_staff_getsa
%David Sundkvist, IRFU, 2005
%Last change: 2005-06-16


if nargin==0 | nargin==1
	help c_staff_saspec;return;
end

if nargin>2
	error('Too many arguments');help c_staff_saspec;return;
end


%fix data since staff data comes in k x l x m matrices..
[k,l,m]=size(dat);
if l==3
	%we have B data
	dataBx=dat(:,1,:);
	dataBy=dat(:,2,:);
	dataBz=dat(:,3,:);
	dat=squeeze(dataBx)+squeeze(dataBy)+squeeze(dataBz);
elseif l==2
	%we have E data
	dataEx=dat(:,1,:);
	dataEy=dat(:,2,:);
	dat=squeeze(dataEx)+squeeze(dataEy);
else 
	error('unknown format of data matrix')
end
	


	
%Construct frequency axis, the values are from a iscmd query
freq=[8.77    11.05   13.92   17.54   22.1    27.84   35.08   44.19   55.68   70.15   88.39   111.36  140.31  176.78  222.72  280.62  353.55  445.45  561.23  707.11  890.9   1122.46    1414.21 1781.8  2244.92 2828.43 3563.59];
freq=freq/1E3;

%plot spectrogram
t0=time(1,1);h=pcolor(time-t0,freq,log10(dat));

%set(gca,'yscale','log') %messes up the y-axis values even if spaceing is fine..


c_bar=colorbar;shading flat;add_timeaxis(gca,t0);



ylabel('Frequency [kHz]');
set(gca,'tickdir','out')

hh=get(c_bar,'ylabel');

if l==3
	set(hh,'string',['(pT)^2/Hz']);
elseif l==2
	set(hh,'string',['(mV/m)^2/Hz']);
else
	set(hh,'string',['(X)^2/Hz']);
end




