function [time,varargout]=c_staff_getsa(time,dt,SC,EBflag,plotflag)

%First simple file to get staff SA data and plot the spectrograms. 
%Usage: [time,varargout]=c_staff_getsa(starttime,dt,[SC],[EBflag],[plotflag])
%
%where
%       starttime - start of interval on format [yyyy mm dd hh mm ss]
%       duration - duration in seconds
%	SC - spacecraft (can only be one so far (def. 1)
%	EBflag - 1 for B, 2 for E, 3 for both (see examples) (def. 1)
%	plotflag - 0 = no plot, 1 = plot (def. 0)
% 
%Examples:
%	[timeB,datB]=c_staff_getsa([2002 03 03 22 23 00],60,1);
%	[timeB,datB]=c_staff_getsa([2002 03 03 22 23 00],60,1,1,1);
%	[timeE,datE]=c_staff_getsa([2002 03 03 22 23 00],60,1,2);
%	[timeB,datB,timeE,datE]=c_staff_getsa([2002 03 03 22 23 00],60,1,3);
%	[timeE,datE]=c_staff_getsa([2002 03 03 22 23 00],60,1,2,1);
%
%
%It only downloads the 8-4096Hz channel so far.. This is written in 30 min!
%
%See also: c_staff_saspec
%
%David Sundkvist, IRFU, 2005.
%Last change: 2005-06-16.

if nargin<5
	plotflag=0;	%default no plot
end
if nargin<4
	EBflag=1;	%default load B
end
if nargin<3
	SC=1;		%so far only on spacecraft is downloaded
end
if nargin<1
	help c_staffgetsa;
	return;
end

DB=Mat_DbOpen('disco:10');


%time=[2002 03 03 22 23 00];
%dt=60;

if EBflag==1 | EBflag==3
	[timeB,dataB]=isGetDataLite(DB,time,dt,'Cluster',num2str(SC),'staff','B_SA','a_B','8-4096Hz','');
end
if EBflag==2 | EBflag==3
[timeE,dataE]=isGetDataLite(DB,time,dt,'Cluster',num2str(SC),'staff','E_SA','a_E','8-4096Hz','');
end




%-------Output---------------------------
if EBflag==1
	time=timeB;
	varargout(1)={dataB};
elseif EBflag==2
	time=timeE;
	varargout(1)={dataE};
elseif EBflag==3
	time=timeB;
	varargout(1)={dataB};
	varargout(2)={timeE};
	varargout(3)={dataE};
end	


%-------Plot---------------------------

if plotflag
	if EBflag==1
		c_staffsaspec(time,dataB);
	elseif EBflag==2
		c_staffsaspec(time,dataE);
	else
		error('no support for multiplot yet')
	end
end


%-------End-------------------------------------------
