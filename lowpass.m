function [newdata]=lowpass(data,fcut,fhz)
%function [newdata]=lowpass(data,fcut,fhz)
% filter the data through low or highpas filter with max frequency fcut
% and subtract from the original
% norder=4 for elliptic filter
%
% $Id$
%
% see also ELLIP, FILTFILT

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_lowpass')


fnyq=fhz/2;

rp=0.1;
rs=60;
norder=4;
dedata=detrend(data);
rest=data-dedata;

[b,a]= ellip(norder,rp,rs,fcut/fnyq);
newdata=filtfilt(b,a,dedata) +rest;
