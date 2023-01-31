function [newdata]=irf_lowpass(data,fcut,fhz)
%IRF_LOWPASS lowpass filter
%
% [newdata]=irf_lowpass(data,fcut,fhz)
% filter the data through low or highpas filter with max frequency fcut
% and subtract from the original
% norder=4 for elliptic filter
%
% see also ELLIP, FILTFILT
%


fnyq=fhz/2;

rp=0.1;
rs=60;
norder=4;
dedata=detrend(data);
rest=data-dedata;

[b,a]= ellip(norder,rp,rs,fcut/fnyq);
newdata=filtfilt(b,a,dedata) +rest;
