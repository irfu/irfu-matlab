function [newdata]=lowpass(data,fcut,fhz)
%function [newdata]=lowpass(data,fcut,fhz)
% filter the data through low or highpas filter with max frequency fcut
% and subtract from the original
% norder=4 for elliptic filter
%
% $Id$
%
% see also ELLIP, FILTFILT


fnyq=fhz/2;

rp=0.1;
rs=60;
norder=4;
dedata=detrend(data);
rest=data-dedata;

[b,a]= ellip(norder,rp,rs,fcut/fnyq);
newdata=filtfilt(b,a,dedata) +rest;
