function [out] = irf_filt(inp,fmin,fmax,Fs,order)
%IRF_FILT   Filter time series
%
% [out] = irf_filt(inp,fmin,fmax,[Fs],[order])
%
% inp, out   - column vectors
%              if inp has more than 1 column, 
%              assume that the first column is time
%              calculate frequency from the first time step
%              assume that all time steps are the same length
%              if inp is a TSeries, out is also a TSeries
% fmin,fmax  - filter frequencies
%              if fmin = 0 do lowpass filter
%              if fmax = 0 do highpass filter
% Fs         - sampling frequency if given as [] then Fs is 
%              calculated from time series 
% order      - the order of filter (elliptical IIR type filter is used) 
%              choose uneven order, 3 or 5 is OK. 
% 
% It is not better to have higher order filter. 
% With orders above 10 for IIR filters like
% cheby1, ellip the result becomes wrong. For high pass filtering with very low 
% passband frequency (<.05) one can use [B2,A2] = cheby1(4,.3,.1,'high'); 
%
% Example: 
%    def=irf_filt(de,0,.1,25,3); 
%    lowpass filter E at .1Hz
%
% Revised to accept TS format. 

isaTSeries = isa(inp,'TSeries');
if isaTSeries
    inptemp = inp;
end

if ((nargin < 4) || (isempty(Fs))) 
    if isaTSeries
        Fs = 1/(inp.time(2)-inp.time(1));
    else
        Fs=1/(inp(2,1)-inp(1,1));
    end
    irf_log('proc',['Using sampling frequency ',num2str(Fs),' Hz']);
end % estimate sampling frequency
if nargin > 4
    %irf_log('proc',['You have specified '  num2str(order) '-th filter order (use uneven order)']);
    n=order; % use this order for filters
end
fmin=fmin/(Fs/2);
fmax=fmax/(Fs/2);if (fmax > 1);fmax=1;end

if isaTSeries
    out=double(inp.data);
else
    out=inp;
end

Rp=.5;Rs=60;fact=1.1; % fact defines the width between stopband and passband
if fmin==0
  if fmax == 1, return;end
    if nargin < 5
        n=ellipord(fmax,min(fmax*fact,0.9999),Rp,Rs);
    end
	irf_log('proc',['using ' num2str(n) '-th order ellip lowpass filter']);
	[B,A] = ellip(n,Rp,Rs,fmax);
elseif fmax ==0
    if nargin < 5
        n=ellipord(fmin,min(fmin*1.1,0.9999),Rp,Rs);
    end 
	[B,A] = ellip(n,Rp,Rs,fmin,'high');
	%irf_log('proc',['using ' num2str(n) '-th highpass order filter']);
else
	%[n wn]=ellipord(fmax,fmax*1.1,Rp,Rs);
	%sprintf('using %d-th order ellip irf_lowpass filter',n)
	%[B1,A1] = ellip(n,Rp,Rs,fmax);
	if nargin < 5
    	n=ellipord(fmax,min(fmax*1.3,0.9999),Rp,Rs);
	end
	irf_log('proc',['using ' num2str(n) '-th order ellip lowpass filter']);
	[B1,A1] = ellip(n,Rp,Rs,fmax);
	if nargin < 5
    	n=ellipord(fmin,fmin*.75,Rp,Rs);
	end
	irf_log('proc',['using ' num2str(n) '-th order ellip high pass filter']);
	[B2,A2] = ellip(n,Rp,Rs,fmin,'high');
end

if isaTSeries
    inp = double(inp.data);
end

% find NaN and put to zero (in output set back to NaN
ind_NaN=find(isnan(inp)); 
inp(ind_NaN)=0;

nColumnsToFilter=size(inp,2);
    nColumnsToFilter_3=size(inp,3);             % for pressure or temperature
    nColumnsToFilter = nColumnsToFilter * nColumnsToFilter_3;
iStartColumn=1; % from which column start filtering
if nColumnsToFilter>1 % assume that first column is time
    iStartColumn=2;
end

if isaTSeries
    iStartColumn=1;
end

if ((fmin ~= 0) && (fmax ~= 0))
	for iCol=iStartColumn:nColumnsToFilter
	out(:,iCol) = filtfilt(B1,A1,inp(:,iCol)); 
	out(:,iCol) = filtfilt(B2,A2,out(:,iCol)); 
	end
else
	for iCol=iStartColumn:nColumnsToFilter
	out(:,iCol) = filtfilt(B,A,inp(:,iCol)); 
	end
end
out(ind_NaN)=NaN;

if isaTSeries
    inptemp.data = out;
    out = inptemp;
end
    

