function [timmin,freq,fftmat]= c_fftspec(data,nfft,overlap)
%C_FFTSPEC computes fft periodogram
%
% [timmin,freq,fftmat]= c_fftspec(data,nfft,overlap)
%
% data - AV cluster format
% nfft - size of nfft (power of 2)
% overlap - in percent
%
% $Id$

ndata = length(data);
tsec = data(:,1);

sf = 1/(tsec(2)-tsec(1));   % sampling frequency Hz
if sf<1.3*25 & sf>.7*25, sf = 25;
elseif sf<1.3*450 & sf>.7*450, sf = 450;
else, disp('cannot guess sampling frequency')
end

np= fix(nfft*(1-overlap/100));   % length of time step
n_of_timestep= fix((ndata-fix(nfft*overlap/100))/np);
nfr=nfft/2-1;
iend=n_of_timestep*np;
dt=nfft/2/sf;

timmin=linspace(tsec(1)+dt,tsec(1)+dt+iend/sf,n_of_timestep);
fftmat=zeros(nfr,n_of_timestep);

ii=1:nfft;

% loop over time steps
for i= 1: n_of_timestep
	dd=data(ii +(i-1)*np,2);
	[power,freq]=powerfft(dd,nfft,sf,1,0);
	lp=length(power);
	fftmat(:,i)=power(2:lp);
end

freq=freq(2:lp);
