function vph = estimate_phase_speed(varargin) 
%
% Simple function to estimate the phase speed from the frequency wave
% number power spectrum. Fits f = v k/2 pi to the power spectrum.
% N.B. Draft version but seems to work well. Does not yet handle multiple
% modes in the same power spectrum.
% See also mms.fk_powerspectrum.m and mms.probe_align_times.m.
% Written by D. B. Graham
%
% vph = mms.estimate_phase_speed(fkpower,freq,wavenumber,[fmin])
% Inputs: (Uses output from mms.fk_powerspectrum.m)
%   power - 2D array of powers
%   freq - 1D array of frequencies
%   wavenumber - 1D array of wavenumbers
%   fmin - Set low frequency threshold of points used to estimate the
%   speed. Optional parameter; default is 100 Hz.
%
% Output:
%   vph - estimated phase speed by fitting linear dispersion relation to
%         data.

if (numel(varargin) < 3)
    help mms.estimate_phase_speed;
    return;
end

fkpower = varargin{1};
freq = varargin{2};
wavenumber = varargin{3};

if (numel(varargin)==3)
    fmin = 100;
else 
    fmin = varargin{4};
end

% Remove spurious points; specifically at large k. 
kmax = 2.0*max(wavenumber)/3.0;
powertemp = fkpower;
rmk = find(abs(wavenumber) > kmax);
rmf = find(freq < fmin);
powertemp(:,rmk) = 0.0;
powertemp(rmf,:) = 0.0;
powertemp(isnan(powertemp)) = 0.0;

% Find position of maximum power to guess vph
[wavenumbers,freqs] = meshgrid(wavenumber,freq);
maxpos = find(powertemp == max(max(powertemp)));
maxf = freqs(maxpos);
maxk = wavenumbers(maxpos);

vphguess = maxf/maxk;
if (vphguess > 0.0)
    powertemp(:,wavenumber < 0.0) = 0;
else
    powertemp(:,wavenumber > 0.0) = 0;
end

vphrange = [vphguess/3 vphguess*3];

% Find all points around this guess vph
highppos = find(powertemp > 0.3*max(max(powertemp)));
ppowers = powertemp(highppos);
fpowers = freqs(highppos);
kpowers = wavenumbers(highppos);

ppowers2 = [];
fpowers2 = [];
kpowers2 = [];
elnum = 1;

for ii = 1:length(ppowers)
    if (abs(fpowers(ii)/kpowers(ii)) > abs(vphrange(1)) && abs(fpowers(ii)/kpowers(ii)) < abs(vphrange(2)))
        ppowers2(elnum) = ppowers(ii);
        fpowers2(elnum) = fpowers(ii);
        kpowers2(elnum) = kpowers(ii);
        elnum = elnum+1;
    end
end

weights = 1+log10(ppowers2/max(ppowers2));

fun = fittype('a*x');
[fit1,gof,fitinfo] = fit(kpowers2',fpowers2',fun,'Weight', weights,'StartPoint',vphguess);
vph = fit1.a*2*pi;

end
