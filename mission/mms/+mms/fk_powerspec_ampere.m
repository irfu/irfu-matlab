function [xvariable,yvariable,powerxy] = fk_powerspec_ampere(varargin)
%
% [fkpower,freq,wavenumber] = mms.fk_powerspec_ampere(J,B,Bback,Tints,options)
%
% Function to calculate the frequency-wave number power spectrum using
% Amperes law on a single spacecraft (based on method from Bellan, JGR, 2016). 
% Written by D. B. Graham.
%
% Input: (All data must be in TSeries format)
%       J -         Fluctuating Current density (in nA m^{-2})
%       B -         Fluctuating Magnetic field (in nT)
%       Bback -     Background magnetic field (in nT)
%       Tints -     Time interval over which the power spectrum is calculated. 
%                   To avoid boundary effects use a longer time interval
%                   for J and B. 
%
% Options:
%       numk -      Set number of wave numbers used in spectrogram. 
%       linear -    Linearly spaced frequencies. Set number to df (default is logarithmic spacing).
%       numf -      Set number of frequencies used in spectrogram.
%       wwidth -    Multiplier for Morlet wavelet width. Default is 1.
%       linear -    Use linear spacing between frequencies. df argument is
%                   required.
%       plotEpower - pass E or another quantity (TSeries) to plot power in spectrograms instead of B
%       frange -    select frequency range for k-k plots. [minf maxf]
%
% Output:
%       powerxy    - array of powers as a function of frequency and
%                    wavenumber. Power is normalized to the maximum value.
%       yvariable       - array of frequencies f in Hz.
%       xvariable - array of wavenumbers k in m^{-1}. Positive values are
%                    aligned with B and negative values are anti-aligned with B.
%
%
%
% Example: 
%   [xvecs,yvecs,Power] = mms.fk_powerspec_ampere(J,B,Bback,Tint,'linear',0.2,'numk',300,'wwidth',2,'frange',[5 10],'plotEpower',E);
%
% Example to plot:
%   pcolor(xvecs.kmag,yvecs.fkmag,log10(Power.Powerkmagf)); 
%   shading('flat');
%   xlabel('|k| (m^{-1})');
%   ylabel('f (Hz)');
%


if (nargin < 4)
	help mms.fk_powerspec_ampere;
	powerxy = NaN;
	xvariable = NaN;
	yvariable = NaN;
	return;
end

J = varargin{1};
B = varargin{2};
Bback = varargin{3};
Tint = varargin{4};

time = J.time;

B = B.resample(J);

numk = 500;
numf = 200;
uselinear = 0;
wwidth = 1;
plotEpower = 0;
frange = 0;

args=varargin(5:end);
if numel(args)>0
	haveoptions=1;
  irf.log('notice','Options were passed.');
else
	haveoptions=0;
end

while haveoptions
	l = 2;
	switch(lower(args{1}))
    case 'numk'
      if numel(args)>1 && isnumeric(args{2})
        numk = floor(args{2});
      end
    case 'numf'
      if numel(args)>1 && isnumeric(args{2})
        numf = floor(args{2});
      end
    case 'linear'
      if numel(args)>1 && isnumeric(args{2})
        df = args{2}; 
        uselinear = 1;
        irf.log('notice','Using linearly spaced frequencies');
      end
    case 'wwidth'
      if numel(args)>1 && isnumeric(args{2})
        wwidth = args{2};
      end
    case 'plotepower'
      if numel(args)>1
        E = args{2};
        E = E.resample(J);
        plotEpower = 1;
      end
    case 'frange' 
       if numel(args)>1 && isnumeric(args{2})
        frange = args{2};      
      end
    otherwise
      irf.log('warning',['Unknown flag: ' args{1}]);
      l=1;
      break;
  end
    args = args(l+1:end);
	if isempty(args), haveoptions=0; end
end

idx = tlim(B.time,Tint);

% If odd, remove last data point (as is done in irf_wavelet)
if mod(length(idx),2)
    idx(end)=[];
end

if ~uselinear
	WB = irf_wavelet(B,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*wwidth);
  WJ = irf_wavelet(J,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*wwidth);
  if plotEpower
    WE = irf_wavelet(E,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*wwidth);
  end   
else
  WB = irf_wavelet(B,'returnpower',0,'cutedge',0,'linear',df,'wavelet_width',5.36*wwidth);
  WJ = irf_wavelet(J,'returnpower',0,'cutedge',0,'linear',df,'wavelet_width',5.36*wwidth);
  if plotEpower
    WE = irf_wavelet(E,'returnpower',0,'cutedge',0,'linear',df,'wavelet_width',5.36*wwidth);
  end 
end
numf = length(WB.f);

L = length(idx);
times = time(idx);

WBx = WB.p{1,1}(idx,:);
WBy = WB.p{1,2}(idx,:);
WBz = WB.p{1,3}(idx,:);
WJx = WJ.p{1,1}(idx,:);
WJy = WJ.p{1,2}(idx,:);
WJz = WJ.p{1,3}(idx,:);

fkPower = WBx.*conj(WBx)+WBy.*conj(WBy)+WBz.*conj(WBz);
if plotEpower
  fkPowerE = WE.p{1,1}(idx,:).*conj(WE.p{1,1}(idx,:))+WE.p{1,2}(idx,:).*conj(WE.p{1,2}(idx,:))+WE.p{1,3}(idx,:).*conj(WE.p{1,3}(idx,:));
end
 
Bback = Bback.resample(times); 

Units = irf_units;
mu0 = Units.mu0;

% Negative sign needed for consistency with four spacecraft timing.
% Compared with FFT method, there a no negative frequencies used here. 
kx = -1i*mu0*(WJy.*conj(WBz) - WJz.*conj(WBy))./fkPower;
ky = -1i*mu0*(WJz.*conj(WBx) - WJx.*conj(WBz))./fkPower;
kz = -1i*mu0*(WJx.*conj(WBy) - WJy.*conj(WBx))./fkPower;

% Wave vectors are generally complex. Here only the real part is taken. 
kx = real(kx);
ky = real(ky);
kz = real(kz);

kmag = sqrt(kx.*kx + ky.*ky + kz.*kz);

Bbackxmat = Bback.x.data*ones(1,numf);
Bbackymat = Bback.y.data*ones(1,numf);
Bbackzmat = Bback.z.data*ones(1,numf);
Bbackabsmat = Bback.abs.data*ones(1,numf);

kpar = (kx.*Bbackxmat + ky.*Bbackymat + kz.*Bbackzmat)./Bbackabsmat;
kperp = sqrt(kmag.^2 - kpar.^2);

kmax = max(max(kmag))*1.1;
kmin = -kmax;
kvec = linspace(-kmax,kmax,numk);
kmagvec = linspace(0,kmax,numk);

dkmag = kmax/numk;
dk = 2*kmax/numk;

if plotEpower
  Powerp = fkPowerE;
else
  Powerp = fkPower;
end

% Sort power into frequency and wave vector
irf.log('notice','Computing power versus kmag,f')
powerkmagf = zeros(numf,numk);
for mm = 1:L
	for nn = 1:numf
    knumber = floor((kmag(mm,nn))/dkmag)+1;
    powerkmagf(nn,knumber) = powerkmagf(nn,knumber) + Powerp(mm,nn);
  end
end

powerkmagf(powerkmagf == 0) = NaN;
%powerkmagf = powerkmagf/(N+1); % Normalization. This should correspond to FFT PSD when summed over k.
powerkmagf = powerkmagf/max(max(powerkmagf)); % Normalization to Max value for plotting. 

xvec1 = kmagvec;
yvec1 = WJ.f;

idxf = 1:numf;

if numel(frange)==2
  freqind = yvec1 > min(frange) & yvec1 < max(frange);
  idxf = idxf(freqind);
end

irf.log('notice','Computing power versus kperp,kpar')
powerkperpkpar = zeros(numk,numk);
for mm = 1:L
	for nn = idxf
    kparnumber = floor((kpar(mm,nn)-kmin)/dk)+1;
    kperpnumber = floor((kperp(mm,nn))/dkmag)+1;
    powerkperpkpar(kparnumber,kperpnumber) = powerkperpkpar(kparnumber,kperpnumber) + Powerp(mm,nn);
  end
end

powerkperpkpar(powerkperpkpar == 0) = NaN;
powerkperpkpar = powerkperpkpar/max(max(powerkperpkpar));

xvec2 = kmagvec;
yvec2 = kvec;

irf.log('notice','Computing power versus kx,ky')
powerkxky = zeros(numk,numk);
for mm = 1:L
	for nn = idxf
    kxnumber = floor((kx(mm,nn)-kmin)/dk)+1;
    kynumber = floor((ky(mm,nn)-kmin)/dk)+1;
    powerkxky(kynumber,kxnumber) = powerkxky(kynumber,kxnumber) + Powerp(mm,nn);
  end
end

powerkxky(powerkxky == 0) = NaN;
powerkxky = powerkxky/max(max(powerkxky));

xvec3 = kvec;
yvec3 = kvec;

irf.log('notice','Computing power versus kx,ky')
powerkxkz = zeros(numk,numk);
for mm = 1:L
	for nn = idxf
    kxnumber = floor((kx(mm,nn)-kmin)/dk)+1;
    kznumber = floor((kz(mm,nn)-kmin)/dk)+1;
    powerkxkz(kznumber,kxnumber) = powerkxkz(kznumber,kxnumber) + Powerp(mm,nn);
  end
end

powerkxkz(powerkxkz == 0) = NaN;
powerkxkz= powerkxkz/max(max(powerkxkz));

xvec4 = kvec;
yvec4 = kvec;

irf.log('notice','Computing power versus ky,kz')
powerkykz = zeros(numk,numk);
for mm = 1:L
	for nn = idxf
    kynumber = floor((ky(mm,nn)-kmin)/dk)+1;
    kznumber = floor((kz(mm,nn)-kmin)/dk)+1;
    powerkykz(kznumber,kynumber) = powerkykz(kznumber,kynumber) + Powerp(mm,nn);
  end
end

powerkykz(powerkykz == 0) = NaN;
powerkykz= powerkykz/max(max(powerkykz));

xvec5 = kvec;
yvec5 = kvec;

irf.log('notice','Computing power versus kx,f')
powerkxf = zeros(numf,numk);
for mm = 1:L
	for nn = 1:numf
    kxnumber = floor((kx(mm,nn)-kmin)/dk)+1;
    powerkxf(nn,kxnumber) = powerkxf(nn,kxnumber) + Powerp(mm,nn);
  end
end

powerkxf(powerkxf == 0) = NaN;
powerkxf = powerkxf/max(max(powerkxf));

xvec6 = kvec;
yvec6 = WJ.f;

irf.log('notice','Computing power versus ky,f')
powerkyf = zeros(numf,numk);
for mm = 1:L
	for nn = 1:numf
    kynumber = floor((ky(mm,nn)-kmin)/dk)+1;
    powerkyf(nn,kynumber) = powerkyf(nn,kynumber) + Powerp(mm,nn);
  end
end

powerkyf(powerkyf == 0) = NaN;
powerkyf = powerkyf/max(max(powerkyf));

xvec7 = kvec;
yvec7 = WJ.f;

irf.log('notice','Computing power versus kz,f')
powerkzf = zeros(numf,numk);
for mm = 1:L
	for nn = 1:numf
    kznumber = floor((kz(mm,nn)-kmin)/dk)+1;
    powerkzf(nn,kznumber) = powerkzf(nn,kznumber) + Powerp(mm,nn);
  end
end

powerkzf(powerkzf == 0) = NaN;
powerkzf = powerkzf/max(max(powerkzf));

xvec8 = kvec;
yvec8 = WJ.f;


xvariable = struct('kmag',xvec1,'kperp',xvec2,'kxkxky',xvec3,'kxkxkz',xvec4,'kykykz',xvec5,'kxf',xvec6,'kyf',xvec7,'kzf',xvec8);
yvariable = struct('fkmag',yvec1,'kpar',yvec2,'kykxky',yvec3,'kzkxkz',yvec4,'kzkykz',yvec5,'fkxf',yvec6,'fkyf',yvec7,'fkzf',yvec8);
powerxy = struct('Powerkmagf',powerkmagf,'Powerkperpkpar',powerkperpkpar,'Powerkxky',powerkxky,'Powerkxkz',powerkxkz,'Powerkykz',powerkykz,'Powerkxf',powerkxf,'Powerkyf',powerkyf,'Powerkzf',powerkzf);

end


