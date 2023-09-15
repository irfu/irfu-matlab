function [xvariable,yvariable,powerxy] = fk_powerspec4SC(varargin)
%
% [fkpower,freq,wavenumber] = mms.fk_powerspec4SC('E?','R?','B?',Tints,options...)
%
% Function to calculate the frequency-wave number power spectrum using
% the four MMS spacecraft. Uses a generalization of mms.fk_powerspectrum.
% Wavelet based cross-spectral analysis is used to calculate the phase
% difference each spacecraft pair and determine 3D wave vector. A
% generalization of the method used in mms.fk_powerspectrum to four point
% measurements.
% Written by D. B. Graham.
%
% Input: (All data must be in TSeries format)
%       E? -        Fields to apply 4SC cross-spectral analysis to. E.g., E
%                   or B fields (if multiple components only the first is used).
%       R? -        Positions of the four spacecraft
%       B? -        Background magnetic field in the same coordinates as R?.
%                   Used to determine the parallel and perpendicular wave numebers using 4SC average.
%       Tints -     Time interval over which the power spectrum is calculated.
%                   To avoid boundary effects use a longer time interval
%                   for E? and B?.
%
% Options:
%       cav -       Number of points in timeseries used to estimate phase.
%                   Optional parameter. Default is cav = 8;
%       numk -      Set number of wave numbers used in spectrogram.
%       linear -    Linearly spaced frequencies. Set number to df (default is logarithmic spacing).
%       numf -      Set number of frequencies used in spectrogram.
%       wwidth -    Multiplier for Morlet wavelet width. Default is 1.
%       linear -    Use linear spacing between frequencies. df argument is
%                   required.
%       frange -    select frequency range for k-k plots. [minf maxf]
%
% Output:
%       powerxy    - array of powers as a function of frequency and
%                    wavenumber. Power is normalized to the maximum value.
%       yvariable       - array of frequencies f in Hz.
%       xvariable - array of wavenumbers k in m^{-1}. Positive values are
%                    aligned with B and negative values are anti-aligned with B.
%
% Notes:
%   Wavelength must be larger than twice the spacecraft separations,
%   otherwise spatial aliasing will occur.
%
%
% Example:
%   [xvecs,yvecs,Power] = fk_powerspec4SC('Epar?','Rxyz?','Bxyz?',Tints)
%   [xvecs,yvecs,Power] = fk_powerspec4SC('Bscmfac?.x','Rxyz?','Bxyz?',Tints,'linear',10,'numk',500,'cav',4,'wwidth',2);
%
% Example to plot:
%   pcolor(xvecs.kmag,yvecs.fkmag*1e-3,log10(Power.Powerkmagf));
%   shading('flat');
%   xlabel('|k| (m^{-1})');
%   ylabel('f (kHz)');
%


if (nargin < 4)
  help mms.fk_powerspec4SC;
  powerxy = NaN;
  xvariable = NaN;
  yvariable = NaN;
  return;
end

ic = 1:4;

c_eval('E?=evalin(''base'',irf_ssub(varargin{1},?));',ic);
c_eval('R?=evalin(''base'',irf_ssub(varargin{2},?));',ic);
c_eval('B?=evalin(''base'',irf_ssub(varargin{3},?));',ic);
Tint = varargin{4};
c_eval('B? = B?.resample(B1);',2:4);
Bav = irf.ts_vec_xyz(B1.time,(B1.data+B2.data+B3.data+B4.data)/4);

c_eval('E? = E?.resample(E1);',2:4);
c_eval('R? = R?.resample(E1);',1:4);
time = E1.time;

cav = 8;
numk = 500;
numf = 200;
uselinear = 0;
wwidth = 1;
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
    case 'cav'
      if numel(args)>1 && isnumeric(args{2})
        cav = args{2};
      end
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

idx = tlim(E1.time,Tint);

% If odd, remove last data point (as is done in irf_wavelet)
if mod(length(idx),2)
  idx(end)=[];
end

if ~uselinear
  c_eval('W? = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''nf'',numf,''wavelet_width'',5.36*wwidth);',ic);
else
  c_eval('W? = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''linear'',df,''wavelet_width'',5.36*wwidth);',ic);
end
numf = length(W1.f);

L = length(idx);
times = time(idx);

c_eval('W?.p = {W?.p{1,1}(idx,:)};',ic);
c_eval('W?.t = times;',ic);

fkPower = 0.25*(cell2mat(W1.p).*conj(cell2mat(W1.p)) + cell2mat(W2.p).*conj(cell2mat(W2.p)) ...
  + cell2mat(W3.p).*conj(cell2mat(W3.p)) + cell2mat(W4.p).*conj(cell2mat(W4.p)));

N = floor(L/cav)-1;
posav = cav/2 + (0:1:N)*cav;
avtimes = times(posav);

Bav = Bav.resample(avtimes);
c_eval('R? = R?.resample(avtimes);',ic);

cx12 = zeros(N+1,numf);
cx13 = zeros(N+1,numf);
cx14 = zeros(N+1,numf);
cx23 = zeros(N+1,numf);
cx24 = zeros(N+1,numf);
cx34 = zeros(N+1,numf);
Powerav = zeros(N+1,numf);

for m = 1:N+1
  cx12(m,:) = irf.nanmean(W1.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W2.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  cx13(m,:) = irf.nanmean(W1.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W3.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  cx14(m,:) = irf.nanmean(W1.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W4.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  cx23(m,:) = irf.nanmean(W2.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W3.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  cx24(m,:) = irf.nanmean(W2.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W4.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  cx34(m,:) = irf.nanmean(W3.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W4.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  Powerav(m,:) = irf.nanmean(fkPower([posav(m)-cav/2+1:posav(m)+cav/2],:));
end

% Compute phase differences between each spacecraft pair
th12 = atan2(imag(cx12),real(cx12));
th13 = atan2(imag(cx13),real(cx13));
th14 = atan2(imag(cx14),real(cx14));
th23 = atan2(imag(cx23),real(cx23));
th24 = atan2(imag(cx24),real(cx24));
th34 = atan2(imag(cx34),real(cx34));

wmat = 2*pi*ones(N+1,1)*(W1.f)';

% Convert phase difference to time delay
dt12 = th12./wmat;
dt13 = th13./wmat;
dt14 = th14./wmat;
dt23 = th23./wmat;
dt24 = th24./wmat;
dt34 = th34./wmat;

% Weighted averaged time delay using all spacecraft pairs
dt2 = 0.5*dt12 + 0.2*(dt13 - dt23) + 0.2*(dt14 - dt24) + 0.1*(dt14 - dt34 - dt23);
dt3 = 0.5*dt13 + 0.2*(dt12 + dt23) + 0.2*(dt14 - dt34) + 0.1*(dt12 + dt24 - dt34);
dt4 = 0.5*dt14 + 0.2*(dt12 + dt24) + 0.2*(dt13 + dt34) + 0.1*(dt12 + dt23 + dt34);
%dt2 = dt12;
%dt3 = dt13;
%dt4 = dt14;

% Compute phase speeds

R1 = R1.data;
R2 = R2.data;
R3 = R3.data;
R4 = R4.data;

kx = zeros(N+1,numf);
ky = zeros(N+1,numf);
kz = zeros(N+1,numf);

for ii = 1:N+1
  dR = [R2(ii,:);R3(ii,:);R4(ii,:)]-[R1(ii,:);R1(ii,:);R1(ii,:)];
  for jj = 1:numf
    % Volumetric tensor with SC1 as center.
    m = dR\[dt2(ii,jj);dt3(ii,jj);dt4(ii,jj)]; % "1/v vector"

    kx(ii,jj) = 2*pi*W1.f(jj)*m(1);
    ky(ii,jj) = 2*pi*W1.f(jj)*m(2);
    kz(ii,jj) = 2*pi*W1.f(jj)*m(3);
  end
end

kx = kx/1e3;
ky = ky/1e3;
kz = kz/1e3;
kmag = sqrt(kx.*kx + ky.*ky + kz.*kz);

Bavxmat = Bav.x.data*ones(1,numf);
Bavymat = Bav.y.data*ones(1,numf);
Bavzmat = Bav.z.data*ones(1,numf);
Bavabsmat = Bav.abs.data*ones(1,numf);

kpar = (kx.*Bavxmat + ky.*Bavymat + kz.*Bavzmat)./Bavabsmat;
kperp = sqrt(kmag.^2 - kpar.^2);

kmax = max(max(kmag))*1.1;
kmin = -kmax;
kvec = linspace(-kmax,kmax,numk);
kmagvec = linspace(0,kmax,numk);

dkmag = kmax/numk;
dk = 2*kmax/numk;

% Sort power into frequency and wave vector
irf.log('notice','Computing power versus kmag,f')
powerkmagf = zeros(numf,numk);
for mm = 1:N+1
  for nn = 1:numf
    knumber = floor((kmag(mm,nn))/dkmag)+1;
    powerkmagf(nn,knumber) = powerkmagf(nn,knumber) + Powerav(mm,nn);
  end
end

powerkmagf(powerkmagf == 0) = NaN;
%powerkmagf = powerkmagf/(N+1); % Normalization. This should correspond to FFT PSD when summed over k.
powerkmagf = powerkmagf/max(max(powerkmagf)); % Normalization to Max value for plotting.
%powerkmagf(powerkmagf < 1.0e-6) = 1e-6;

xvec1 = kmagvec;
yvec1 = W1.f;

idxf = 1:numf;

if numel(frange)==2
  freqind = yvec1 > min(frange) & yvec1 < max(frange);
  idxf = idxf(freqind);
end
idxf

irf.log('notice','Computing power versus kperp,kpar')
powerkperpkpar = zeros(numk,numk);
for mm = 1:N+1
  for nn = idxf
    kparnumber = floor((kpar(mm,nn)-kmin)/dk)+1;
    kperpnumber = floor((kperp(mm,nn))/dkmag)+1;
    powerkperpkpar(kparnumber,kperpnumber) = powerkperpkpar(kparnumber,kperpnumber) + Powerav(mm,nn);
  end
end

powerkperpkpar(powerkperpkpar == 0) = NaN;
powerkperpkpar = powerkperpkpar/max(max(powerkperpkpar));
%powerkperpkpar(powerkperpkpar < 1.0e-6) = 1e-6;

xvec2 = kmagvec;
yvec2 = kvec;

irf.log('notice','Computing power versus kx,ky')
powerkxky = zeros(numk,numk);
for mm = 1:N+1
  for nn = idxf
    kxnumber = floor((kx(mm,nn)-kmin)/dk)+1;
    kynumber = floor((ky(mm,nn)-kmin)/dk)+1;
    powerkxky(kynumber,kxnumber) = powerkxky(kynumber,kxnumber) + Powerav(mm,nn);
  end
end

powerkxky(powerkxky == 0) = NaN;
powerkxky = powerkxky/max(max(powerkxky));
%powerkxky(powerkxky < 1.0e-6) = 1e-6;

xvec3 = kvec;
yvec3 = kvec;

irf.log('notice','Computing power versus kx,ky')
powerkxkz = zeros(numk,numk);
for mm = 1:N+1
  for nn = idxf
    kxnumber = floor((kx(mm,nn)-kmin)/dk)+1;
    kznumber = floor((kz(mm,nn)-kmin)/dk)+1;
    powerkxkz(kznumber,kxnumber) = powerkxkz(kznumber,kxnumber) + Powerav(mm,nn);
  end
end

powerkxkz(powerkxkz == 0) = NaN;
powerkxkz= powerkxkz/max(max(powerkxkz));
%powerkxkz(powerkxkz < 1.0e-6) = 1e-6;

xvec4 = kvec;
yvec4 = kvec;

irf.log('notice','Computing power versus ky,kz')
powerkykz = zeros(numk,numk);
for mm = 1:N+1
  for nn = idxf
    kynumber = floor((ky(mm,nn)-kmin)/dk)+1;
    kznumber = floor((kz(mm,nn)-kmin)/dk)+1;
    powerkykz(kznumber,kynumber) = powerkykz(kznumber,kynumber) + Powerav(mm,nn);
  end
end

powerkykz(powerkykz == 0) = NaN;
powerkykz= powerkykz/max(max(powerkykz));
%powerkykz(powerkykz < 1.0e-6) = 1e-6;

xvec5 = kvec;
yvec5 = kvec;

irf.log('notice','Computing power versus kx,f')
powerkxf = zeros(numf,numk);
for mm = 1:N+1
  for nn = 1:numf
    kxnumber = floor((kx(mm,nn)-kmin)/dk)+1;
    powerkxf(nn,kxnumber) = powerkxf(nn,kxnumber) + Powerav(mm,nn);
  end
end

powerkxf(powerkxf == 0) = NaN;
powerkxf = powerkxf/max(max(powerkxf));
%powerkxf(powerkxf < 1.0e-6) = 1e-6;

xvec6 = kvec;
yvec6 = W1.f;

irf.log('notice','Computing power versus ky,f')
powerkyf = zeros(numf,numk);
for mm = 1:N+1
  for nn = 1:numf
    kynumber = floor((ky(mm,nn)-kmin)/dk)+1;
    powerkyf(nn,kynumber) = powerkyf(nn,kynumber) + Powerav(mm,nn);
  end
end

powerkyf(powerkyf == 0) = NaN;
powerkyf = powerkyf/max(max(powerkyf));
%powerkyf(powerkyf < 1.0e-6) = 1e-6;

xvec7 = kvec;
yvec7 = W1.f;

irf.log('notice','Computing power versus kz,f')
powerkzf = zeros(numf,numk);
for mm = 1:N+1
  for nn = 1:numf
    kznumber = floor((kz(mm,nn)-kmin)/dk)+1;
    powerkzf(nn,kznumber) = powerkzf(nn,kznumber) + Powerav(mm,nn);
  end
end

powerkzf(powerkzf == 0) = NaN;
powerkzf = powerkzf/max(max(powerkzf));
%powerkzf(powerkzf < 1.0e-6) = 1e-6;

xvec8 = kvec;
yvec8 = W1.f;


xvariable = struct('kmag',xvec1,'kperp',xvec2,'kxkxky',xvec3,'kxkxkz',xvec4,'kykykz',xvec5,'kxf',xvec6,'kyf',xvec7,'kzf',xvec8);
yvariable = struct('fkmag',yvec1,'kpar',yvec2,'kykxky',yvec3,'kzkxkz',yvec4,'kzkykz',yvec5,'fkxf',yvec6,'fkyf',yvec7,'fkzf',yvec8);
powerxy = struct('Powerkmagf',powerkmagf,'Powerkperpkpar',powerkperpkpar,'Powerkxky',powerkxky,'Powerkxkz',powerkxkz,'Powerkykz',powerkykz,'Powerkxf',powerkxf,'Powerkyf',powerkyf,'Powerkzf',powerkzf);

end


