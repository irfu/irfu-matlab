function [fkpower,freq,wavenumber] = fk_powerspectrum(varargin)
%
% [fkpower,freq,wavenumber] = mms.fk_powerspectrum(probecomb,trange,V6,Bxyz,zphase)
%
% Function to calculate the frequency-wave number power spectrum using
% MMS's four spin plane probes. Follows the same procedure as
% c_wk_powerspec. Function uses L2 scpot probe potential data to construct
% electric fields aligned with B separated by 60 m. Wavelet based
% cross-spectral analysis is used to calculate the phase difference between
% and the fields, and hence the wave number. The power is then binned
% according to frequency and wave number (Graham et al., JGR, 2016).
% Written by D. B. Graham.
%
% Input: (All data must be in TSeries format)
%       probecomb - Probe combination to use (1 or 3 or 5). 1 for B aligned with
%                   probes 1 and 2. 3 for B aligned with probes 3 and 4. 5
%                   for B aligned with 5 and 6.
%       trange -    time interval over which the power spectrum is calculated.
%                   B should be closely aligned with one probe pair over this time.
%                   See mms.probe_align_times
%       V6 -        L2 probe potentials. Timing corrections are applied in this
%                   function.   Do not apply them before running this function.
%       Bxyz -      Magnetic field in DMPA coordinates.
%       zphase -    Spacecraft phase (zphase). Obtained from ancillary_defatt.
%                   Needed only if probecomb = 1 or 3.
%
% Options:
%       cav -       Number of points in timeseries used to estimate phase.
%                   Optional parameter. Default is cav = 128;
%       field -     Set to 1 (default) to use electric fields calculated from
%                   opposing probes and SCpot. Set to 0 to use only
%                   opposing probe potentials.
%       numk -      Set number of wave numbers used in spectrogram.
%       linear -    Linearly spaced frequencies. Set number to df (default is logarithmic spacing).
%       numf -     Set number of frequencies used in spectrogram.
%
% Output:
%       fkpower    - array of powers as a function of frequency and
%                    wavenumber. Power is normalized to the maximum value.
%       freq       - array of frequencies f in Hz.
%       wavenumber - array of wavenumbers k in m^{-1}. Positive values are
%                    aligned with B and negative values are anti-aligned with B.
%
% Example:
%   [fkpower,freq,wavenumber] = mms.fk_powerspectrum(probecomb,trange,V6,Bxyz,zphase,'cav',256,'field',0);
%   [fkpower,freq,wavenumber] = mms.fk_powerspectrum(probecomb,trange,V6,Bxyz,zphase,'cav',256,'linear',50,'numk',500);
%
% Directions and speeds are consistent with expectations based on time
% delays. Work still in progress. Time corrections for the probes need to be
% revised.

use_56=(varargin{1}==5);

if ((numel(varargin) < 5 && ~use_56)||(numel(varargin) < 4))
  help mms.fk_powerspectrum;
  return;
end

probe = varargin{1};
ts2=varargin{2};
SCpot = varargin{3};
Bxyz = varargin{4};

if ~use_56
  zphase = varargin{5};
  args=varargin(6:end);
else
  args=varargin(5:end);
end

cav = 128;
numk = 101;
numf = 200;
fieldflag = 1;
uselinear = 0;


if numel(args)>0
  flag_have_options=1;
  irf.log('notice','have options');
else
  flag_have_options=0;
end

while flag_have_options
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
        df = args{2}; %#ok<NASGU>
        uselinear = 1;
        irf.log('notice','Using linearly spaced frequencies');
      end
    case 'field'
      if numel(args)>1 && isnumeric(args{2})
        if args{2}
          fieldflag = 1;
          irf.log('notice','Using reconstructed electric field for computation.')
        else
          fieldflag = 0;
          irf.log('notice','Using probe potentials for computation.')
        end
      end
    otherwise
      irf.log('warning',['Unknown flag: ' args{1}]);
      l=1;
      break;
  end
  args = args(l+1:end);
  if isempty(args), flag_have_options=0; end
end


% Correct for timing in spacecraft potential data.
E12 = TSeries(SCpot.time,(SCpot.data(:,1)-SCpot.data(:,2))/0.120);
E34 = TSeries(SCpot.time,(SCpot.data(:,3)-SCpot.data(:,4))/0.120);
E56 = TSeries(SCpot.time,(SCpot.data(:,5)-SCpot.data(:,6))/0.0292);
V1 = TSeries(SCpot.time,SCpot.data(:,1));
V3 = TSeries(SCpot.time+ 7.629e-6,SCpot.data(:,3));
V5 = TSeries(SCpot.time+15.259e-6,SCpot.data(:,5));
E12.time = E12.time + 26.703e-6;
E34.time = E34.time + 30.518e-6;
E56.time = E56.time + 34.332e-6;

test_resamp=1;
if test_resamp
  V3 = mms.dft_timeshift(V3,-7.629e-6);
  V5 = mms.dft_timeshift(V5,-15.259e-6);
  E12 = mms.dft_timeshift(E12,-26.703e-6);
  E34 = mms.dft_timeshift(E34,-30.518e-6);
  E56 = mms.dft_timeshift(E56,-34.332e-6);
else
  V3 = V3.resample(V1.time); %#ok<UNRCH>
  V5 = V5.resample(V1.time);
  E12 = E12.resample(V1.time); %These resamples need to be changed.
  E34 = E34.resample(V1.time);
  E56 = E56.resample(V1.time);

end

V2 = V1 - E12 * 0.120;
V4 = V3 - E34 * 0.120;
V6 = V5 - E56 * 0.0292;
% Make new SCpot with corrections
SCpot = irf.ts_scalar(V1.time,[V1.data V2.data V3.data V4.data V5.data V6.data]);



ts2l = ts2+[-1 1];
SCpot = SCpot.tlim(ts2l);
Bxyz = Bxyz.tlim(ts2l);

if ~use_56 %Unnecessary if probes 5 & 6 are used.

  zphase = zphase.tlim(ts2l);

  norepeat = ones(length(zphase.time),1);
  nph = length(zphase.data);
  for ii=2:nph
    if(zphase.time(ii) > zphase.time(ii-1))
      if(zphase.data(ii) < zphase.data(ii-1))
        zphase.data(ii:end) = zphase.data(ii:end)+double(360.0);
      end
    else
      norepeat(ii) = 0;
    end
  end

  zphasetime = zphase.time(norepeat == 1);
  zphasedata = zphase.data(norepeat == 1);

  zphase = TSeries(zphasetime,zphasedata,'to',1);

  zphase = zphase.resample(SCpot);
end
Bxyz = Bxyz.resample(SCpot);

%Construct fields or potentials (N.B. Amplitude is not important).
time = SCpot.time;
SCV12 = (SCpot.data(:,1)+SCpot.data(:,2))/2;
SCV34 = (SCpot.data(:,3)+SCpot.data(:,4))/2;

if ~fieldflag
  SCV12 = zeros(size(SCV12));
  SCV34 = zeros(size(SCV34));
end

E1 = (SCpot.data(:,1)-SCV34)*1e3/60; %#ok<NASGU>
E2 = (SCV34-SCpot.data(:,2))*1e3/60; %#ok<NASGU>
E3 = (SCpot.data(:,3)-SCV12)*1e3/60; %#ok<NASGU>
E4 = (SCV12-SCpot.data(:,4))*1e3/60; %#ok<NASGU>
E5 = (SCpot.data(:,5)-(SCV34+SCV12)/2)*1e3/14.6; %#ok<NASGU> %Added
E6 = ((SCV34+SCV12)/2-SCpot.data(:,6))*1e3/14.6; %#ok<NASGU> %Added

c_eval('E? = TSeries(time,E?,''to'',1);',1:6); %Generalized

if ~use_56 %Added third dimension to work with 5-6, does not affect 12,34.
  %Get spacecraft phase
  phase_p1=zphase.data/180*pi + pi/6; %#ok<NASGU>
  phase_p3=zphase.data/180*pi + 2*pi/3; %#ok<NASGU>
  phase_p2=zphase.data/180*pi + 7*pi/6; %#ok<NASGU>
  phase_p4=zphase.data/180*pi + 5*pi/3; %#ok<NASGU>
  c_eval('rp?=[60*cos(phase_p?) 60*sin(phase_p?), zeros(length(phase_p?),1)];', 1:4);
  probe_nr = [1 3];
else
  probe_nr = 5;
  phase_p5=ones(length(Bxyz.data(:,1)),1); %#ok<NASGU>
  phase_p6=-1*ones(length(Bxyz.data(:,1)),1); %#ok<NASGU>
  c_eval('rp?=[zeros(length(Bxyz.data(:,1)),1) zeros(length(Bxyz.data(:,1)),1) phase_p?*14.6];', 5:6);
end

c_eval('thetap?b = (rp?(:,1).*Bxyz.data(:,1)+rp?(:,2).*Bxyz.data(:,2)+rp?(:,3).*Bxyz.data(:,3))./(sqrt(rp?(:,1).^2+rp?(:,2).^2+rp?(:,3).^2).*Bxyz.abs.data);',probe_nr);
idx = tlim(SCpot.time,ts2);

% If odd, remove last data point (as is done in irf_wavelet)
if mod(length(idx),2)
  idx(end)=[];
end

c_eval('thetatest = irf.nanmean(thetap?b(idx));',probe);
c_eval('thetap?b = abs(thetap?b);',probe_nr);
c_eval('thetap?b = TSeries(Bxyz.time,thetap?b,''to'',1);',probe_nr);

if ~uselinear
  if (thetatest > 0)
    c_eval('W1c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''nf'',numf);',probe+1);
    c_eval('W2c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''nf'',numf);',probe);
  else
    c_eval('W1c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''nf'',numf);',probe);
    c_eval('W2c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''nf'',numf);',probe+1);
  end
else
  if (thetatest > 0)
    c_eval('W1c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''linear'',df);',probe+1);
    c_eval('W2c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''linear'',df);',probe);
  else
    c_eval('W1c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''linear'',df);',probe);
    c_eval('W2c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0,''linear'',df);',probe+1);
  end
  numf = length(W1c.f); %#ok<NODEF>
end

L = length(idx);
times = time(idx);

W1c.p = {W1c.p{1,1}(idx,:)};
W1c.t = times;
W2c.p = {W2c.p{1,1}(idx,:)}; %#ok<NODEF>
W2c.t = times;

fkPower = 0.5*(cell2mat(W1c.p).*conj(cell2mat(W1c.p)) + cell2mat(W2c.p).*conj(cell2mat(W2c.p)));

N = floor(L/cav)-1;
posav = cav/2 + (0:1:N)*cav;
avtimes = times(posav);
Bs = Bxyz.resample(avtimes);
c_eval('thetap?b = thetap?b.resample(avtimes);',probe_nr);
if ~use_56
  c_eval('rcos = 60.0*thetap?b.data;',probe);
else
  c_eval('rcos = 14.6*thetap?b.data;',probe); %Added for 5-6.
end

c34x = zeros(N+1,numf);
Powerav = zeros(N+1,numf);

for m = 1:N+1
  c34x(m,:) = irf.nanmean(W1c.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W2c.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
  Powerav(m,:) = irf.nansum(fkPower([posav(m)-cav/2+1:posav(m)+cav/2],:));
end

cross34x = W1c;
cross34x.p = abs(c34x);
cross34x.t = avtimes;
cross34x.f = W1c.f;

tempr = real(c34x);
tempi = imag(c34x);
th = (atan2(tempi,tempr));
kval = zeros(N+1,numf);

for q = 1:numf
  kval(:,q) = th(:,q)./rcos;
end

disprel = zeros(numk,numf);
mink = -pi/min(rcos);
maxk = pi/min(rcos);
dk = (maxk - mink)/numk;
kvec = mink + [0:1:numk-1]*dk;

if 1
  for m = 1:N+1
    for q = 1:numf
      knumber = floor((kval(m,q)-mink)/dk)+1;
      disprel(knumber,q) = disprel(knumber,q) + Powerav(m,q);
    end
  end
end

disprel(disprel == 0) = NaN;
disprel = disprel/max(max(disprel));
disprel(disprel < 1.0e-3) = 1e-3;

wavenumber = kvec;
freq = cross34x.f;
fkpower = disprel';

end

