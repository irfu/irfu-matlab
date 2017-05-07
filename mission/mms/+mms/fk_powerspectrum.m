function [fkpower,freq,wavenumber] = fk_powerspectrum(varargin) 
%
% [fkpower,freq,wavenumber] = mms.fk_powerspectrum(SCpot,Bxyz,zphase,trange,SCnum,probecomb) 
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
%       SCpot -     L2 probe potentials. Timing corrections are applied in this
%       function. Do not apply them before running this function.
%       Bxyz -      Magnetic field in DMPA coordinates.
%       trange -    time interval over which the power spectrum is calculated. 
%                   B should be closely aligned with one probe pair over this time. 
%                   See mms.probe_align_times
%       zphase -    Spacecraft phase (zphase). Obtained from ancillary_defatt.
%       probecomb - Probe combination to use (1 or 3). 1 for B aligned with
%                   probes 1 and 2. 3 for B aligned with probes 3 and 4.
%
% Options: 
%       cav -       Number of points in timeseries used to estimate phase.
%                   Optional parameter. Default is cav = 128;
%       field -     Set to 1 (default) to use electric fields calculated from
%                   opposing probes and SCpot. Set to 0 to use only
%                   opposing probe potentials.
%
% Output: 
%       fkpower    - array of powers as a function of frequency and
%                    wavenumber. Power is normalized to the maximum value.
%       freq       - array of frequencies f in Hz.
%       wavenumber - array of wavenumbers k in m^{-1}. Positive values are
%                    aligned with B and negative values are anti-aligned with B.
% 
% Example:
%   [fkpower,freq,wavenumber] = mms.fk_powerspectrum(SCpot,Bxyz,zphase,trange,SCnum,probecomb,'cav',256,'field',0); 
%
% Directions and speeds are consistent with expectations based on time
% delays. Work still in progress. Time corrections for the probes need to be
% revised. 


if (numel(varargin) < 5)
    help mms.fk_powerspectrum;
    return;
end

SCpot = varargin{1};
Bxyz = varargin{2};
zphase = varargin{3};
ts2=varargin{4};
probe = varargin{5};
cav = 128;
fieldflag = 1;

args=varargin(6:end);
if numel(args)>0,
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
        case 'field'
            if numel(args)>1 && isnumeric(args{2})
                if args{2},
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
V3 = V3.resample(V1.time);
V5 = V5.resample(V1.time);
E12 = E12.resample(V1.time);
E34 = E34.resample(V1.time);
E56 = E56.resample(V1.time);
V2 = V1 - E12 * 0.120;
V4 = V3 - E34 * 0.120;
V6 = V5 - E56 * 0.0292;
% Make new SCpot with corrections
SCpot = irf.ts_scalar(V1.time,[V1.data V2.data V3.data V4.data V5.data V6.data]);

ts2l = ts2+[-1 1];
SCpot = SCpot.tlim(ts2l);
Bxyz = Bxyz.tlim(ts2l);

zphase = zphase.tlim(ts2l);

norepeat = ones(length(zphase.time),1);
nph = length(zphase.data);
for ii=[2:nph]
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
Bxyz = Bxyz.resample(SCpot);

%Construct fields or potentials (N.B. Amplitude is not important).
time = SCpot.time;
SCV12 = (SCpot.data(:,1)+SCpot.data(:,2))/2;
SCV34 = (SCpot.data(:,3)+SCpot.data(:,4))/2;
if ~fieldflag
    SCV12 = zeros(size(SCV12));
    SCV34 = zeros(size(SCV34));
end
  
E1 = (SCpot.data(:,1)-SCV34)*1e3/60;
E2 = (SCV34-SCpot.data(:,2))*1e3/60;
E3 = (SCpot.data(:,3)-SCV12)*1e3/60;
E4 = (SCV12-SCpot.data(:,4))*1e3/60;
c_eval('E? = TSeries(time,E?,''to'',1);',[1:4]);

%Get spacecraft phase
phase_p1=zphase.data/180*pi + pi/6;
phase_p3=zphase.data/180*pi + 2*pi/3;
phase_p2=zphase.data/180*pi + 7*pi/6;
phase_p4=zphase.data/180*pi + 5*pi/3;
c_eval('rp?=[60*cos(phase_p?) 60*sin(phase_p?)];',[1:4]);

c_eval('thetap?b = (rp?(:,1).*Bxyz.data(:,1)+rp?(:,2).*Bxyz.data(:,2))./(sqrt(rp?(:,1).^2+rp?(:,2).^2).*Bxyz.abs.data);',[1 3]);
idx = tlim(SCpot.time,ts2);
c_eval('thetatest = irf.nanmean(thetap?b(idx));',probe);
c_eval('thetap?b = abs(thetap?b);',[1 3]);
c_eval('thetap?b = TSeries(Bxyz.time,thetap?b,''to'',1);',[1 3]);

if (thetatest > 0)
    c_eval('W1c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0);',probe+1);
    c_eval('W2c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0);',probe);  
else
  c_eval('W1c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0);',probe);
  c_eval('W2c = irf_wavelet(E?,''returnpower'',0,''cutedge'',0);',probe+1);
end
   
L = length(idx);
times = time(idx);

W1c.p = {W1c.p{1,1}(idx,:)};
W1c.t = times;
W2c.p = {W2c.p{1,1}(idx,:)};
W2c.t = times;

fkPower = 0.5*(cell2mat(W1c.p).*conj(cell2mat(W1c.p)) + cell2mat(W2c.p).*conj(cell2mat(W2c.p)));

numf = 200;
N = floor(L/cav)-1;
posav = cav/2 + [0:1:N]*cav;
avtimes = times(posav);
Bs = Bxyz.resample(avtimes);
c_eval('thetap?b = thetap?b.resample(avtimes);',[1 3]);
c_eval('rcos = 60.0*thetap?b.data;',probe);

c34x = zeros(N+1,numf);
Powerav = zeros(N+1,numf);

for m = [1:1:N+1]
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
     
for q = [1:1:numf]
     kval(:,q) = th(:,q)./rcos;
end

numk = 101;

disprel = zeros(numk,numf);
mink = -pi/min(rcos);
maxk = pi/min(rcos);
dk = (maxk - mink)/numk; 
kvec = mink + [0:1:numk-1]*dk;

if 1
for m = [1:1:N+1]
    for q = [1:1:numf]
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

