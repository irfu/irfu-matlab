function [power,freq,wavenumber] = c_wk_powerspec(VL1,B,trange, SCnum, probecomb) 
%
% [power,freq,wavenumber] = wkpowerspec(VL1,B,trange, SCnum, probecomb) 
%
% Function to calculate the frequency-wave number power spectrum using
% EFW's internal burst mode data. Written by D. B. Graham.
%
% Function uses L1 EFW data to construct electric fields aligned with B
% separated by 44m. Wavelet based cross-spectral analysis is used to
% calculate the phase difference between and the fields, and hence the wave
% number. The power is then binned according to frequency and wave number.
%
% Input: 
%       VL1 -       L1 probe potentials from EFW. Must be in the form 
%                   [tint Vp1 Vp2 Vp3 Vp4].
%       B -         Magnetic field data in GSE coordinates.
%       trange -    time interval over which the power spectrum is calculated. 
%                   B should be closely aligned with one probe pair over this time. 
%       SCnum -     Spacecraft number: 1-4.
%       probecomb - Probe combination to use (1 or 3). 1 for B aligned with
%                   probes 1 and 2. 3 for B aligned with probes 3 and 4.
%                   See also c_pl_sc_orient.
%
% Output: 
%       power      - array of powers as a function of frequency and
%                    wavenumber. Power is normalized to the maximum value.
%       freq       - array of frequencies f in Hz.
%       wavenumber - array of wavenumbers k in m^{-1}. Positive values are
%                    aligned with B and negative values are anti-aligned with B.
%
% Example:  
% tint = irf.tint('2005-02-25T10:37:11.5Z/2005-02-25T10:37:11.7Z');
% VL1 = c_caa_var_get('Data__C4_CP_EFW_L1_IB','caa','mat');
% B = c_caa_var_get('B_vec_xyz_gse__C4_CP_FGM_FULL','caa','ts');
% [power,freq,wavenumber] = wkpowerspec(VL1, B, tint, 4, 3);
%
% To Do: Update for MMS.

% If data is TSeries, data is converted to the older format

% Begin temporary fix to convert TS format to older format
if isa(VL1,'TSeries') 
    ttemp = VL1.time.epochUnix;
    datatemp = double(VL1.data);
    VL1 = [ttemp, double(datatemp)];
end
if isa(B,'TSeries') 
    ttemp = B.time.epochUnix;
    datatemp = double(B.data);
    B = [ttemp, datatemp];
end
if isa(trange,'EpochTT') 
    trange = trange.epochUnix;
end
% End of temporary fix

ic = SCnum;
probe = probecomb;

BDSC=c_coord_trans('GSE','DSC',B,'cl_id',ic);

ts2=trange;

%Construct fields
resc = 0.002165;
VL1 = [VL1(:,1) VL1(:,2)*resc VL1(:,3)*resc VL1(:,4)*resc VL1(:,5)*resc];
SCV12 = [VL1(:,1) (VL1(:,2)+VL1(:,3))/2];
SCV34 = [VL1(:,1) (VL1(:,4)+VL1(:,5))/2];
E1 = VL1(:,2)-SCV34(:,2);
E2 = SCV34(:,2)-VL1(:,3);
E3 = VL1(:,4)-SCV12(:,2);
E4 = SCV12(:,2)-VL1(:,5);
time = VL1(:,1);

%Get spacecraft phase
[tt,phase_data] = irf_isdat_get([irf_ssub('Cluster/?/ephemeris/phase_2',ic)], VL1(1,1)-10, VL1(end,1)-VL1(1,1)+20);
idx = find(phase_data == 0);
phasedata = [tt(idx) [0:1:length(idx)-1]'*360];
probett = VL1(:,1);
BDSC = irf_resamp(BDSC,probett);
phasedata = irf_resamp(phasedata,probett);
phasedata = [phasedata(:,1) mod(phasedata(:,2),360)];
phase_p1=phasedata(:,2)/180*pi + 3*pi/4 ;
phase_p3=phase_p1     - pi/2   ;
phase_p2=phase_p1     + pi     ;
phase_p4=phase_p1     + pi/2   ;
rp1=[44*cos(phase_p1) 44*sin(phase_p1)]; % in DSC reference frame
rp2=[44*cos(phase_p2) 44*sin(phase_p2)];
rp3=[44*cos(phase_p3) 44*sin(phase_p3)];
rp4=[44*cos(phase_p4) 44*sin(phase_p4)];
thetap1b = (rp1(:,1).*BDSC(:,2)+rp1(:,2).*BDSC(:,3))./(sqrt(rp1(:,1).^2+rp1(:,2).^2).*sqrt(BDSC(:,2).^2+BDSC(:,3).^2+BDSC(:,4).^2));
thetap3b = (rp3(:,1).*BDSC(:,2)+rp3(:,2).*BDSC(:,3))./(sqrt(rp3(:,1).^2+rp3(:,2).^2).*sqrt(BDSC(:,2).^2+BDSC(:,3).^2+BDSC(:,4).^2));
[~,idx]=irf_tlim(BDSC,ts2);
c_eval('thetatest = thetap?b(idx);',probe);
if (thetatest > 0)
    c_eval('W1c = irf_wavelet([time,E?],''returnpower'',0,''cutedge'',0);',probe+1);
    c_eval('W2c = irf_wavelet([time,E?],''returnpower'',0,''cutedge'',0);',probe);  
else
  c_eval('W1c = irf_wavelet([time,E?],''returnpower'',0,''cutedge'',0);',probe);
  c_eval('W2c = irf_wavelet([time,E?],''returnpower'',0,''cutedge'',0);',probe+1);
end
    
thetap1b = [VL1(:,1) abs(thetap1b)];
thetap3b = [VL1(:,1) abs(thetap3b)];



[~,idx]=irf_tlim(VL1,ts2);
L = length(idx);
times = time(idx);
E3s = E3(idx);
E4s = E4(idx);

W1c.p = {W1c.p{1,1}(idx,:)};
W1c.t = times;
W2c.p = {W2c.p{1,1}(idx,:)};
W2c.t = times;

Power = 0.5*(cell2mat(W1c.p).*conj(cell2mat(W1c.p)) + cell2mat(W2c.p).*conj(cell2mat(W2c.p)));

numf = 200;
cav = 128;
N = floor(L/cav)-1;
posav = cav/2 + (0:1:N)*cav;
avtimes = times(posav);
Bs = irf_resamp(B,avtimes);
thetap1b = irf_resamp(thetap1b,avtimes);
thetap3b = irf_resamp(thetap3b,avtimes);
c_eval('rcos = 44.0*thetap?b(:,2);',probe);

c34x = zeros(N+1,numf);
Powerav = zeros(N+1,numf);

for m = [1:1:N+1]
    c34x(m,:) = irf.nanmean(W1c.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W2c.p{1,1}([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    Powerav(m,:) = irf.nansum(Power([posav(m)-cav/2+1:posav(m)+cav/2],:));
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

numk = 101;

disprel = zeros(numk,numf);
mink = -pi/min(rcos);
maxk = pi/min(rcos);
dk = (maxk - mink)/numk; 
kvec = mink + [0:1:numk-1]*dk;

if 1
for m = [1:1:N+1]
    for q = 1:numf
        knumber = floor((kval(m,q)-mink)/dk)+1;
        disprel(knumber,q) = disprel(knumber,q) + Powerav(m,q);
    end
end
end       
 
disprel(find(disprel == 0)) = NaN;
disprel = disprel/max(max(disprel));
disprel(find(disprel < 1.0e-3)) = 1e-3;

wavenumber = kvec;
freq = cross34x.f;
power = disprel';

end

