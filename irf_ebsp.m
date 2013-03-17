function [outTime,frequencyVec,BB_xxyyzz_fac,...
    EESum_xxyy_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
    k_thphSVD_fac,polSVD_fac,ellipticity]=...
    irf_ebsp(e,dB,fullB,B0,xyz,freq_int,varargin)
%IRF_EBSP   Calculates E&B wavelet spectra, Poynting flux, polarization
%
% irf_ebsp(E,dB,fullB,B0,xyz,freq_int,[OPTIONS])
% modified from irf_pl_ebs
% assumes equidistant time spacing
%
% It uses a Morlet wavelet.
%
% Input:
%
% E = wave electric field, columns (t ex ey ez)
% dB = wave magnetic field, columns (t bx by bz)
% fullB = high resolution background magnetic field, columns (t bx by bz)
% B0 = background magnetic field, columns (t bx by bz)
% xyz = position vector of spacecraft, columns (t x y z)
% freq_int = frequency interval: either 'pc12', 'pc35' or numeric [fmin fmax]
%
% Options:
%   'noresamp' - no resampling, E and dB are given at the same timeline
%   'fac'      - use FAC coordinate system, otherwise no coordinate system 
%                transformation is performed
%   'dEdotB=0' - compute dEz from dB dot B = 0
%   'fullB=dB' - dB contains DC field
% 
% Output:
%
%     timeVector
%     frequencyVector
%     BVector
%     BB_xxyyzz_fac
%     EESum_xxyy_ISR2
%     EE_xxyyzz_FAC,...
%     Poynting_xyz_FAC
%     Poynting_rThetaPhi_FAC
%     k_thphSVD_fac
%     polSVD_fac
%     ellipticity
%
% Examples:
%
%    [timeVec,frequencyVector,BB_xxyyzz_fac]=...
%        irf_ebsp(e,b,B,xyz,'pc12');
%
%    [timeVec,frequencyVec,BB_xxyyzz_fac,...
%        EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,Poynting_xyz_FAC]=...
%        irf_ebsp(e,b,[],B,xyz,'pc35','fullB=dB');
%
%    [timeVec,frequencyVec,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,...
%        EE_xxyyzz_FAC,Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
%        k_thphSVD_fac,polSVD_fac,ellipticity]=...
%        irf_ebsp(e,b,[],B,xyz,'pc12','fullB=dB','dEdotB=0');
%
%  See also: IRF_PL_EBS, IRF_PL_EBSP

%% Check the input
nWavePeriodToAverage = 4; % Number of wave periods to average
angleBElevationMax = 15;  % Below which we cannot apply E*B=0

wantPolarization = 0;
if nargout==3,
    wantEE = 0;
elseif nargout==7,
    wantEE = 1;
elseif nargout==10,
    wantEE = 1;
    wantPolarization = 1;
else
	error('irf_ebsp: unknown number of output parameters');
end

flag_no_resamp = 0; flag_want_fac = 0; flag_dEdotB0 = 0; flag_fullB_dB = 0;
for i=1:length(varargin)
    switch lower(varargin{i})
        case 'noresamp'
            flag_no_resamp = 1;
        case 'fac'
            flag_want_fac = 1; % Use FAC coordinate system
        case 'dedotb=0'
            flag_dEdotB0 = 1;
        case 'fullb=db'
            flag_fullB_dB = 1;
        otherwise
            irf_log('fcal',['Option ''' varargin{i} '''not recognized'])
    end
end

if flag_want_fac && (isempty(B0) || isempty(xyz))
    error('B0 and XYZ must be given for option FAC')
end
B0 = irf_resamp(B0,dB);
if flag_fullB_dB
    fullB = dB;
    dB(:,2:4) = dB(:,2:4) - B0(:,2:4);
end
if flag_dEdotB0 && isempty(fullB)
    error('fullB must be given for option dEdotB=0')
end
Bx = fullB(:,2); By = fullB(:,3); Bz = fullB(:,4); % Needed for parfor
angleBElevation=atan2d(Bz,sqrt(Bx.^2+By.^2));
idxBparSpinPlane= abs(angleBElevation)<angleBElevationMax;
        
pc12_range=0; pc35_range=0; default_range=0;
if ischar(freq_int)
    switch lower(freq_int)
        case {'pc12'}
            freq_int=[.1 5];
            pc12_range=1;
            deltaT = 1;
            tint = round(dB([1 end],1));
        case {'pc35'}
            freq_int=[.002 .1];
            pc35_range=1;
            deltaT = 60;
            tint = round(dB([1 end],1)/60)*60;
        otherwise
            error('Must choose either pc12 or pc35.');
    end
    outSampling = 1/deltaT;
    outTime = (tint(1):deltaT:tint(2))' + deltaT/2; outTime(end) = [];
else
    if freq_int(2)<freq_int(1)
        error('FREQ_INT must be [f_min f_max], f_min<f_max')
    end
    outSampling = freq_int(2)/5;
    deltaT = 1/outSampling;
    outTime = (dB(1,1):deltaT:dB(end,1))' + deltaT/2; outTime(end) = [];
end
if wantEE % Check the sampling rate
    if isempty(e)
        error('E cannot be empty for the chosen output parameters')
    end
    sampl_e=1/(e(2,1)-e(1,1));
    sampl_b=1/(dB(2,1)-dB(1,1));
    if flag_no_resamp
        if sampl_e ~= sampl_b
            error('E and B must have the same sampling for NORESAMP')
        elseif size(e,1)~=size(dB,1)
            error('E and B must have the same number of points for NORESAMP')
        end
        inSampling=sampl_e;
    else
        if     sampl_b > 1.5*sampl_e, e=irf_resamp(e,dB); B0=irf_resamp(B0,dB);...
                inSampling=sampl_b; disp('irf_pl_ebs: interpolating e to b');
        elseif sampl_e > 1.5*sampl_b, dB=irf_resamp(dB,e); B0=irf_resamp(B0,e);...
                inSampling=sampl_e; disp('irf_pl_ebs: interpolating b to e');
        elseif sampl_e == sampl_b && size(e,1)==size(dB,1),   inSampling=sampl_e;
        else   inSampling=2*sampl_e;
            t=max(e(1,1),dB(1,1)):1/inSampling:min(e(end,1),dB(end,1)); t=t';
            e=irf_resamp(e,t); dB=irf_resamp(dB,t); B0=irf_resamp(B0,t);
            irf_log('proc','interpolating b and e to 2x e sampling');
        end
        disp(['Fs=' num2str(inSampling) ', Fs_e=' num2str(sampl_e)...
            ', Fs_b=' num2str(sampl_b)]);
    end
else
    inSampling=1/(dB(2,1)-dB(1,1));
    e = [];
end
if inSampling/2<freq_int(2)
    error('F_MAX must be lower than the Nyquist frequecy')
end
if wantEE && size(e,2) <4 && flag_dEdotB0==0
    error('E must have all 3 components or flag ''dEdotdB=0'' must be given')
end
    
% Remove the last sample if the total number of samples is odd
if size(dB,1)/2 ~= floor(size(dB,1)/2)
    e=e(1:end-1,:);
    dB=dB(1:end-1,:);
    B0=B0(1:end-1,:);
    xyz=xyz(1:end-1,:);
    Bx = Bx(1:end-1,:); By = By(1:end-1,:); Bz = Bz(1:end-1,:);
end
inTime = dB(:,1);
  
% If E has all three components, transform E and B waveforms to a 
%  magnetic field aligned coordinate (FAC) and save eISR for computation 
%  of ESUM. Ohterwise we compute Ez within the main loop and do the 
%  transformation to FAC there.
if flag_want_fac
     xyz = irf_resamp(xyz,dB);
    if ~flag_dEdotB0
        eISR2=e(:,1:3);
        [dB,e]=irf_convert_fac(xyz,B0,dB,e);
    else dB = irf_convert_fac(xyz,B0,dB);
    end
else % Keep B direction for || Poynting flux
    bn = irf_norm(B0); bn(:,1) = [];
end

%% Find the frequencies for an FFT of all data and set important parameters
nd2=size(e,1)/2;
nyq=1/2;
freq=inSampling*(1:nd2)/(nd2)*nyq;
w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

Morlet_width=5.36;
freq_number=ceil((log10(freq_int(2)) - log10(freq_int(1)))*12); %to get proper overlap for Morlet
amin=log10(0.5*inSampling/freq_int(2));amax=log10(0.5*inSampling/freq_int(1));anumber=freq_number;
%  amin=0.01; % The highest frequency to consider is 0.5*sampl/10^amin
%  amax=2; % The lowest frequency to consider is 0.5*sampl/10^amax
%  anumber=400; % The number of frequencies
a=logspace(amin,amax,anumber);
%  a=logspace(0.01,2.4,100);
w0=inSampling/2; % The maximum frequency
%  sigma=5.36/w0; % The width of the Morlet wavelet
sigma=Morlet_width/w0; % The width of the Morlet wavelet

%% Make the FFT of all data
dB(:,1) = []; idxNanB = isnan(dB); dB(idxNanB) = 0; Swb=fft(dB,[],1); 
Swe = []; idxNanE = []; SweISR2 = []; idxNanEISR2 = [];% Needed for parfor
if wantEE
    fprintf('irf_ebsp ... calculate E and B wavelet transform ... ');
    e(:,1) = []; idxNanE = isnan(e); e(idxNanE)=0; Swe=fft(e,[],1);
    if flag_want_fac && ~flag_dEdotB0
        eISR2(:,1) = []; idxNanEISR2 = isnan(e); eISR2(idxNanEISR2)=0;
        SweISR2=fft(eISR2,[],1);
    end
else
    fprintf('irf_ebsp ... calculate B wavelet transform ....');
end

%% Loop through all frequencies
ndata = length(inTime); nfreq = length(a); ndataOut=length(outTime);
powerEx_plot = zeros(ndata,nfreq);
powerEy_plot = zeros(ndata,nfreq);
powerEz_plot = zeros(ndata,nfreq);
power2E_plot = zeros(ndata,nfreq);
power2E_ISR2_plot = zeros(ndata,nfreq);
power2B_plot = zeros(ndata,nfreq);
powerBx_plot = zeros(ndata,nfreq);
powerBy_plot = zeros(ndata,nfreq);
powerBz_plot = zeros(ndata,nfreq);
S_plot_x = zeros(ndata,nfreq);
S_plot_y = zeros(ndata,nfreq);
Spar_plot_z = zeros(ndata,nfreq);
polarizationEllipseRatio = zeros(ndataOut,nfreq);
polarizationSign = zeros(ndataOut,nfreq);
degreeOfPolarization = zeros(ndataOut,nfreq);
thetaSVD_fac = zeros(ndataOut,nfreq);
phiSVD_fac = zeros(ndataOut,nfreq);

% Get the correct frequencies for the wavelet transform
frequencyVec=w0./a;
censur = floor(2*a*outSampling/inSampling*nWavePeriodToAverage);
for ind_a=1:length(a), % Main loop over frequencies
  %disp([num2str(ind_a) '. frequency, ' num2str(newfreq(ind_a)) ' Hz.']);
  
  %% resample to 1 second sampling for Pc1-2 or 1 minute sampling for Pc3-5
  % average top frequencies to 1 second/1 minute
  % below will be an average over 4 wave periods. first find where one
  % sample is less than four wave periods
  if frequencyVec(ind_a)/nWavePeriodToAverage > outSampling
      avWindow = 1/outSampling;
  else
      avWindow = nWavePeriodToAverage/frequencyVec(ind_a);
  end
 
  %% Get the wavelet transform by IFFT of the FFT
  mWexp = exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  mWexp2 = repmat(mWexp,1,2); mWexp = repmat(mWexp,1,3);
  Wb = ifft(sqrt(1).*Swb.*mWexp,[],1); Wb(idxNanB) = NaN;
  We = []; WeISR2 = [];
  if wantEE
      We = ifft(sqrt(1).*Swe.*mWexp,[],1); We(idxNanE) = NaN;
      if flag_want_fac && ~flag_dEdotB0
          WeISR2 = ifft(sqrt(1).*SweISR2.*mWexp2,[],1);
          WeISR2(idxNanEISR2) = NaN;
      end
  end
  
  newfreqmat=w0/a(ind_a);
  %% Power spectrum of E and Poynting flux
  if wantEE
      % Power spectrum of E
      if flag_want_fac && ~flag_dEdotB0
          % power = (2*pi)*conj(W).*W./newfreqmat;
          SUMpowerEISR2 = sum( 2*pi*(WeISR2.*conj(WeISR2))./newfreqmat ,2);
      else SUMpowerEISR2 = sum( 2*pi*(We.*conj(We))./newfreqmat ,2);
      end
      power2E_ISR2_plot(:,ind_a) = SUMpowerEISR2;
      
      if flag_dEdotB0 % Compute Ez from dE * dB = 0
          rWe = real(We); iWe = imag(We);
          wEz = -(rWe(:,1).*Bx+rWe(:,2).*By)./Bz-...
              1j*(iWe(:,1).*Bx+iWe(:,2).*By)./Bz;
          wEz(idxBparSpinPlane,3) = NaN;
          if flag_want_fac, We = irf_convert_fac(xyz,B0,[We(:,1:2); wEz]); end
      end
      powerE = 2*pi*(We.*conj(We))./newfreqmat;
      powerE(:,4) = sum(powerE,2);
      powerEx_plot(:,ind_a) = powerE(:,1);
      powerEy_plot(:,ind_a) = powerE(:,2);
      powerEz_plot(:,ind_a) = powerE(:,3);
      power2E_plot(:,ind_a) = powerE(:,4);
      
      % Poynting flux calculations, assume E and b units mV/m and nT, get  S in uW/m^2
      coef_poynt=10/4/pi*(1/4)*(4*pi); % 4pi from wavelets, see A. Tjulins power estimates a few lines above
      S = zeros(ndata,3);
      Wex=We(:,1);Wey=We(:,2);Wez=We(:,3);
      Wbx=Wb(:,1);Wby=Wb(:,2);Wbz=Wb(:,3);
      S(:,1)= coef_poynt*real(Wey.*conj(Wbz)+conj(Wey).*Wbz-Wez.*conj(Wby)-conj(Wez).*Wby)./newfreqmat;
      S(:,2)= coef_poynt*real(Wez.*conj(Wbx)+conj(Wez).*Wbx-Wex.*conj(Wbz)-conj(Wex).*Wbz)./newfreqmat;
      S(:,3)= coef_poynt*real(Wex.*conj(Wby)+conj(Wex).*Wby-Wey.*conj(Wbx)-conj(Wey).*Wbx)./newfreqmat;
      %For some reason this code works 20% slower than the above one
      %conjWe = conj(We); conjWb = conj(Wb);
      %S = coef_poynt*real( ...
      %    + We(:,[2 3 1]).*conjWb(:,[3 1 2]) + conjWe(:,[2 3 1]).*Wb(:,[3 1 2])...
      %    - We(:,[3 1 2]).*conjWb(:,[2 3 1]) - conjWe(:,[3 1 2]).*Wb(:,[2 3 1])...
      %    )./newfreqmat;
      if flag_want_fac, Sx=S(:,1); Sy=S(:,2); Spar=S(:,3);
      else Sx=S(:,1); Sy=S(:,2); Spar = sum(S.*bn,2);
      end
      Spar_plot_z(:,ind_a) = Spar;
      S_plot_x(:,ind_a) = Sx;
      S_plot_y(:,ind_a) = Sy;
  end
  
  %% Power spectrum of B
  powerB = 2*pi*(Wb.*conj(Wb))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
  powerBx_plot(:,ind_a) = powerB(:,1);
  powerBy_plot(:,ind_a) = powerB(:,2);
  powerBz_plot(:,ind_a) = powerB(:,3);
  power2B_plot(:,ind_a) = powerB(:,4);
  
  if wantPolarization % Polarization parameters
      %% Construct spectral matrix and average it
      SM = zeros(3,3,ndata);
      SM(1,1,:) = 2*pi*(Wb(:,1).*conj(Wb(:,1)))./newfreqmat;
      SM(1,2,:) = 2*pi*(Wb(:,1).*conj(Wb(:,2)))./newfreqmat;
      SM(1,3,:) = 2*pi*(Wb(:,1).*conj(Wb(:,3)))./newfreqmat;
      SM(2,1,:) = 2*pi*(Wb(:,2).*conj(Wb(:,1)))./newfreqmat;
      SM(2,2,:) = 2*pi*(Wb(:,2).*conj(Wb(:,2)))./newfreqmat;
      SM(2,3,:) = 2*pi*(Wb(:,2).*conj(Wb(:,3)))./newfreqmat;
      SM(3,1,:) = 2*pi*(Wb(:,3).*conj(Wb(:,1)))./newfreqmat;
      SM(3,2,:) = 2*pi*(Wb(:,3).*conj(Wb(:,2)))./newfreqmat;
      SM(3,3,:) = 2*pi*(Wb(:,3).*conj(Wb(:,3)))./newfreqmat;
      SM = permute(SM,[3,1,2]);
      
      avSM = zeros(ndataOut,3,3); % Averaged SM
      for comp=1:3
          avSM(:,:,comp) = averageData(SM(:,:,comp),...
              inTime,outTime,avWindow,1);
      end
      % Remove data possibly influenced by edge effects
      censurIdx=[1:min(censur(ind_a),length(outTime))...
          max(1,length(outTime)-censur(ind_a)):length(outTime)];
      avSM(censurIdx,:,:) = NaN;
      
      %% compute singular value decomposition
      A = zeros(6,3,ndataOut); %real matrix which is superposition of real part of spectral matrix over imaginary part
      U = zeros(6,3,ndataOut);
      W = zeros(3,3,ndataOut);
      V = zeros(3,3,ndataOut);
      %wSingularValues = zeros(3,ndata2);
      %R = zeros(3,3,ndata2); %spectral matrix in coordinate defined by V axes
      A(1:3,:,:) = real(permute(avSM,[2,3,1]));
      A(4:6,:,:) = -imag(permute(avSM,[2,3,1]));   
      for i = 1:ndataOut,
          if any(any(isnan(A(:,:,i))))
              U(:,:,i) = NaN; W(:,:,i) = NaN; V(:,:,i) = NaN;
          else [U(:,:,i),W(:,:,i),V(:,:,i)] = svd(A(:,:,i),0);
          end
          %wSingularValues(:,i) = svd(A(:,:,i),0);
      end
      
      %% compute direction of propogation  
      theta=atan(sqrt(V(1,3,:).*V(1,3,:)+V(2,3,:).*V(2,3,:))./V(3,3,:));
      %phi=zeros(ndata);
      if V(1,3,:) >= 0, phi=atan(V(2,3,:)./V(1,3,:));
      elseif V(1,3,:) < 0 & V(2,3,:) < 0, phi=atan(V(2,3,:)./V(1,3,:))-pi;
      else phi=atan(V(2,3,:)./V(1,3,:))+pi;
      end
      
      %% Calculate polarization parameters    
      %planarity(:,ind_a) = 1 - sqrt(W(3,3,:)./W(1,1,:)); %planarity of polarization
      %  Lp(:,ind_a) = W(2,2,:)./W(1,1,:); %ratio of two axes of polarization ellipse
      %  planarity(:,ind_a) = 1 - sqrt(wSingularValues(3,:)./wSingularValues(1,:)); %planarity of polarization
      %    polarizationEllipseRatio(:,ind_a) = wSingularValues(2,:)./wSingularValues(1,:); %ratio of two axes of polarization ellipse
      polElliRat = W(2,2,:)./W(1,1,:); %ratio of two axes of polarization ellipse
      polElliRat(censurIdx) = NaN;
      polarizationEllipseRatio(:,ind_a) = polElliRat;
      % XXX FIXME: this loop can be optimized/eliminated
      %dop = zeros(1,ndataOut);
      %for i = 1:ndataOut,
      %    SMsqueeze = reshape(avSM(i,:,:),3,3);
      %    dop(:,i) = (3/2.*trace(real(SMsqueeze)^2)./(trace(SMsqueeze))^2 - 1/2);
      %end
      rA= real(avSM);
      dop = (3/2*(...
          rA(:,1,1).*rA(:,1,1)+rA(:,2,1).*rA(:,1,2)+rA(:,3,1).*rA(:,1,3)+...
          rA(:,1,2).*rA(:,2,1)+rA(:,2,2).*rA(:,2,2)+rA(:,3,2).*rA(:,2,3)+...
          rA(:,1,3).*rA(:,3,1)+rA(:,2,3).*rA(:,3,2)+rA(:,3,3).*rA(:,3,3))./...
          ((avSM(:,1,1)+avSM(:,2,2)+avSM(:,3,3)).^2) - 1/2); % XXX : Need a reference to this formula
      dop(censurIdx) = NaN;
      thetaSVD_fac(:,ind_a) = theta;
      phiSVD_fac(:,ind_a) = phi;
      polarizationSign(:,ind_a) = sign(imag(avSM(:,1,2))); %sign of polarization
      degreeOfPolarization(:,ind_a) = dop;
  end
end
fprintf('Done.\n');

%% set data gaps to NaN and remove edge effects
censur = floor(2*a);
for ind_a=1:length(a)
    censurIdx=[1:min(censur(ind_a),length(inTime))...
        max(1,length(inTime)-censur(ind_a)):length(inTime)];
    powerBx_plot(censurIdx,ind_a) = NaN;
    powerBx_plot(censurIdx,ind_a) = NaN;
    powerBx_plot(censurIdx,ind_a) = NaN;
    power2B_plot(censurIdx,ind_a) = NaN;
    powerEx_plot(censurIdx,ind_a) = NaN;
    powerEy_plot(censurIdx,ind_a) = NaN;
    powerEz_plot(censurIdx,ind_a) = NaN;
    power2E_plot(censurIdx,ind_a) = NaN;
    power2E_ISR2_plot(censurIdx,ind_a) = NaN;
    S_plot_x(censurIdx,ind_a) = NaN;
    S_plot_y(censurIdx,ind_a) = NaN;
    Spar_plot_z(censurIdx,ind_a) = NaN;
end
powerBx_plot = averageData(powerBx_plot,inTime,outTime);
powerBy_plot = averageData(powerBy_plot,inTime,outTime);
powerBz_plot = averageData(powerBz_plot,inTime,outTime);
power2B_plot = averageData(power2B_plot,inTime,outTime);
powerEx_plot = averageData(powerEx_plot,inTime,outTime);
powerEy_plot = averageData(powerEy_plot,inTime,outTime);
powerEz_plot = averageData(powerEz_plot,inTime,outTime);
power2E_plot = averageData(power2E_plot,inTime,outTime);
power2E_ISR2_plot = averageData(power2E_ISR2_plot,inTime,outTime);
S_plot_x = averageData(S_plot_x,inTime,outTime);
S_plot_y = averageData(S_plot_y,inTime,outTime);
Spar_plot_z = averageData(Spar_plot_z,inTime,outTime);

if wantEE,
    [S_azimuth,S_elevation,S_r]=cart2sph(S_plot_x,S_plot_y,Spar_plot_z);
    %EtoB_plot=sqrt(power2E_plot./power2B_plot);
end

if pc12_range || default_range,
  ind_lowPower = find(abs(power2B_plot) < .025);
end
if pc35_range,
  ind_lowPower = find(abs(power2B_plot) < nanmean(nanmean(log10(abs(power2B_plot)))));
  %ind_lowPower = find(abs(power2B_SM_plot) < .005);
end
ind_single_eq = find(abs(polarizationEllipseRatio) < .4);
ind2_lowPower = find(abs(degreeOfPolarization) < .4);
if wantPolarization,
    thetaSVD_fac(ind_lowPower) = NaN;
    phiSVD_fac(ind_lowPower) = NaN;
    polarizationEllipseRatio(ind_lowPower) = NaN;
    polarizationSign(ind_lowPower) = NaN;
    degreeOfPolarization(ind_lowPower) = NaN;

    thetaSVD_fac(ind2_lowPower) = NaN;
    phiSVD_fac(ind2_lowPower) = NaN;
    polarizationEllipseRatio(ind2_lowPower) = NaN;
    polarizationSign(ind2_lowPower) = NaN;
    
    %if the polarization is linear, there is a single equation and the
    %direction of the wave vector is not a unique solution
    thetaSVD_fac(ind_single_eq) = NaN;
    phiSVD_fac(ind_single_eq) = NaN;
end

%% Output
BB_xxyyzz_fac = powerBx_plot;
BB_xxyyzz_fac(:,:,2) = powerBy_plot;
BB_xxyyzz_fac(:,:,3) = powerBz_plot;
BB_xxyyzz_fac(:,:,4) = power2B_plot;
if nargout>3 % E and Poyinting Flux
    EESum_xxyy_ISR2 = power2E_ISR2_plot;
    EE_xxyyzz_FAC(:,:,4) = power2E_plot;
    EE_xxyyzz_FAC(:,:,1) = powerEx_plot;
    EE_xxyyzz_FAC(:,:,2) = powerEy_plot;
    EE_xxyyzz_FAC(:,:,3) = powerEz_plot;
    Poynting_xyz_FAC = S_plot_x;
    Poynting_xyz_FAC(:,:,2) = S_plot_y;
    Poynting_xyz_FAC(:,:,3) = Spar_plot_z;
    Poynting_rThetaPhi_FAC = S_r;
    Poynting_rThetaPhi_FAC(:,:,2) = pi/2-S_elevation;
    Poynting_rThetaPhi_FAC(:,:,3) = S_azimuth;
end
if nargout>7 % Polarization parameters
    k_thphSVD_fac = thetaSVD_fac;
    k_thphSVD_fac(:,:,2) = phiSVD_fac;
    polSVD_fac = degreeOfPolarization;
    ellipticity = polarizationEllipseRatio.*polarizationSign;
end
end

function out = averageData(data,x,y,avWindow,flagSerial)
% average data with time x to time y using window
    dtx = median(diff(x)); dty = median(diff(y));
    if nargin<4, avWindow = dty; end
    if nargin<5, flagSerial = 0; end
    dt2 = avWindow/2;
    ndataOut = length(y);
    
    % Pad data with NaNs from each side
    nPointToAdd = ceil(dt2/dtx);
    padNan = zeros(nPointToAdd,size(data,2))*NaN;
    data = [padNan; data; padNan];
    padTime = dtx*(1:nPointToAdd);
    x = [x(1)-fliplr(padTime)'; x; x(end)+padTime'];
    
    out = zeros(ndataOut,size(data,2));
    if flagSerial % Serial execution
        for i=1:length(y)
            out(i,:) = FastNanMean(data,x>=y(i)-dt2 & x<y(i)+dt2);
        end
    else % Parallel execution
        parfor i=1:length(y)
        out(i,:) = FastNanMean(data,x>=y(i)-dt2 & x<y(i)+dt2);
        end
    end
end
function m = FastNanMean(x,idx)
% Faster version of nanmean()
    xx = x(idx,:);
    % Find NaNs and set them to zero
    nans = isnan(xx); xx(nans) = 0;
    % Count up non-NaNs.
    n = sum(~nans,1);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(xx,1) ./ n;
    m(n<size(xx,1)*0.75) = NaN; % minDataFrac = .075
end

function m = nanmean(x,dim,minDataFrac)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   Revision: 1.1.8.1   Date: 2010/03/16 00:15:50 

if nargin < 3, minDataFrac = 0; end

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    if n<numel(x)*minDataFrac, m = NaN;
    else
        n(n==0) = NaN; % prevent divideByZero warnings
        % Sum up non-NaNs, and divide by the number of non-NaNs.
        m = sum(x) ./ n;
    end
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
    m(n<size(x,dim)*minDataFrac) = NaN;
end
end