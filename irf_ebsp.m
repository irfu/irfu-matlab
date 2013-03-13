function [outTime,frequencyVector,BB_xxyyzz_fac,...
    EESum_xxyy_ISR2,EE_xxyyzz_FAC,...
    Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,...
    k_thphSVD_fac,polSVD_fac,ellipticity]=...
    irf_ebsp(e,b,B0,xyz,freq_int,varargin)
%IRF_EBSP   Calculates E&B wavelet spectra, Poynting flux, polarization
%
% irf_ebsp(e,b,B0,xyz,freq_int,[OPTIONS])
% modified from irf_pl_ebs
% assumes equidistant time spacing
%
% It uses a Morlet wavelet.
% e = wave electric field, columns (t ex ey ez)
% b = wave magnetic field, columns (t bx by bz)
% B0 = background magnetic field, columns (t bx by bz)
% xyz = position vector of spacecraft, columns (t x y z)
%
% freq_int = frequency interval: either 'pc12', 'pc35' or numeric [fmin fmax]
% 
% Returns calculated parameters:
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
%    [timeVector,frequencyVector,BVector,BB_xxyyzz_fac]=...
%        irf_ebsp(e,b,B,xyz,'pc12');
%   [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
%        Poynting_xyz_FAC]=irf_ebsp(e,b,B,xyz,'pc35');
%   [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
%        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity]=...
%        irf_ebsp(e,b,B,xyz,'pc12');
%
%  See also: IRF_PL_EBS, IRF_PL_EBSP

nWavePeriodToAverage = 4; % Number of wave periods to average

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

flag_no_resamp = 0;
flag_want_fac = 0;
flag_dedotdb0 = 0;
for i=1:length(varargin)
	switch lower(varargin{i})
	case 'noresamp'
		flag_no_resamp = 1;
    case 'fac'
		flag_want_fac = 1; % Use FAC coordinate system
    case 'dedotdb=0'
		flag_dedotdb0 = 1;
	otherwise
		irf_log('fcal',['Option ''' varargin{i} '''not recognized'])
	end
end

if flag_want_fac && (isempty(B0) || isempty(xyz))
    error('B0 and XYZ must be given for option FAC')
end

pc12_range=0;
pc35_range=0;
default_range=0;
if ischar(freq_int)
    switch lower(freq_int)
        case {'pc12'}
            freq_int=[.1 5];
            pc12_range=1;
            deltaT = 1;
            tint = round(b([1 end],1));
        case {'pc35'}
            freq_int=[.002 .1];
            pc35_range=1;
            deltaT = 60;
            tint = round(b([1 end],1)/60)*60;
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
    outTime = (b(1,1):deltaT:b(end,1))' + deltaT/2; outTime(end) = [];
end
if wantEE % Check the sampling rate
    if isempty(e)
        error('E cannot be empty for the chosen output parameters')
    end
    sampl_e=1/(e(2,1)-e(1,1));
    sampl_b=1/(b(2,1)-b(1,1));
    if flag_no_resamp
        if sampl_e ~= sampl_b
            error('E and B must have the same sampling for NORESAMP')
        elseif size(e,1)~=size(b,1)
            error('E and B must have the same number of points for NORESAMP')
        end
        inSampling=sampl_e;
    else
        if     sampl_b > 1.5*sampl_e, e=irf_resamp(e,b); B0=irf_resamp(B0,b);...
                inSampling=sampl_b; disp('irf_pl_ebs: interpolating e to b');
        elseif sampl_e > 1.5*sampl_b, b=irf_resamp(b,e); B0=irf_resamp(B0,e);...
                inSampling=sampl_e; disp('irf_pl_ebs: interpolating b to e');
        elseif sampl_e == sampl_b && size(e,1)==size(b,1),   inSampling=sampl_e;
        else   inSampling=2*sampl_e;
            t=max(e(1,1),b(1,1)):1/inSampling:min(e(end,1),b(end,1)); t=t';
            e=irf_resamp(e,t); b=irf_resamp(b,t); B0=irf_resamp(B0,t);
            irf_log('proc','interpolating b and e to 2x e sampling');
        end
        disp(['Fs=' num2str(inSampling) ', Fs_e=' num2str(sampl_e)...
            ', Fs_b=' num2str(sampl_b)]);
    end
else
    inSampling=1/(b(2,1)-b(1,1));
    e = [];
end
if inSampling/2<freq_int(2)
    error('F_MAX must be lower than the Nyquist frequecy')
end
if wantEE && size(e,2) <4 && flag_dedotdb0==0
    error('E must have all 3 components or flag ''dEdotdB=0'' must be given')
end
    
%% Remove the last sample if the total number of samples is odd
if size(b,1)/2 ~= floor(size(b,1)/2)
    e=e(1:end-1,:);
    b=b(1:end-1,:);
    B0=B0(1:end-1,:);
    xyz=xyz(1:end-1,:);
end
inTime = b(:,1);
  
%% If E has all three components, transform E and B waveforms to a 
%  magnetic field aligned coordinate (FAC) and save eISR for computation 
%  of ESUM. Ohterwise we compute Ez within the main loop and do the 
%  transformation to FAC there.
B0 = irf_resamp(B0,b);
if flag_want_fac
     xyz = irf_resamp(xyz,b);
    if ~flag_dedotdb0
        eISR2=e(:,1:3);
        [b,e]=irf_convert_fac(xyz,B0,b,e);
    end
else % Keep B direction for || Poynting flux
    bn = irf_norm(B0); bn = irf_resamp(bn,outTime); bn(:,1) = [];
end
  
% set to zero NaNs
ind_nan_e = isnan(e); e(ind_nan_e)=0; ind_nan_b = isnan(b); b(ind_nan_b)=0;

%% Find the frequencies for an FFT of all data
nd2=size(e,1)/2;
nyq=1/2;
freq=inSampling*(1:nd2)/(nd2)*nyq;
w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

%% Set some important parameters
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
Swb=fft(b(:,2:4),[],1); Swe = []; SweISR2 = [];
if wantEE
    fprintf('irf_ebsp ... calculate E and B wavelet transform ... ');
    Swe=fft(e(:,2:end),[],1);
    if flag_want_fac && ~flag_dedotdb0
        SweISR2=fft(eISR2(:,2:3),[],1);
    end
else
    fprintf('irf_ebsp ... calculate B wavelet transform ....');
end

%% Get the correct frequencies for the wavelet transform
newfreq=w0./a;

%% Loop through all frequencies
ndata = size(b,1); nfreq = length(a);
 
ndataOut=length(outTime);
ind_nan_b = interp1(b(:,1),sum(ind_nan_b,2),outTime,'linear','extrap');
ind_nan_e = interp1(b(:,1),sum(ind_nan_e,2),outTime,'linear','extrap');

powerEx_plot = zeros(ndataOut,nfreq);
powerEy_plot = zeros(ndataOut,nfreq);
powerEz_plot = zeros(ndataOut,nfreq);
power2E_plot = zeros(ndataOut,nfreq);
power2E_ISR2_plot = zeros(ndataOut,nfreq);
power2B_plot = zeros(ndataOut,nfreq);
powerBx_SM_plot = zeros(ndataOut,nfreq);
powerBy_SM_plot = zeros(ndataOut,nfreq);
powerBz_SM_plot = zeros(ndataOut,nfreq);
power2B_SM_plot = zeros(ndataOut,nfreq);
polarizationEllipseRatio = zeros(ndataOut,nfreq);
polarizationSign = zeros(ndataOut,nfreq);
degreeOfPolarization = zeros(ndataOut,nfreq);
Spar_plot_z = zeros(ndataOut,nfreq);
S_plot_x = zeros(ndataOut,nfreq);
S_plot_y = zeros(ndataOut,nfreq);
thetaSVD_fac = zeros(ndataOut,nfreq);
phiSVD_fac = zeros(ndataOut,nfreq);
censur = floor(2*a*outSampling/inSampling*nWavePeriodToAverage);
censurPower = floor(2*a*outSampling/inSampling*2);
for ind_a=1:length(a), % Main loop over frequencies
  %disp([num2str(ind_a) '. frequency, ' num2str(newfreq(ind_a)) ' Hz.']);
  
  %% resample to 1 second sampling for Pc1-2 or 1 minute sampling for Pc3-5
  % average top frequencies to 1 second/1 minute
  % below will be an average over 4 wave periods. first find where one
  % sample is less than four wave periods
  if newfreq(ind_a)/nWavePeriodToAverage > outSampling
      avWindow = 1/outSampling;
  else
      avWindow = nWavePeriodToAverage/newfreq(ind_a);
  end
  
  % Remove data possibly influenced by edge effects
  censurIdx=[1:min(censur(ind_a),length(outTime))...
      max(1,length(outTime)-censur(ind_a)):length(outTime)];
  censurPwrIdx=[1:min(censurPower(ind_a),length(outTime))...
      max(1,length(outTime)-censurPower(ind_a)):length(outTime)];
 
  %% Get the wavelet transform by IFFT of the FFT
  mWexp = exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  mWexp2 = repmat(mWexp,1,2); mWexp = repmat(mWexp,1,3);
  Wb = ifft(sqrt(1).*Swb.*mWexp,[],1); We = []; WeISR2 = []; 
  if wantEE
      Wwe = sqrt(1).*Swe.*mWexp;
      We = ifft(Wwe,[],1);
      if flag_want_fac && ~flag_dedotdb0
          WeISR2 = ifft(sqrt(1).*SweISR2.*mWexp2,[],1);
      end
  end
  
  newfreqmat=w0/a(ind_a);
  Sx = []; Sy = []; Spar = [];
  if wantEE
      %% Power spectrum of E
      if flag_want_fac && ~flag_dedotdb0
          % power = (2*pi)*conj(W).*W./newfreqmat;
          SUMpowerEISR2 = sum( 2*pi*(WeISR2.*conj(WeISR2))./newfreqmat ,2);
      else SUMpowerEISR2 = sum( 2*pi*(We.*conj(We))./newfreqmat ,2);
      end
      SUMpowerEISR2 = irf_resamp([inTime SUMpowerEISR2],outTime);
      SUMpowerEISR2(:,1) = [];
      SUMpowerEISR2(censurPwrIdx) = NaN;
      if flag_dedotdb0
          % Compute Ez from dE * dB = 0
          We(:,3) = -(We(:,1).*Wb(:,1)+We(:,2).*Wb(:,2))./Wb(:,3);
          if flag_want_fac
              [Wb,We] = irf_convert_fac(xyz,B0,Wb,We);
          end
      end
      powerE = 2*pi*(We.*conj(We))./newfreqmat;
      powerE(:,4) = sum(powerE,2);
      % XXX: Maybe we want always to resample power to outTime?
      powerE = irf_resamp([inTime powerE],outTime);
      powerE(:,1) = [];
      powerE(censurPwrIdx) = NaN;
      
      %% Poynting flux calculations, assume E and b units mV/m and nT, get  S in uW/m^2
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
      % XXX: Maybe we want always to resample power to outTime?
      S = irf_resamp([inTime S],outTime); S(:,1) = [];
      S(censurPwrIdx,:) = NaN;
      if flag_want_fac, Sx=S(:,1); Sy=S(:,2); Spar=S(:,3);
      else Sx=S(:,1); Sy=S(:,2); Spar = sum(S.*bn,2);
      end 
  end
  
  %% Power spectrum of B
  powerB = 2*pi*(Wb.*conj(Wb))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
  % XXX: Maybe we want always to resample power to outTime?
  powerB = irf_resamp([inTime powerB],outTime); powerB(:,1) = [];
  powerB(censurPwrIdx,:) = NaN;
  
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
      
      avSM = zeros(ndataOut,4,3); % Averaged SM
      for comp=1:3
          avSM(:,:,comp) = irf_resamp([inTime squeeze(SM(:,:,comp))],...
              outTime,'window',avWindow);
      end
      avSM(:,1,:) = []; % get rid of the time columns
      
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
          [U(:,:,i),W(:,:,i),V(:,:,i)] = svd(A(:,:,i),0);
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
      dop = zeros(1,ndataOut);
      for i = 1:ndataOut,
          SMsqueeze = squeeze(avSM(i,:,:));
          dop(:,i) = (3/2.*trace(real(SMsqueeze)^2)./(trace(SMsqueeze))^2 - 1/2);
      end
      dop(censurIdx) = NaN;
      thetaSVD_fac(:,ind_a) = theta;
      phiSVD_fac(:,ind_a) = phi;
      polarizationSign(:,ind_a) = sign(imag(avSM(:,1,2))); %sign of polarization
      degreeOfPolarization(:,ind_a) = dop;
  end
  
  %% save power
  if wantEE,
      powerEx_plot(:,ind_a) = powerE(:,1);
      powerEy_plot(:,ind_a) = powerE(:,2);
      powerEz_plot(:,ind_a) = powerE(:,3);
      power2E_plot(:,ind_a) = powerE(:,4);
      power2E_ISR2_plot(:,ind_a) = SUMpowerEISR2;
      Spar_plot_z(:,ind_a) = Spar;
      S_plot_x(:,ind_a)=Sx;
      S_plot_y(:,ind_a)=Sy;
  end
  powerBx_SM_plot(:,ind_a) = powerB(:,1);
  powerBy_SM_plot(:,ind_a) = powerB(:,2);
  powerBz_SM_plot(:,ind_a) = powerB(:,3);
  power2B_SM_plot(:,ind_a) = powerB(:,4);
end
fprintf('Done.\n');


%% set data gaps to NaN and remove edge effects
idx_nan_e = sum(ind_nan_e,2)>0;
idx_nan_b = sum(ind_nan_b,2)>0;
powerBx_SM_plot(idx_nan_b,:) = NaN;
powerBy_SM_plot(idx_nan_b,:) = NaN;
powerBz_SM_plot(idx_nan_b,:) = NaN;
power2B_SM_plot(idx_nan_b,:) = NaN;
thetaSVD_fac(idx_nan_b,:) = NaN;
phiSVD_fac(idx_nan_b,:) = NaN;
polarizationSign(idx_nan_b,:) = NaN;
degreeOfPolarization(idx_nan_b,:) = NaN;

ndataOut=size(power2B_SM_plot,1);
if pc12_range || default_range,
  censur3=floor(.4*a);
end
if pc35_range,
  censur3=floor(.005*a);
end
censur2=floor(.4*a);

for i=1:length(idx_nan_b)-1,
    if idx_nan_b(i) < idx_nan_b(i+1),
        display('front edge');
        for j=1:length(a),
            censur_index_front=[max(i-censur3(j),1):i];
            power2B_SM_plot(censur_index_front,j) = NaN;
            thetaSVD_fac(censur_index_front,j) = NaN;
            phiSVD_fac(censur_index_front,j) = NaN;
            polarizationSign(censur_index_front,j) = NaN;
            degreeOfPolarization(censur_index_front,j) = NaN;
        end
    end
    if idx_nan_b(i) > idx_nan_b(i+1),
        display('back edge');
        for j=1:length(a),
            censur_index_back=[i:min(i+censur3(j),ndataOut)];
            power2B_SM_plot(censur_index_back,j) = NaN;
            thetaSVD_fac(censur_index_back,j) = NaN;
            phiSVD_fac(censur_index_back,j) = NaN;
            polarizationSign(censur_index_back,j) = NaN;
            degreeOfPolarization(censur_index_back,j) = NaN;
        end
    end

end

% power2B_SM_plot(censur_index_front,:) = NaN;
% thetaSVD_fac(censur_index_front,:) = NaN;
% phiSVD_fac(censur_index_front,:) = NaN;
% polarizationSign(censur_index_front,:) = NaN;
% degreeOfPolarization(censur_index_front,:) = NaN;
% power2B_SM_plot(censur_index_front,:) = NaN;
% thetaSVD_fac(censur_index_back,:) = NaN;
% phiSVD_fac(censur_index_back,:) = NaN;
% polarizationSign(censur_index_back,:) = NaN;
% degreeOfPolarization(censur_index_back,:) = NaN;


if wantEE,
    powerEx_plot(idx_nan_e,:) = NaN;
    powerEy_plot(idx_nan_e,:) = NaN;
    powerEz_plot(idx_nan_e,:) = NaN;
    power2E_plot(idx_nan_e,:) = NaN;
    power2E_ISR2_plot(idx_nan_e,:) = NaN;
    Spar_plot_z(idx_nan_e,:) = NaN;
    S_plot_x(idx_nan_e,:) = NaN;
    S_plot_y(idx_nan_e,:) = NaN;
    [S_azimuth,S_elevation,S_r]=cart2sph(S_plot_x,S_plot_y,Spar_plot_z);
    EtoB_plot=sqrt(power2E_plot./power2B_plot);
end

%ind_lowpower = find(abs(power2B_plot) < .07); %for matching Pickett's 30 march 2002
%ind_lowPower = find(abs(power2B_plot) < .05);
if pc12_range || default_range,
  ind_lowPower = find(abs(power2B_SM_plot) < .025);
end
if pc35_range,
  ind_lowPower = find(abs(power2B_SM_plot) < nanmean(nanmean(log10(abs(power2B_SM_plot)))));
  %ind_lowPower = find(abs(power2B_SM_plot) < .005);
end
ind_single_eq = find(abs(polarizationEllipseRatio) < .4);
ind2_lowPower = find(abs(degreeOfPolarization) < .4);
if wantPolarization,
    thetaSVD_fac(ind_lowPower) = NaN;
    phiSVD_fac(ind_lowPower) = NaN;
    polarizationEllipseRatio(ind_lowPower) = NaN;
    polarizationSign(ind_lowPower) = NaN;
    %Ls4(ind_lowpower) = NaN;
    degreeOfPolarization(ind_lowPower) = NaN;

    thetaSVD_fac(ind2_lowPower) = NaN;
    phiSVD_fac(ind2_lowPower) = NaN;
    polarizationEllipseRatio(ind2_lowPower) = NaN;
    polarizationSign(ind2_lowPower) = NaN;
    %Ls4(ind2_lowpower) = NaN;
    %degreeOfPolarization(ind2_lowPower) = NaN;
    
    %if the polarization is linear, there is a single equation and the
    %direction of the wave vector is not a unique solution
    thetaSVD_fac(ind_single_eq) = NaN;
    phiSVD_fac(ind_single_eq) = NaN;
    
end

if pc35_range,
    outTime=fix(outTime/60);
    outTime=outTime*60;
end

if nargout==3,
	%timeVector = e(:,1);
    outTime = fix(outTime);
	frequencyVector = newfreq;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
elseif nargout==7,
	%timeVector = e(:,1);
    outTime = fix(outTime);
	frequencyVector = newfreq;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
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
else
 	%timeVector = e(:,1);
    outTime = fix(outTime);
	frequencyVector = newfreq;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
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
    k_thphSVD_fac = thetaSVD_fac;
    k_thphSVD_fac(:,:,2) = phiSVD_fac;
    polSVD_fac = degreeOfPolarization;
    ellipticity = polarizationEllipseRatio.*polarizationSign;
end

end

function m = nanmean(x,dim)
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

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end
end