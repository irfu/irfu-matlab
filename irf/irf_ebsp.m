function res = irf_ebsp(e,dB,fullB,B0,xyz,freq_int,varargin)
%IRF_EBSP   Calculates wavelet spectra, Poynting flux, polarization params
%
%  irf_ebsp(E,dB,fullB,B0,xyz,freq_int,[OPTIONS])
%
%  Calculates wavelet spectra of E&B and Poynting flux using wavelets
%  (Morlet wavelet). Also computes polarization parameters of B using SVD.
%  SVD is performed on spectral matrices computed from the time series of B
%  using wavelets and then averaged over a number of wave periods.
%
%  Input:
%
%    E        - wave electric field, columns (t ex ey ez)
%    dB       - wave magnetic field, columns (t bx by bz)
%    fullB    - high resolution background magnetic field used for E*B=0,
%               columns (t bx by bz)
%    B0       - background magnetic field used for field-aligned
%               coordinates (FAC Z), columns (t bx by bz)
%    xyz      - position vector of spacecraft used for field-aligned
%               coordinates (FAC X, Y), columns (t x y z)
%    freq_int - frequency interval: either 'pc12', 'pc35' or
%               arbitrary interval [fmin fmax]
%
%  Output is a structure containing the following fields:
%
%     t           - Time
%     f           - Frequency
%     bb          - B power spectrum (xx, yy, zz)
%     ee_ss       - E power spectrum (xx+yy spacecraft coords, e.g. ISR2)
%     ee          - E power spectrum (xx, yy, zz)
%     pf_xyz      - Poynting flux (xyz)
%     pf_rtp      - Poynting flux (r, theta, phi) [angles in degrees]
%     dop         - 3D degree of polarization
%     dop2d       - 2D degree of polarization in the polarization plane
%     planarity   - planarity of polarization
%     ellipticity - ellipticity of polarization ellipse
%     k           - k-vector (theta, phi FAC) [angles in degrees]
%
%  Options:
%   'polarization' - compute polarization parameters
%   'noresamp'     - no resampling, E and dB are given at the same timeline
%   'fac'          - use FAC coordinate system (defined by B0 and
%                    optionally xyz), otherwise no coordinate system
%                    transformation is performed
%   'dEdotB=0'     - compute dEz from dB dot B = 0, uses fullB
%   'fullB=dB'     - dB contains DC field
%   'nAv'          - number of wave periods to average (default=8)
%   'facMatrix'    - specify rotation matrix to FAC system
%   'mwidthcoef'   - specify coefficient to multiple Morlet wavelet width by. 1
%   corresponds to standard Morlet wavelet.
%   'returncomplex'- set to 1 to return the complex amplitudes of the
%   wavelet transforms of both electric and magnetic field. Default is 0.
%   'downsample'   - set to 1 to downsample the output (default).
%
%  Examples:
%
%    res = irf_ebsp(e,b,B,B0,xyz,'pc12');
%
%    res = irf_ebsp(e,b,[],B0,xyz,'polarization','pc35','fullB=dB');
%
%    res = irf_ebsp(e,b,[],B0,xyz,'pc12','fullB=dB','dEdotB=0');
%
%  See also: IRF_PL_EBSP, IRF_CONVERT_FAC

% This software was developed as part of the MAARBLE (Monitoring,
% Analyzing and Assessing Radiation Belt Energization and Loss)
% collaborative research project which has received funding from the
% European Community's Seventh Framework Programme (FP7-SPACE-2011-1)
% under grant agreement n. 284520.

% Begin temporary fix to convert TS format to older format (Must include spacecraft position)
if isa(e,'TSeries')
  ttemp = e.time.epochUnix;
  datatemp = double(e.data);
  e = [ttemp, double(datatemp)];
end
if isa(dB,'TSeries')
  ttemp = dB.time.epochUnix;
  datatemp = double(dB.data);
  dB = [ttemp, datatemp];
end
if isa(fullB,'TSeries')
  ttemp = fullB.time.epochUnix;
  datatemp = double(fullB.data);
  fullB = [ttemp, datatemp];
end
if isa(B0,'TSeries')
  ttemp = B0.time.epochUnix;
  datatemp = double(B0.data);
  B0 = [ttemp, datatemp];
end
if isa(xyz,'TSeries')
  ttemp = xyz.time.epochUnix;
  datatemp = double(xyz.data);
  xyz = [ttemp, datatemp];
end
% End of temporary fix

%% Check the input
nWavePeriodToAverage = 8; % Number of wave periods to average
angleBElevationMax = 15;  % Below which we cannot apply E*B=0
facMatrix = []; % matrix for totation to FAC
mwidthcoef = 1;

wantPolarization = 0;
if isempty(e), wantEE = 0; else, wantEE = 1; end

res = struct('t',[],'f',[],'flagFac',0,...
  'bb_xxyyzzss',[],'ee_xxyyzzss',[],'ee_ss',[],...
  'pf_xyz',[],'pf_rtp',[],...
  'dop',[],'dop2d',[],'planarity',[],'ellipticity',[],'k_tp',[],...
  'fullB',fullB,'B0',B0,'r',xyz);

flag_no_resamp = 0; flag_want_fac = 0; flag_dEdotB0 = 0; flag_fullB_dB = 0;
flag_return_complex = 0;downsample = 1;
args = varargin;
while 1
  l = 1;
  if isempty(args), break, end
  switch lower(args{1})
    case 'polarization'
      wantPolarization = 1;
    case 'mwidthcoef'
      if numel(args)>1 && isnumeric(args{2})
        mwidthcoef = args{2}; l = 2;
      else
        error('parameter ''mwidthcoef'' without parameter value')
      end
    case 'noresamp'
      flag_no_resamp = 1;
    case 'fac'
      flag_want_fac = 1; % Use FAC coordinate system
    case 'dedotb=0'
      flag_dEdotB0 = 1;
    case 'fullb=db'
      flag_fullB_dB = 1;
    case 'nav'
      if length(args)==1 || ~isnumeric(args{2}) || uint8(args{2})~=args{2}
        error('NAV requires a second integer argument')
      end
      nWavePeriodToAverage = args{2}; l = 2;
    case 'facmatrix'
      if length(args)==1 || ~isstruct(args{2}) ||...
          ~isfield(args{2},'t') || ~isfield(args{2},'rotMatrix')
        error('FACMATRIX requires a second argument struct(t,rotMatrix)')
      end
      facMatrix = args{2}; l = 2;
    case 'returncomplex'
      flag_return_complex = args{2}; l = 2;
    case 'downsample'
      downsample = args{2}; l = 2;
    otherwise
      irf_log('fcal',['Option ''' args{1} '''not recognized'])
  end
  args = args(l+1:end);
end

if flag_want_fac && isempty(facMatrix)
  if isempty(B0)
    error('irf_ebsp(): at least B0 should be given for option FAC');
  end
  if isempty(xyz)
    irf_log('fcal','assuming s/c position [1 0 0] for estimating FAC');
    xyz=[0 1 0 0];
  end
  xyz = irf_resamp(xyz,dB);
end
B0 = irf_resamp(B0,dB);
if flag_fullB_dB
  fullB = dB;
  res.fullB = fullB;
  dB(:,2:4) = dB(:,2:4) - B0(:,2:4);
end
if flag_dEdotB0 && isempty(fullB)
  error('fullB must be given for option dEdotB=0')
end

pc12_range=0; pc35_range=0; other_range=0;
if ischar(freq_int)
  switch lower(freq_int)
    case {'pc12'}
      pc12_range=1;
      freq_int=[.1 5];
      deltaT = 1;
      tint = round(dB([1 end],1));
    case {'pc35'}
      pc35_range=1;
      freq_int=[.002 .1];
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
  other_range=1;
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
    else,   inSampling=2*sampl_e;
      t=max(e(1,1),dB(1,1)):1/inSampling:min(e(end,1),dB(end,1)); t=t';
      e=irf_resamp(e,t); dB=irf_resamp(dB,t); B0=irf_resamp(B0,t);
      fullB = irf_resamp(fullB,t);
      irf_log('proc','interpolating b and e to 2x e sampling');
    end
    disp(['Fs=' num2str(inSampling) ', Fs_e=' num2str(sampl_e)...
      ', Fs_b=' num2str(sampl_b)]);
  end
else
  inSampling=1/(dB(2,1)-dB(1,1));
  e = [];
end
if inSampling/2 < freq_int(2)
  error('F_MAX must be lower than the Nyquist frequency')
end
if wantEE && size(e,2) <4 && flag_dEdotB0==0
  error('E must have all 3 components or flag ''dEdotdB=0'' must be given')
end

if size(dB,1)/2 ~= floor(size(dB,1)/2)
  dB=dB(1:end-1,:);
  B0=B0(1:end-1,:);
  if isempty(facMatrix)
    xyz=xyz(1:end-1,:);
  else
    facMatrix.t = facMatrix.t(1:end-1,:);
    facMatrix.rotMatrix = facMatrix.rotMatrix(1:end-1,:,:);
  end
  if wantEE, e=e(1:end-1,:); end
end
inTime = dB(:,1);
if ~downsample
  outTime = inTime;
end
Bx = []; By = []; Bz = []; idxBparSpinPlane = [];
if flag_dEdotB0
  Bx = fullB(:,2); By = fullB(:,3); Bz = fullB(:,4); % Needed for parfor

  % Remove the last sample if the total number of samples is odd
  if size(fullB,1)/2 ~= floor(size(fullB,1)/2)
    Bx = Bx(1:end-1,:); By = By(1:end-1,:); Bz = Bz(1:end-1,:);
  end

  angleBElevation=atand(Bz./sqrt(Bx.^2+By.^2));
  idxBparSpinPlane= abs(angleBElevation)<angleBElevationMax;
end

% If E has all three components, transform E and B waveforms to a
%  magnetic field aligned coordinate (FAC) and save eISR for computation
%  of ESUM. Ohterwise we compute Ez within the main loop and do the
%  transformation to FAC there.
timeB0 = 0;
if flag_want_fac
  res.flagFac = 1;
  timeB0 = B0(:,1);
  if wantEE
    if ~flag_dEdotB0
      eISR2=e(:,1:3);
      if size(e,2)<4
        error('E must be a 3D vector to be rotated to FAC')
      end
      if isempty(facMatrix), e=irf_convert_fac(e,B0,xyz);
      else, e=irf_convert_fac(e,facMatrix);
      end
    end
  end
  if isempty(facMatrix), dB = irf_convert_fac(dB,B0,xyz);
  else, dB=irf_convert_fac(dB,facMatrix);
  end
end

%% Find the frequencies for an FFT of all data and set important parameters
nd2=length(inTime)/2; nyq=1/2; freq=inSampling*(1:nd2)/(nd2)*nyq;
w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

Morlet_width=5.36*mwidthcoef;
freq_number=ceil((log10(freq_int(2)) - log10(freq_int(1)))*12*mwidthcoef); %to get proper overlap for Morlet
amin=log10(0.5*inSampling/freq_int(2));
amax=log10(0.5*inSampling/freq_int(1));
anumber=freq_number;
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
    eISR2(:,1) = []; idxNanEISR2 = isnan(eISR2); eISR2(idxNanEISR2)=0;
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
if flag_return_complex
  Ex_complex = zeros(ndata,nfreq);
  Ey_complex = zeros(ndata,nfreq);
  Ez_complex = zeros(ndata,nfreq);
  Bx_complex = zeros(ndata,nfreq);
  By_complex = zeros(ndata,nfreq);
  Bz_complex = zeros(ndata,nfreq);
end
S_plot_x = zeros(ndata,nfreq);
S_plot_y = zeros(ndata,nfreq);
S_plot_z = zeros(ndata,nfreq);
planarity = zeros(ndataOut,nfreq);
ellipticity = zeros(ndataOut,nfreq);
degreeOfPolarization3D = zeros(ndataOut,nfreq);
degreeOfPolarization2D = zeros(ndataOut,nfreq);
thetaSVD_fac = zeros(ndataOut,nfreq);
phiSVD_fac = zeros(ndataOut,nfreq);

% Get the correct frequencies for the wavelet transform
frequencyVec=w0./a;
censur = floor(2*a*outSampling/inSampling*nWavePeriodToAverage);
parfor ind_a=1:length(a) % Main loop over frequencies
  %disp([num2str(ind_a) '. frequency, ' num2str(newfreq(ind_a)) ' Hz.']);

  %% resample to 1 second sampling for Pc1-2 or 1 minute sampling for Pc3-5
  % average top frequencies to 1 second/1 minute
  % below will be an average over 8 wave periods. first find where one
  % sample is less than eight wave periods
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
    if size(Swe,2) == 2, We = ifft(sqrt(1).*Swe.*mWexp2,[],1);
    else, We = ifft(sqrt(1).*Swe.*mWexp,[],1);
    end
    We(idxNanE) = NaN;
    if flag_want_fac && ~flag_dEdotB0
      WeISR2 = ifft(sqrt(1).*SweISR2.*mWexp2,[],1);
      WeISR2(idxNanEISR2) = NaN;
    end
  end

  newfreqmat=w0/a(ind_a);
  %% Power spectrum of E and Poynting flux
  if wantEE
    % Power spectrum of E, power = (2*pi)*conj(W).*W./newfreqmat
    if flag_want_fac && ~flag_dEdotB0
      SUMpowerEISR2 = sum( 2*pi*(WeISR2.*conj(WeISR2))./newfreqmat ,2);
    else, SUMpowerEISR2 = sum( 2*pi*(We.*conj(We))./newfreqmat ,2);
    end
    power2E_ISR2_plot(:,ind_a) = SUMpowerEISR2;

    if flag_dEdotB0 % Compute Ez from dE * B = 0
      rWe = real(We); iWe = imag(We);
      wEz = -(rWe(:,1).*Bx+rWe(:,2).*By)./Bz-...
        1j*(iWe(:,1).*Bx+iWe(:,2).*By)./Bz;
      wEz(idxBparSpinPlane) = NaN;
      if flag_want_fac
        if isempty(facMatrix)
          We = irf_convert_fac([timeB0 We(:,1:2) wEz],B0,xyz);
        else
          We = irf_convert_fac([timeB0 We(:,1:2) wEz],facMatrix);
        end
        We(:,1) = [];
      else, We = [We(:,1:2) wEz];
      end
    end
    powerE = 2*pi*(We.*conj(We))./newfreqmat;
    powerE(:,4) = sum(powerE,2);
    powerEx_plot(:,ind_a) = powerE(:,1);
    powerEy_plot(:,ind_a) = powerE(:,2);
    powerEz_plot(:,ind_a) = powerE(:,3);
    power2E_plot(:,ind_a) = powerE(:,4);
    if flag_return_complex
      Ex_complex(:,ind_a) = sqrt(2*pi).*We(:,1)./sqrt(newfreqmat);
      Ey_complex(:,ind_a) = sqrt(2*pi).*We(:,2)./sqrt(newfreqmat);
      Ez_complex(:,ind_a) = sqrt(2*pi).*We(:,3)./sqrt(newfreqmat);
    end

    % Poynting flux calculations, assume E and b units mV/m and nT, get  S in uW/m^2
    coef_poynt=10/4/pi*(1/4)*(4*pi); % 4pi from wavelets, see A. Tjulins power estimates a few lines above
    S = zeros(ndata,3);
    Wex=We(:,1);Wey=We(:,2);Wez=We(:,3);
    Wbx=Wb(:,1);Wby=Wb(:,2);Wbz=Wb(:,3);
    S(:,1)= coef_poynt*real(Wey.*conj(Wbz)+conj(Wey).*Wbz-Wez.*conj(Wby)-conj(Wez).*Wby)./newfreqmat;
    S(:,2)= coef_poynt*real(Wez.*conj(Wbx)+conj(Wez).*Wbx-Wex.*conj(Wbz)-conj(Wex).*Wbz)./newfreqmat;
    S(:,3)= coef_poynt*real(Wex.*conj(Wby)+conj(Wex).*Wby-Wey.*conj(Wbx)-conj(Wey).*Wbx)./newfreqmat;

    S_plot_x(:,ind_a) = S(:,1);
    S_plot_y(:,ind_a) = S(:,2);
    S_plot_z(:,ind_a) = S(:,3);
  end

  %% Power spectrum of B
  powerB = 2*pi*(Wb.*conj(Wb))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
  powerBx_plot(:,ind_a) = powerB(:,1);
  powerBy_plot(:,ind_a) = powerB(:,2);
  powerBz_plot(:,ind_a) = powerB(:,3);
  power2B_plot(:,ind_a) = powerB(:,4);
  if flag_return_complex
    Bx_complex(:,ind_a) = sqrt(2*pi).*Wb(:,1)./sqrt(newfreqmat);
    By_complex(:,ind_a) = sqrt(2*pi).*Wb(:,2)./sqrt(newfreqmat);
    Bz_complex(:,ind_a) = sqrt(2*pi).*Wb(:,3)./sqrt(newfreqmat);
  end
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
      avSM(:,:,comp) = AverageData(SM(:,:,comp),...
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
    for i = 1:ndataOut
      if any(any(isnan(A(:,:,i))))
        U(:,:,i) = NaN; W(:,:,i) = NaN; V(:,:,i) = NaN;
      else, [U(:,:,i),W(:,:,i),V(:,:,i)] = svd(A(:,:,i),0);
      end
      %wSingularValues(:,i) = svd(A(:,:,i),0);
    end

    %% compute direction of propogation
    signKz = sign(V(3,3,:));
    V(3,3,:) = V(3,3,:).*signKz;
    V(2,3,:) = V(2,3,:).*signKz;
    V(1,3,:) = V(1,3,:).*signKz;
    thetaSVD_fac(:,ind_a) = ...
      abs(squeeze(atand(sqrt(V(1,3,:).*V(1,3,:)+V(2,3,:).*V(2,3,:))./V(3,3,:)))); %#ok<PFOUS>
    phiSVD_fac(:,ind_a) = squeeze(atan2d(V(2,3,:),V(1,3,:))); %#ok<PFOUS>

    %% Calculate polarization parameters
    planarityLocal = squeeze(1 - sqrt(W(3,3,:)./W(1,1,:)));
    planarityLocal(censurIdx) = NaN;
    planarity(:,ind_a) = planarityLocal;

    %ellipticity: ratio of axes of polarization ellipse axes*sign of polarization
    ellipticityLocal = ...
      squeeze(W(2,2,:)./W(1,1,:)).*sign(imag(avSM(:,1,2)));
    ellipticityLocal(censurIdx) = NaN;
    ellipticity(:,ind_a) = ellipticityLocal;

    % DOP = sqrt[(3/2.*trace(SM^2)./(trace(SM))^2 - 1/2)]; Samson, 1973, JGR
    dop = sqrt((3/2*(...
      avSM(:,1,1).*avSM(:,1,1)+avSM(:,2,1).*avSM(:,1,2)+avSM(:,3,1).*avSM(:,1,3)+...
      avSM(:,1,2).*avSM(:,2,1)+avSM(:,2,2).*avSM(:,2,2)+avSM(:,3,2).*avSM(:,2,3)+...
      avSM(:,1,3).*avSM(:,3,1)+avSM(:,2,3).*avSM(:,3,2)+avSM(:,3,3).*avSM(:,3,3))./...
      ((avSM(:,1,1)+avSM(:,2,2)+avSM(:,3,3)).^2) - 1/2));
    dop(censurIdx) = NaN;
    degreeOfPolarization3D(:,ind_a) = dop;

    % DOP in 2D = sqrt[2*trace(rA^2)/trace(rA)^2 - 1)]; Ulrich
    Vnew = permute(V,[3,1,2]);
    avSM2dim = mult_mat(Vnew,mult_mat(avSM,transpose_mat(Vnew)));
    avSM2dim = avSM2dim(:,1:2,1:2);
    avSM = avSM2dim;
    dop2dim = sqrt((2*(...
      avSM(:,1,1).*avSM(:,1,1)+avSM(:,2,1).*avSM(:,1,2)+...
      avSM(:,1,2).*avSM(:,2,1)+avSM(:,2,2).*avSM(:,2,2))./...
      ((avSM(:,1,1)+avSM(:,2,2)).^2) - 1));
    dop2dim(censurIdx) = NaN;
    degreeOfPolarization2D(:,ind_a) = dop2dim;

  end % wantPolarization
end % main parfor loop
fprintf('Done.\n');
%% set data gaps to NaN and remove edge effects
censur = floor(2*a);
for ind_a=1:length(a)
  censurIdx=[1:min(censur(ind_a),length(inTime))...
    max(1,length(inTime)-censur(ind_a)):length(inTime)];
  powerBx_plot(censurIdx,ind_a) = NaN;
  powerBy_plot(censurIdx,ind_a) = NaN;
  powerBz_plot(censurIdx,ind_a) = NaN;
  power2B_plot(censurIdx,ind_a) = NaN;
  if flag_return_complex
    Bx_complex(censurIdx,ind_a) = NaN;
    By_complex(censurIdx,ind_a) = NaN;
    Bz_complex(censurIdx,ind_a) = NaN;
  end
  if wantEE
    powerEx_plot(censurIdx,ind_a) = NaN;
    powerEy_plot(censurIdx,ind_a) = NaN;
    powerEz_plot(censurIdx,ind_a) = NaN;
    power2E_plot(censurIdx,ind_a) = NaN;
    power2E_ISR2_plot(censurIdx,ind_a) = NaN;
    if flag_return_complex
      Ex_complex(censurIdx,ind_a) = NaN;
      Ey_complex(censurIdx,ind_a) = NaN;
      Ez_complex(censurIdx,ind_a) = NaN;
    end
    S_plot_x(censurIdx,ind_a) = NaN;
    S_plot_y(censurIdx,ind_a) = NaN;
    S_plot_z(censurIdx,ind_a) = NaN;
  end
end

%% remove edge effects from data gaps
idxNanE = sum(idxNanE,2)>0;
idxNanB = sum(idxNanB,2)>0;
idxNanEISR2 = sum(idxNanEISR2,2)>0;

ndata2=size(power2B_plot,1);
if pc12_range || other_range
  censur3=floor(1.8*a);
end
if pc35_range
  censur3=floor(.4*a);
end
for i=1:length(idxNanB)-1
  if idxNanB(i) < idxNanB(i+1)
    for j=1:length(a)
      censur_index_front = max(i-censur3(j),1):i;
      powerBx_plot(censur_index_front,j) = NaN;
      powerBy_plot(censur_index_front,j) = NaN;
      powerBz_plot(censur_index_front,j) = NaN;
      power2B_plot(censur_index_front,j) = NaN;
      if flag_return_complex
        Bx_complex(censur_index_front,j) = NaN;
        By_complex(censur_index_front,j) = NaN;
        Bz_complex(censur_index_front,j) = NaN;
      end
      S_plot_x(censur_index_front,j) = NaN;
      S_plot_y(censur_index_front,j) = NaN;
      S_plot_z(censur_index_front,j) = NaN;
    end
  end
  if idxNanB(i) > idxNanB(i+1)
    for j=1:length(a)
      censur_index_back = i:min(i+censur3(j),ndata2);
      powerBx_plot(censur_index_back,j) = NaN;
      powerBy_plot(censur_index_back,j) = NaN;
      powerBz_plot(censur_index_back,j) = NaN;
      power2B_plot(censur_index_back,j) = NaN;
      if flag_return_complex
        Bx_complex(censur_index_back,j) = NaN;
        By_complex(censur_index_back,j) = NaN;
        Bz_complex(censur_index_back,j) = NaN;
      end
      S_plot_x(censur_index_back,j) = NaN;
      S_plot_y(censur_index_back,j) = NaN;
      S_plot_z(censur_index_back,j) = NaN;
    end
  end
end

ndata3=size(power2E_plot,1);
for i=1:length(idxNanE)-1
  if idxNanE(i) < idxNanE(i+1)
    for j=1:length(a)
      censur_index_front = max(i-censur3(j),1):i;
      powerEx_plot(censur_index_front,j) = NaN;
      powerEy_plot(censur_index_front,j) = NaN;
      powerEz_plot(censur_index_front,j) = NaN;
      power2E_plot(censur_index_front,j) = NaN;
      power2E_ISR2_plot(censur_index_front,j) = NaN;
      if flag_return_complex
        Ex_complex(censur_index_front,j) = NaN;
        Ey_complex(censur_index_front,j) = NaN;
        Ez_complex(censur_index_front,j) = NaN;
      end
      S_plot_x(censur_index_front,j) = NaN;
      S_plot_y(censur_index_front,j) = NaN;
      S_plot_z(censur_index_front,j) = NaN;
    end
  end
  if idxNanE(i) > idxNanE(i+1)
    for j=1:length(a)
      censur_index_back = i:min(i+censur3(j),ndata3);
      powerEx_plot(censur_index_back,j) = NaN;
      powerEy_plot(censur_index_back,j) = NaN;
      powerEz_plot(censur_index_back,j) = NaN;
      power2E_plot(censur_index_back,j) = NaN;
      power2E_ISR2_plot(censur_index_back,j) = NaN;
      if flag_return_complex
        Ex_complex(censur_index_back,j) = NaN;
        Ey_complex(censur_index_back,j) = NaN;
        Ez_complex(censur_index_back,j) = NaN;
      end
      S_plot_x(censur_index_back,j) = NaN;
      S_plot_y(censur_index_back,j) = NaN;
      S_plot_z(censur_index_back,j) = NaN;
    end
  end
end

ndata4=size(power2E_ISR2_plot,1);
for i=1:length(idxNanEISR2)-1
  if idxNanEISR2(i) < idxNanEISR2(i+1)
    for j=1:length(a)
      censur_index_front = max(i-censur3(j),1):i;
      power2E_ISR2_plot(censur_index_front,j) = NaN;
    end
  end
  if idxNanEISR2(i) > idxNanEISR2(i+1)
    for j=1:length(a)
      censur_index_back = i:min(i+censur3(j),ndata4);
      power2E_ISR2_plot(censur_index_back,j) = NaN;
    end
  end
end



%%
powerBx_plot = AverageData(powerBx_plot,inTime,outTime);
powerBy_plot = AverageData(powerBy_plot,inTime,outTime);
powerBz_plot = AverageData(powerBz_plot,inTime,outTime);
power2B_plot = AverageData(power2B_plot,inTime,outTime);
if flag_return_complex
  Bx_complex = AverageData(Bx_complex,inTime,outTime);
  By_complex = AverageData(By_complex,inTime,outTime);
  Bz_complex = AverageData(Bz_complex,inTime,outTime);
end
bb_xxyyzzss = powerBx_plot;
bb_xxyyzzss(:,:,2) = powerBy_plot;
bb_xxyyzzss(:,:,3) = powerBz_plot;
bb_xxyyzzss(:,:,4) = power2B_plot;


% Output
res.t = outTime;
res.f = frequencyVec;
res.bb_xxyyzzss = bb_xxyyzzss;
if flag_return_complex
  B_complex = Bx_complex;
  B_complex(:,:,2) = By_complex;
  B_complex(:,:,3) = Bz_complex;
  res.B_complex = B_complex;
end
if wantEE
  powerEx_plot = AverageData(powerEx_plot,inTime,outTime);
  powerEy_plot = AverageData(powerEy_plot,inTime,outTime);
  powerEz_plot = AverageData(powerEz_plot,inTime,outTime);
  power2E_plot = AverageData(power2E_plot,inTime,outTime);
  if flag_return_complex
    Ex_complex = AverageData(Ex_complex,inTime,outTime);
    Ey_complex = AverageData(Ey_complex,inTime,outTime);
    Ez_complex = AverageData(Ez_complex,inTime,outTime);
  end
  power2E_ISR2_plot = AverageData(power2E_ISR2_plot,inTime,outTime);
  S_plot_x = AverageData(S_plot_x,inTime,outTime);
  S_plot_y = AverageData(S_plot_y,inTime,outTime);
  S_plot_z = AverageData(S_plot_z,inTime,outTime);
  [S_azimuth,S_elevation,S_r]=cart2sph(S_plot_x,S_plot_y,S_plot_z);

  ee_xxyyzzss(:,:,4) = power2E_plot;
  ee_xxyyzzss(:,:,1) = powerEx_plot;
  ee_xxyyzzss(:,:,2) = powerEy_plot;
  ee_xxyyzzss(:,:,3) = powerEz_plot;
  Poynting_XYZ = S_plot_x;
  Poynting_XYZ(:,:,2) = S_plot_y;
  Poynting_XYZ(:,:,3) = S_plot_z;
  Poynting_RThPh = S_r;
  Poynting_RThPh(:,:,2) = pi/2-S_elevation;
  Poynting_RThPh(:,:,3) = S_azimuth;
  Poynting_RThPh(:,:,2:3) = Poynting_RThPh(:,:,2:3)*180/pi;

  % Output
  res.ee_ss = power2E_ISR2_plot;
  res.ee_xxyyzzss = ee_xxyyzzss;
  if flag_return_complex
    E_complex = Ex_complex;
    E_complex(:,:,2) = Ey_complex;
    E_complex(:,:,3) = Ez_complex;
    res.E_complex = E_complex;
  end
  res.pf_xyz = Poynting_XYZ;
  res.pf_rtp = Poynting_RThPh;
end

if wantPolarization
  % Define parameters for which we cannot compute the wave vector
  indLowPlanarity = planarity < 0.5;
  indLowEllipticity = abs(ellipticity) < .2;

  thetaSVD_fac(indLowPlanarity) = NaN;
  phiSVD_fac(indLowPlanarity) = NaN;

  thetaSVD_fac(indLowEllipticity) = NaN;
  phiSVD_fac(indLowEllipticity) = NaN;

  k_ThPhSVD_fac = thetaSVD_fac;
  k_ThPhSVD_fac(:,:,2) = phiSVD_fac;

  % Output
  res.dop = degreeOfPolarization3D;
  res.dop2d = degreeOfPolarization2D;
  res.planarity = planarity;
  res.ellipticity = ellipticity;
  res.k_tp = k_ThPhSVD_fac;
end
end % Main function

function out = AverageData(data,x,y,avWindow,flagSerial)
if length(x) == length(y);out = data;return;end
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
function out=transpose_mat(inp)
out = inp;
if numel(size(inp))==2 || (size(inp,2)~=size(inp,3))
  error('not impemented');
end
for ii=1:size(inp,2)
  for jj=ii+1:size(inp,2)
    out(:,ii,jj) = inp(:,jj,ii);
    out(:,jj,ii) = inp(:,ii,jj);
  end
end
end
function out = mult_mat(inp1,inp2)
dimInp1 = numel(size(inp1));
dimInp2 = numel(size(inp2));
if (dimInp1==dimInp2)
  numOfMult=size(inp1,3);
  T = zeros(size(inp1,1),size(inp1,2),size(inp2,3));
  for ii=1:size(inp1,2)
    for jj=1:size(inp2,3)
      for kk=1:numOfMult
        T(:,ii,jj)=T(:,ii,jj)+inp1(:,ii,kk).*inp2(:,kk,jj);
      end
    end
  end
elseif (dimInp1==3) && (dimInp2==2)
  numOfOutp=size(inp1,2);
  numOfInp=size(inp2,2);
  T = inp2(:,1)*zeros(1,numOfOutp);
  for ii=1:numOfOutp
    for jj=1:numOfInp
      T(:,ii)=T(:,ii)+inp1(:,ii,jj).*inp2(:,jj);
    end
  end
end
out = T;
end
