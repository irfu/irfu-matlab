function [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
    Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity]=irf_ebsp(e,b,B,xyz,varargin)
%IRF_EBSP   Calculates E wavelet spectra, B wavelet spectra, Poynting flux,
% E/B, ellipticity, polarization; 
% Returns values to be used with irf_pl_ebsp!!
%
% irf_ebsp(e,b,pos,freq_int)
% modified from irf_pl_ebs
% assumes equidistant time spacing
%
% It uses a Morlet wavelet.
% e = wave electric field, columns (t ex ey ez)
% b = wave magnetic field, columns (t bx by bz)
% B = background magnetic field, columns (t bx by bz)
% xyz = position vector of spacecraft, columns (t x y z)
%
% frequency interval chosen as input to varargin, or default [.01 5]
% 
% Returns calculated parameters
%
% Examples:
%    [timeVector,frequencyVector,BVector,BB_xxyyzz_fac]=...
%        irf_ebsp(e,b,B,xyz,'freq',[.01,1]);
%   [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
%        Poynting_xyz_FAC]=irf_ebsp(e,b,B,xyz,'freq',[.01,1]);
%   [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
%        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity]=...
%        irf_ebsp(e,b,B,xyz,'freq',[.01,1]);
% should also be possible to use default frequency range of [.01 5] **need
% to check this**
%   [timeVector,frequencyVector,BVector,BB_xxyyzz_fac,EESum_xxyyzz_ISR2,EE_xxyyzz_FAC,...
%        Poynting_xyz_FAC,Poynting_rThetaPhi_FAC,k_thphSVD_fac,polSVD_fac,ellipticity]=...
%        irf_ebsp(e,b,B,xyz);
%
% $Id$ 

%  b0=b;
  
%   sampl_low=5;
%   t_low=b(1,1):1/sampl_low:b(end,1); t_low=t_low';
%   b_low=irf_resamp(b,t_low); %low sample frequency to avoid filter issues
%   bf=irf_filt(b_low,1/600,0,[],5);
% %  bf=irf_filt(b,1/300,0,[],5);
%   b0=b_low;
%   b0(:,2:4)=b_low(:,2:4)-bf(:,2:4);
%   B=b0;
% %  b=bf;  %use B with background field removed

% %% resample to 22.5 Hz
%   sampl_b=1/(b(2,1)-b(1,1));
%   sampl=22.5;
%   t=b(1,1):1/sampl:b(end,1); t=t'; 
%   B=irf_resamp(B,t);
%   e=irf_resamp(e,t); b=irf_resamp(b,t); disp('resampling to 22.5 Hz');
%   b(:,2:4)=b(:,2:4)-B(:,2:4);

  [ax,args,nargs] = axescheck(varargin{:});
  

%% Check the sampling rate
  sampl_e=1/(e(2,1)-e(1,1));
  sampl_b=1/(b(2,1)-b(1,1));
  if     sampl_b > 1.5*sampl_e, e=irf_resamp(e,b); B=irf_resamp(B,b);...
          sampl=sampl_b; disp('irf_pl_ebs: interpolating e to b');
  elseif sampl_e > 1.5*sampl_b, b=irf_resamp(b,e); B=irf_resamp(B,e);...
          sampl=sampl_e; disp('irf_pl_ebs: interpolating b to e');
  elseif sampl_e == sampl_b & size(e)==size(b),   sampl=sampl_e;
  else   sampl=2*sampl_e; 
      t=max(e(1,1),b(1,1)):1/sampl:min(e(end,1),b(end,1)); t=t'; 
      e=irf_resamp(e,t); b=irf_resamp(b,t); B=irf_resamp(B,t);
      irf_log('proc','interpolating b and e to 2x e sampling');
  end
  disp(['Fs=' num2str(sampl) 'Fs_e=' num2str(sampl_e) 'Fs_b=' num2str(sampl_b)]);

  %% Remove the last sample if the total number of samples is odd

  if size(e,1)/2 ~= floor(size(e,1)/2)
    e=e(1:end-1,:);
    b=b(1:end-1,:);
    B=B(1:end-1,:);
    xyz=xyz(1:end-1,:);
  end

  eISR2=e; %keep ISR2 values in order to calculate EESUM_xxyy_ISR2
  
  %% transform to a magnetic field aligned coordinate (fac) 
  xyz=irf_resamp(xyz,b);
  [b,e]=irf_convert_fac(xyz,B,b,e);

  % set to zero NaNs
  ind_nan_e=isnan(e); e(ind_nan_e)=0;
  ind_nan_eISR2=isnan(eISR2); eISR2(ind_nan_eISR2)=0;
  ind_nan_b=isnan(b); b(ind_nan_b)=0;
  ind_nan_B=isnan(B); B(ind_nan_B)=0;


%% the direction of background magnetic field
%bn=irf_norm(irf_resamp(B,e));
bn=irf_norm(B);
bn(:,1) = [];
t=e(:,1);
%and magnitude
Btot=B(:,1);
Btot(:,2)=sqrt(B(:,2).*B(:,2)+B(:,3).*B(:,3)+B(:,4).*B(:,4));

  %% Find the frequencies for an FFT of all data

  nd2=size(e,1)/2;
  nyq=1/2;
  freq=sampl*(1:nd2)/(nd2)*nyq;
  w=[0,freq,-freq(end-1:-1:1)];% The frequencies corresponding to FFT

  %% Set some important parameters
  freq_int=[.01 5];
  freq_number=25;
  Morlet_width=5.36;

  if isnumeric(args{end})
      freq_int=args{end};
      args=args(1:end-2);
  elseif ismember({'freq'},args)
      disp('frequency interval values missing. using default')
      args=args(1:end-1);
  end
    
  amin=log10(0.5*sampl/freq_int(2));amax=log10(0.5*sampl/freq_int(1));anumber=freq_number;
%  amin=0.01; % The highest frequency to consider is 0.5*sampl/10^amin
%  amax=2; % The lowest frequency to consider is 0.5*sampl/10^amax
%  anumber=400; % The number of frequencies
  a=logspace(amin,amax,anumber);
%  a=logspace(0.01,2.4,100);
w0=sampl/2; % The maximum frequency
%  sigma=5.36/w0; % The width of the Morlet wavelet
sigma=Morlet_width/w0; % The width of the Morlet wavelet

disp('irf_ebsp ... calculate e and b wavelet transform ....');

%% Make the FFT of all data
Swe=fft(e(:,2:4),[],1);
SweISR2=fft(eISR2(:,2:3),[],1);
Swb=fft(b(:,2:4),[],1);

%% Get the correct frequencies for the wavelet transform
newfreq=w0./a;

wantBB = 1;
wantEE = 0;
wantPolarization = 0;

if nargout==4,
    wantBB = 1;
elseif nargout==8,
    wantEE = 1;
elseif nargout==11,
    wantEE = 1;
    wantPolarization = 1;
else
	error('irf_ebsp: unknown number of output parameters');
end

%% Loop through all frequencies
ndata = size(e,1); nfreq = length(a);
powerEx_plot = zeros(ndata,nfreq);
powerEy_plot = zeros(ndata,nfreq);
powerEz_plot = zeros(ndata,nfreq);
power2E_plot = zeros(ndata,nfreq);
power2E_ISR2_plot = zeros(ndata,nfreq);
powerBx_plot = zeros(ndata,nfreq);
powerBy_plot = zeros(ndata,nfreq);
powerBz_plot = zeros(ndata,nfreq);
power2B_plot = zeros(ndata,nfreq);
powerBx_SM_plot = zeros(ndata,nfreq);
powerBy_SM_plot = zeros(ndata,nfreq);
powerBz_SM_plot = zeros(ndata,nfreq);
power2B_SM_plot = zeros(ndata,nfreq);
polarizationEllipseRatio = zeros(ndata,nfreq);
%Ls4 = zeros(ndata,nfreq);
polarizationSign = zeros(ndata,nfreq);
degreeOfPolarization = zeros(ndata,nfreq);
Spar_plot_z = zeros(ndata,nfreq);
S_plot_x = zeros(ndata,nfreq);
S_plot_y = zeros(ndata,nfreq);
thetaSVD_fac = zeros(ndata,nfreq);
phiSVD_fac = zeros(ndata,nfreq);
parfor ind_a=1:length(a),
 % if debug, disp([num2str(ind_a) '. frequency, ' num2str(newfreq(ind_a)) ' Hz.']);end
  mWexp = exp(-sigma*sigma*((a(ind_a).*w'-w0).^2)/2);
  mWexp2 = repmat(mWexp,1,2);
  mWexp = repmat(mWexp,1,3);
  
  if wantEE, 
      Wwe = sqrt(1).*Swe.*mWexp;
      WweISR2 = sqrt(1).*SweISR2.*mWexp2;
  end
  Wwb = sqrt(1).*Swb.*mWexp;
  
  %% Get the wavelet transform by IFFT of the FFT
  if wantEE, 
      We = ifft(Wwe,[],1);
      WeISR2 = ifft(WweISR2,[],1);
  end
  Wb = ifft(Wwb,[],1);
  
  %% Calculate the power spectrum
  newfreqmat=w0/a(ind_a);
  %  power=(2*pi)*conj(W).*W./newfreqmat;
  if wantEE,
    powerE = 2*pi*(We.*conj(We))./newfreqmat;
    powerE(:,4) = sum(powerE,2);
    powerEISR2 = 2*pi*(WeISR2.*conj(WeISR2))./newfreqmat;
    powerEISR2(:,3) = sum(powerEISR2,2);
  end
  powerB = 2*pi*(Wb.*conj(Wb))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
    
  %% spectral matrix
  spectralMatrix = zeros(3,3,ndata); 
  A = zeros(6,3,ndata); %real matrix which is superposition of real part of spectral matrix over imaginary part
  U = zeros(6,3,ndata);
  W = zeros(3,3,ndata);
  V = zeros(3,3,ndata);
  wSingularValues = zeros(3,ndata);
%  SMdet = zeros(ndata);
%  R = zeros(3,3,ndata); %spectral matrix in coordinate defined by V axes
  %for i = 1:size(Wb,1)
%   for i = 1:ndata,
%       SM(:,:,i) = 2*pi*(transpose(Wb(i,:))*conj(Wb(i,:)))./newfreqmat;
%       A(:,:,i) = cat(1,real(SM(:,:,i)),imag(SM(:,:,i)));
%   end

  spectralMatrix(1,1,:) = 2*pi*(Wb(:,1).*conj(Wb(:,1)))./newfreqmat;
  spectralMatrix(1,2,:) = 2*pi*(Wb(:,1).*conj(Wb(:,2)))./newfreqmat;
  spectralMatrix(1,3,:) = 2*pi*(Wb(:,1).*conj(Wb(:,3)))./newfreqmat;
  spectralMatrix(2,1,:) = 2*pi*(Wb(:,2).*conj(Wb(:,1)))./newfreqmat;
  spectralMatrix(2,2,:) = 2*pi*(Wb(:,2).*conj(Wb(:,2)))./newfreqmat;
  spectralMatrix(2,3,:) = 2*pi*(Wb(:,2).*conj(Wb(:,3)))./newfreqmat;
  spectralMatrix(3,1,:) = 2*pi*(Wb(:,3).*conj(Wb(:,1)))./newfreqmat;
  spectralMatrix(3,2,:) = 2*pi*(Wb(:,3).*conj(Wb(:,2)))./newfreqmat;
  spectralMatrix(3,3,:) = 2*pi*(Wb(:,3).*conj(Wb(:,3)))./newfreqmat;
   
  SMpermute = permute(spectralMatrix,[3,1,2]);

  %% average spectral matrix over one wave period

  sampl_av = fix(sampl/a(ind_a));
  if sampl_av/2 == floor(sampl_av/2)
    sampl_av=sampl_av+1;
  end

  for i = sampl_av:ndata-sampl_av,
      if sampl_av < 2
          SMpermute(i,:,:)=SMpermute(i,:,:);
      elseif sampl_av > 2 && sampl_av < 4
          SMpermute(i,:,:) = (SMpermute(i-1,:,:)+SMpermute(i,:,:)+...
              SMpermute(i+1,:,:))./sampl_av;
      else 
          SMpermute(i,:,:) = sum(SMpermute(i-((sampl_av-1)/2):i+...
              ((sampl_av-1)/2),:,:),1)./sampl_av;
      end
  end            

    %% compute singular value decomposition

  A(1:3,:,:) = real(permute(SMpermute,[2,3,1]));
  A(4:6,:,:) = -imag(permute(SMpermute,[2,3,1]));

  %for i = 1:size(Wb,1)
  for i = 1:ndata,
     [U(:,:,i),W(:,:,i),V(:,:,i)] = svd(A(:,:,i),0);
%     wSingularValues(:,i) = svd(A(:,:,i),0);
  %   SMdet(i) = det(SM(:,:,i))
  end
  censur=floor(2*a);
  censur_indexes=[1:min(censur(ind_a),size(e,1)) max(1,size(e,1)-censur(ind_a)):size(e,1)];
  W(:,:,censur_indexes) = NaN;
  planarity(:,ind_a) = 1 - sqrt(W(3,3,:)./W(1,1,:)); %planarity of polarization
%  Lp(:,ind_a) = W(2,2,:)./W(1,1,:); %ratio of two axes of polarization ellipse
%  planarity(:,ind_a) = 1 - sqrt(wSingularValues(3,:)./wSingularValues(1,:)); %planarity of polarization
  if wantPolarization,
%    polarizationEllipseRatio(:,ind_a) = wSingularValues(2,:)./wSingularValues(1,:); %ratio of two axes of polarization ellipse
    polarizationEllipseRatio(:,ind_a) = W(2,2,:)./W(1,1,:); %ratio of two axes of polarization ellipse
  end
  
  %% compute direction of propogation
  if wantPolarization,
      theta=atan(sqrt(V(1,3,:).*V(1,3,:)+V(2,3,:).*V(2,3,:))./V(3,3,:));
      %phi=zeros(ndata);
      if V(1,3,:) >= 0, phi=atan(V(2,3,:)./V(1,3,:));
      elseif V(1,3,:) < 0 & V(2,3,:) < 0, phi=atan(V(2,3,:)./V(1,3,:))-pi;
      else phi=atan(V(2,3,:)./V(1,3,:))+pi;
      end
  end
  
  %% Poynting flux calculations, assume E and b units mV/m and nT, get  S in uW/m^2
  if wantEE,
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
    %Spar=sum(S.*bn,2);
    %since e and b have been converted to a FAC system, and bn is in...
    %ISR2, we cannot keep Spar as it was originally (above)
    Sx=S(:,1);
    Sy=S(:,2);
    Spar=S(:,3);
  end                             
                             
  %% Remove data possibly influenced by edge effects
  censur=floor(2*a);
  censur_indexes=[1:min(censur(ind_a),size(e,1)) max(1,size(e,1)-censur(ind_a)):size(e,1)];
  if wantEE, 
      powerE(censur_indexes,:) = NaN;
      powerEISR2(censur_indexes,:) = NaN;
  end
  powerB(censur_indexes,:) = NaN;
  if wantEE, 
      Spar(censur_indexes) = NaN;
      Sx(censur_indexes) = NaN;
      Sy(censur_indexes) = NaN;
  end
  SMpermute(censur_indexes,:,:) = NaN;
  if wantPolarization, 
    theta(censur_indexes) = NaN;
    phi(censur_indexes) = NaN;
  end

   
  %% Calculate polarization parameters
  
%   for i = 1:ndata,
%     VV=squeeze(V(:,:,i));
%     R(:,:,i)=transpose(VV)*squeeze(SM(:,:,i))*VV;
% %    R(:,:,i)=transpose(V(:,:,i))*SM(:,:,i)*V(:,:,i);
%   end
  if wantPolarization,
    thetaSVD_fac(:,ind_a) = theta;
    phiSVD_fac(:,ind_a) = phi;
    polarizationSign(:,ind_a) = sign(imag(SMpermute(:,1,2))); %sign of polarization
      %Ls4(:,ind_a) = imag(SMpermute(:,1,2))./sqrt((imag(SMpermute(:,1,2))).^2+(imag(SMpermute(:,1,3))).^2+(imag(SMpermute(:,2,3))).^2);

     % SMdet = SM(1,1,:).*(SM(2,2,:).*SM(3,3,:)-SM(2,3,:).*SM(3,2,:))-SM(2,1,:).*(SM(1,2,:).*SM(3,3,:)-SM(1,3,:).*SM(3,2,:))+SM(1,3,:).*(SM(1,2,:).*SM(2,3,:)-SM(1,3,:).*SM(2,2,:));
     % Ls6(:,ind_a) = real(sqrt(1-(4*SMdet)./(SM(1,1,:)+SM(2,2,:)+SM(3,3,:)).^2));
    degreeOfPolarization(:,ind_a) = sqrt(3/2.*(SMpermute(:,1,1).*...
        SMpermute(:,1,1)+SMpermute(:,2,2).*SMpermute(:,2,2)+...
        SMpermute(:,3,3).*SMpermute(:,3,3))./(SMpermute(:,1,1)+...
        SMpermute(:,2,2)+SMpermute(:,3,3)).^2-1/2); %degree of polarization
  end

    
  %% power
  if wantEE,
      powerEx_plot(:,ind_a) = powerE(:,1);
      powerEy_plot(:,ind_a) = powerE(:,2);
      powerEz_plot(:,ind_a) = powerE(:,3);
      power2E_plot(:,ind_a) = powerE(:,4);
      power2E_ISR2_plot(:,ind_a) = powerEISR2(:,3);
      Spar_plot_z(:,ind_a) = Spar;
      S_plot_x(:,ind_a)=Sx;
      S_plot_y(:,ind_a)=Sy;
  end
  powerBx_plot(:,ind_a) = powerB(:,1);
  powerBy_plot(:,ind_a) = powerB(:,2);
  powerBz_plot(:,ind_a) = powerB(:,3);
  power2B_plot(:,ind_a) = powerB(:,4);
  powerBx_SM_plot(:,ind_a) = SMpermute(:,1,1);
  powerBy_SM_plot(:,ind_a) = SMpermute(:,2,2);
  powerBz_SM_plot(:,ind_a) = SMpermute(:,3,3);
  power2B_SM_plot(:,ind_a) = SMpermute(:,1,1)+SMpermute(:,2,2)+SMpermute(:,3,3);
end
idx_nan_e = sum(ind_nan_e,2)>0;
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
ind_lowPower = find(abs(power2B_plot) < .025);
if wantPolarization,
    thetaSVD_fac(ind_lowPower) = NaN;
    phiSVD_fac(ind_lowPower) = NaN;
    polarizationEllipseRatio(ind_lowPower) = NaN;
    polarizationSign(ind_lowPower) = NaN;
    %Ls4(ind_lowpower) = NaN;
    degreeOfPolarization(ind_lowPower) = NaN;

    ind2_lowPower = find(abs(degreeOfPolarization) < .5);
    thetaSVD_fac(ind2_lowPower) = NaN;
    phiSVD_fac(ind2_lowPower) = NaN;
    polarizationEllipseRatio(ind2_lowPower) = NaN;
    polarizationSign(ind2_lowPower) = NaN;
    %Ls4(ind2_lowpower) = NaN;
    %degreeOfPolarization(ind2_lowPower) = NaN;
end

if nargout==4,
	timeVector = e(:,1);
	frequencyVector = newfreq;
    BVector = Btot;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
elseif nargout==8,
	timeVector = e(:,1);
	frequencyVector = newfreq;
    BVector = Btot;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
    EESum_xxyyzz_ISR2 = power2E_ISR2_plot;
    EE_xxyyzz_FAC = power2E_plot;
    Poynting_xyz_FAC = S_plot_x;
    Poynting_xyz_FAC(:,:,2) = S_plot_y;
    Poynting_xyz_FAC(:,:,3) = Spar_plot_z;
    Poynting_rThetaPhi_FAC = S_r;
    Poynting_rThetaPhi_FAC(:,:,2) = pi/2-S_elevation;
    Poynting_rThetaPhi_FAC(:,:,3) = S_azimuth;
else
 	timeVector = e(:,1);
	frequencyVector = newfreq;
    BVector = Btot;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
    EESum_xxyyzz_ISR2 = power2E_ISR2_plot;
    EE_xxyyzz_FAC = power2E_plot;
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