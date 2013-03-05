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
% frequency interval chosen as input to varargin, either 'pc12' or 'pc35' 
% if something else is used, default [.02 5]
% 
% Returns calculated parameters
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
% $Id$ 


  [ax,args,nargs] = axescheck(varargin{:});

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
  %freq_number=25;
  Morlet_width=5.36;
  pc12_range=0;
  pc35_range=0;
  default_range=0;
  
  switch lower(args{1})
      case {'pc12'}
          freq_range = args{1};
          freq_int=[.1 5];
          pc12_range=1;
      case {'pc35'}
          freq_range = args{1};
          freq_int=[.002 .1];
          pc35_range=1;
      otherwise 
          display('Must choose either pc12 or pc35. Using default [.01 5]');
          freq_int=[.01 5];
          default_range=1;
  end
   freq_number=ceil((log10(freq_int(2)) - log10(freq_int(1)))*12); %to get proper overlap for Morlet
         
%   if isnumeric(args{end})
%       freq_int=args{end};
%       args=args(1:end-2);
%   elseif ismember({'freq'},args)
%       disp('frequency interval values missing. using default')
%       args=args(1:end-1);
%   end
    
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
%display(newfreq);

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
if pc12_range
    sampl1 = 1;
end
if pc35_range
    sampl1 = 1/60;
end
if default_range
    sampl1 = 1;
end
t1 = e(1,1):1/sampl1:e(end,1); t1=t1'; 
ndata2=length(t1);
ind_nan_b = interp1(t,sum(ind_nan_b,2),t1,'linear','extrap');
ind_nan_e = interp1(t,sum(ind_nan_e,2),t1,'linear','extrap');

powerEx_plot = zeros(ndata2,nfreq);
powerEy_plot = zeros(ndata2,nfreq);
powerEz_plot = zeros(ndata2,nfreq);
power2E_plot = zeros(ndata2,nfreq);
power2E_ISR2_plot = zeros(ndata2,nfreq);
powerBx_plot = zeros(ndata2,nfreq);
powerBy_plot = zeros(ndata2,nfreq);
powerBz_plot = zeros(ndata2,nfreq);
power2B_plot = zeros(ndata2,nfreq);
powerBx_SM_plot = zeros(ndata2,nfreq);
powerBy_SM_plot = zeros(ndata2,nfreq);
powerBz_SM_plot = zeros(ndata2,nfreq);
power2B_SM_plot = zeros(ndata2,nfreq);
polarizationEllipseRatio = zeros(ndata2,nfreq);
%Ls4 = zeros(ndata,nfreq);
polarizationSign = zeros(ndata2,nfreq);
degreeOfPolarization = zeros(ndata2,nfreq);
Spar_plot_z = zeros(ndata2,nfreq);
S_plot_x = zeros(ndata2,nfreq);
S_plot_y = zeros(ndata2,nfreq);
thetaSVD_fac = zeros(ndata2,nfreq);
phiSVD_fac = zeros(ndata2,nfreq);
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
  
  %Wwb = interp1(t,Wwb,t1,'linear','extrap')
  Wwb2 = interp1(t,Wwb,t1,'linear','extrap')
  Wwe = interp1(t,Wwe,t1,'linear','extrap')
  WweISR2 = interp1(t,WweISR2,t1,'linear','extrap')
  
  %% Get the wavelet transform by IFFT of the FFT
  if wantEE, 
      We = ifft(Wwe,[],1);
      WeISR2 = ifft(WweISR2,[],1);
  end
  Wb = ifft(Wwb,[],1);
  Wb2 = ifft(Wwb2,[],1);
  
  %% Calculate the power spectrum
  newfreqmat=w0/a(ind_a);
  %  power=(2*pi)*conj(W).*W./newfreqmat;
  if wantEE,
    powerE = 2*pi*(We.*conj(We))./newfreqmat;
    powerE(:,4) = sum(powerE,2);
    powerEISR2 = 2*pi*(WeISR2.*conj(WeISR2))./newfreqmat;
    powerEISR2(:,3) = sum(powerEISR2,2);
  end
  powerB = 2*pi*(Wb2.*conj(Wb2))./newfreqmat;
  powerB(:,4) = sum(powerB,2);
    
  %% spectral matrix
  spectralMatrix = zeros(3,3,ndata); 
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
  
  %% resample to 1 second sampling for Pc1-2 or 1 minute sampling for Pc3-5
  %average top frequencies to 1 second/1 minute 
  %below will be an average over 4 wave periods. first find where one
  %sample is less than four wave periods
  num_wavePeriod = 4;
  sampl_av = int16(num_wavePeriod*sampl1/newfreqmat);
  %SMpermute_av=zeros(size(SMpermute));
  SMpermute_av=SMpermute;
  
  if pc35_range
      if sampl_av < 2,
          for i=1:749,
              SMpermute_av(i,:,:) = mean(SMpermute(1:i+749,:,:),1);
          end
          for i=750:ndata-749,
              SMpermute_av(i,:,:) = mean(SMpermute(i-749:i+749,:,:),1);
          end
          for i=ndata-749:ndata,
              SMpermute_av(i,:,:) = mean(SMpermute(i-749:end,:,:),1);
          end
      end
  end
  if pc12_range || default_range
      if sampl_av < 2,
          for i=1:12,
              SMpermute_av(i,:,:) = mean(SMpermute(1:i+12,:,:),1);
          end
          for i=13:ndata-12,
              SMpermute_av(i,:,:) = mean(SMpermute(i-12:i+12,:,:),1);
          end
          for i=ndata-12:ndata,
              SMpermute_av(i,:,:) = mean(SMpermute(i-12:end,:,:),1);
          end
      end
  end
      

  SMpermute = SMpermute_av;
  
  %t1 is defined outside of parfor
  SMpermute = interp1(t,SMpermute,t1,'linear','extrap')
  sampl = sampl1;
  %ndata2 = int16(size(SMpermute,1));
  ndata2 = size(SMpermute,1);
%display(ndata2); 
%display(size(SMpermute));
  A = zeros(6,3,ndata2); %real matrix which is superposition of real part of spectral matrix over imaginary part
  U = zeros(6,3,ndata2);
  W = zeros(3,3,ndata2);
  V = zeros(3,3,ndata2);
  wSingularValues = zeros(3,ndata2);
  R = zeros(3,3,ndata2); %spectral matrix in coordinate defined by V axes
  dop = zeros(1,ndata2);


  %% average spectral matrix over some number of wave period

% % %   num_wavePeriod = 4;
% % %   sampl_av = fix(num_wavePeriod*sampl/newfreqmat);
% % %   if sampl_av/2 == floor(sampl_av/2)
% % %     sampl_av=sampl_av+1;
% % %   end
% % %   %display(newfreqmat);
% % %   %display(sampl_av);
% % % 
% % % %   for i = sampl_av:ndata-sampl_av,
% % % %       if sampl_av < 2
% % % %           SMpermute(i,:,:)=SMpermute(i,:,:);
% % % %       elseif sampl_av > 2 && sampl_av < 4
% % % %           SMpermute(i,:,:) = (SMpermute(i-1,:,:)+SMpermute(i,:,:)+...
% % % %               SMpermute(i+1,:,:))./sampl_av;
% % % %       else 
% % % %           SMpermute(i,:,:) = sum(SMpermute(i-((sampl_av-1)/2):i+...
% % % %               ((sampl_av-1)/2),:,:),1)./sampl_av;
% % % %       end
% % % %   end     
% % %   for i = 1:(sampl_av-1)/2,
% % %       SMpermute(i,:,:) = mean(SMpermute(1:(sampl_av-1)/2-1,:,:),1);
% % %   end
% % %   for i = ndata-(sampl_av-1)/2:ndata,
% % %       SMpermute(i,:,:) = mean(SMpermute(ndata-(sampl_av-1)/2:ndata,:,:),1);
% % %   end
% % %   for i = (sampl_av-1)/2+1:ndata-(sampl_av-1)/2,
% % %       if sampl_av < 2
% % %           SMpermute(i,:,:)=SMpermute(i,:,:);
% % %       elseif sampl_av > 2 && sampl_av < 4
% % %           SMpermute(i,:,:) = (SMpermute(i-1,:,:)+SMpermute(i,:,:)+...
% % %               SMpermute(i+1,:,:))./sampl_av;
% % %       else 
% % %           SMpermute(i,:,:) = mean(SMpermute(i-((sampl_av-1)/2):i+...
% % %               ((sampl_av-1)/2),:,:),1);
% % %       end
% % %   end            

    %% average spectral matrix over four wave periods
  SMpermute_av = SMpermute;
  if sampl_av == 2,
      for i=1:ndata2-1,
          SMpermute_av(i,:,:) = mean(SMpermute(i:i+1,:,:),1);
      end
  end
  if sampl_av > 2,
    %sampl_half = int16(floor(sampl_av/2));
    sampl_half = int32(floor(sampl_av/2));
    %display(sampl_half);
    if sampl_av < ndata2,
        for i=1:sampl_half,
            SMpermute_av(i,:,:) = mean(SMpermute(1:i+sampl_half,:,:),1);
        end
        for i=sampl_half+1:ndata2-sampl_half-1,
            SMpermute_av(i,:,:) = mean(SMpermute(i-sampl_half:i+sampl_half,:,:),1);
        end
        for i=ndata2-sampl_half:ndata2,
            SMpermute_av(i,:,:) = mean(SMpermute(i-sampl_half:end,:,:),1);
        end
    elseif sampl_av >= ndata2 && sampl_half < ndata2,
        for i=1:sampl_av-ndata2,
            SMpermute_av(i,:,:) = mean(SMpermute(1:i+sampl_half,:,:),1);
        end
        for i=sampl_av-ndata2+1:ndata2-(sampl_av-ndata2),
            SMpermute_av = mean(SMpermute,1);
        end
        for i=ndata2-(sampl_av-ndata2)+1:ndata2,
            %display(i-sampl_half);
            %display(i);
            SMpermute_av(i,:,:) = mean(SMpermute(i-sampl_half+1:end,:,:),1);
        end
    elseif sampl_av > ndata2 && sampl_half >= ndata2,
        SMpermute_av = mean(SMpermute,1);
    end
    
  end  
  SMpermute = SMpermute_av;
  %display(any(isnan(SMpermute)));
  SMpermute(isnan(SMpermute))=0;

    %% compute singular value decomposition
%display(size(A));
%display(size(SMpermute));
  A(1:3,:,:) = real(permute(SMpermute,[2,3,1]));
  A(4:6,:,:) = -imag(permute(SMpermute,[2,3,1]));

  %for i = 1:size(Wb,1)
  for i = 1:ndata2,
     [U(:,:,i),W(:,:,i),V(:,:,i)] = svd(A(:,:,i),0);
%     wSingularValues(:,i) = svd(A(:,:,i),0);
  %   SMdet(i) = det(SM(:,:,i))
  end
  %censur=floor(2*a);
  %censur_indexes=[1:min(censur(ind_a),size(e,1)) max(1,size(e,1)-censur(ind_a)):size(e,1)];
  %W(:,:,censur_indexes) = NaN;

  
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
      S = zeros(ndata2,3);
      Wex=We(:,1);Wey=We(:,2);Wez=We(:,3);
      Wbx=Wb2(:,1);Wby=Wb2(:,2);Wbz=Wb2(:,3);
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
   
  %following lines cause issue - need to work this out
       %% Remove data possibly influenced by edge effects
   censur=floor(2*a);
   censur_indexes=[1:min(censur(ind_a),size(e,1)) max(1,size(e,1)-censur(ind_a)):size(e,1)];
   %note! the following lines are repeated below the parfor
   if pc12_range || default_range,
     censur2=floor(.4*a);
   end
   if pc35_range,
       censur2=floor(.1*a);
   end
   censur_indices2=[1:min(censur2(ind_a),length(t1)) max(1,length(t1)-censur2(ind_a)):length(t1)];
   
%    if wantEE, 
%        powerE(censur_indexes,:) = NaN;
%        powerEISR2(censur_indexes,:) = NaN;
%    end
%    powerB(censur_indexes,:) = NaN;
%    if wantEE, 
%        Spar(censur_indexes) = NaN;
%        Sx(censur_indexes) = NaN;
%        Sy(censur_indexes) = NaN;
%    end
   
   if wantEE, 
       powerE(censur_indices2,:) = NaN;
       powerEISR2(censur_indices2,:) = NaN;
   end
   powerB(censur_indices2,:) = NaN;
   if wantEE, 
       Spar(censur_indices2) = NaN;
       Sx(censur_indices2) = NaN;
       Sy(censur_indices2) = NaN;
   end
   
   SMpermute(censur_indices2,:,:) = NaN;
   if wantPolarization, 
     theta(censur_indices2) = NaN;
     phi(censur_indices2) = NaN;
     W(:,:,censur_indices2) = NaN;
   end
   planarity(:,ind_a) = 1 - sqrt(W(3,3,:)./W(1,1,:)); %planarity of polarization
%  Lp(:,ind_a) = W(2,2,:)./W(1,1,:); %ratio of two axes of polarization ellipse
%  planarity(:,ind_a) = 1 - sqrt(wSingularValues(3,:)./wSingularValues(1,:)); %planarity of polarization
  if wantPolarization,
%    polarizationEllipseRatio(:,ind_a) = wSingularValues(2,:)./wSingularValues(1,:); %ratio of two axes of polarization ellipse
    polarizationEllipseRatio(:,ind_a) = W(2,2,:)./W(1,1,:); %ratio of two axes of polarization ellipse
  end
   
  %% Calculate polarization parameters
  
   for i = 1:ndata2,
   %for i = 3500:3502,
%      VV=squeeze(V(:,:,i));
%      R(:,:,i)=VV*squeeze(spectralMatrix(:,:,i))*transpose(VV);
%      R(1,1,i)=real(R(1,1,i));
%      R(2,2,i)=real(R(2,2,i));
%      R(3,3,i)=real(R(3,3,i));
     SMsqueeze = squeeze(SMpermute(i,:,:));
 %    Rsqueeze = squeeze(R(:,:,i));
     
     dop(:,i) = (3/2.*trace(real(SMsqueeze)^2)./(trace(SMsqueeze))^2 - 1/2);
     %%%dop(:,i) = (2.*trace(real(Rsqueeze(1:2,1:2))^2)./(trace((Rsqueeze(1:2,1:2))))^2 - 1);
     %%%dop(:,i) = (2.*trace(real(SMsqueeze(1:2,1:2))^2)./(trace((SMsqueeze(1:2,1:2))))^2 - 1);
     %%%dop(:,i) = (3/2.*trace(real(Rsqueeze)^2)./(trace(Rsqueeze))^2 - 1/2);
   end
   dop(:,censur_indices2) = NaN;
%     if ind_a == 10,
%         SMpermute(3500,:,:)
%         R(:,:,3500)
% %        V(:,:,3500)
%         dop(:,3500)
%     end
%   if any(SMpermute(:,1,1) ~= real(SMpermute(:,1,1))),
%     display('imag');
%   end
          
  if wantPolarization,
    thetaSVD_fac(:,ind_a) = theta;
    phiSVD_fac(:,ind_a) = phi;
    polarizationSign(:,ind_a) = sign(imag(SMpermute(:,1,2))); %sign of polarization
      %Ls4(:,ind_a) = imag(SMpermute(:,1,2))./sqrt((imag(SMpermute(:,1,2))).^2+(imag(SMpermute(:,1,3))).^2+(imag(SMpermute(:,2,3))).^2);

     % SMdet = SM(1,1,:).*(SM(2,2,:).*SM(3,3,:)-SM(2,3,:).*SM(3,2,:))-SM(2,1,:).*(SM(1,2,:).*SM(3,3,:)-SM(1,3,:).*SM(3,2,:))+SM(1,3,:).*(SM(1,2,:).*SM(2,3,:)-SM(1,3,:).*SM(2,2,:));
     % Ls6(:,ind_a) = real(sqrt(1-(4*SMdet)./(SM(1,1,:)+SM(2,2,:)+SM(3,3,:)).^2));
%    degreeOfPolarization(:,ind_a) = sqrt(3/2.*(SMpermute(:,1,1).*...
%        SMpermute(:,1,1)+SMpermute(:,2,2).*SMpermute(:,2,2)+...
%        SMpermute(:,3,3).*SMpermute(:,3,3))./(SMpermute(:,1,1)+...
%        SMpermute(:,2,2)+SMpermute(:,3,3)).^2-1/2); %degree of polarization
        saveSMpermute=SMpermute;
        SMpermute=real(SMpermute);
        Ssquared = zeros(ndata2,3,3); 
        Ssquared(:,1,1) = SMpermute(:,1,1).*SMpermute(:,1,1)+SMpermute(:,1,2).*SMpermute(:,2,1)+SMpermute(:,1,3).*SMpermute(:,3,1);
        Ssquared(:,2,1) = SMpermute(:,2,1).*SMpermute(:,1,1)+SMpermute(:,2,2).*SMpermute(:,2,1)+SMpermute(:,2,3).*SMpermute(:,3,1);
        Ssquared(:,3,1) = SMpermute(:,3,1).*SMpermute(:,1,1)+SMpermute(:,3,2).*SMpermute(:,2,1)+SMpermute(:,3,3).*SMpermute(:,3,1);
        
        Ssquared(:,1,2) = SMpermute(:,1,1).*SMpermute(:,1,2)+SMpermute(:,1,2).*SMpermute(:,2,2)+SMpermute(:,1,3).*SMpermute(:,3,2);
        Ssquared(:,2,2) = SMpermute(:,2,1).*SMpermute(:,1,2)+SMpermute(:,2,2).*SMpermute(:,2,2)+SMpermute(:,2,3).*SMpermute(:,3,2);
        Ssquared(:,3,2) = SMpermute(:,3,1).*SMpermute(:,1,2)+SMpermute(:,3,2).*SMpermute(:,2,2)+SMpermute(:,3,3).*SMpermute(:,3,2);
        
        Ssquared(:,1,3) = SMpermute(:,1,1).*SMpermute(:,1,3)+SMpermute(:,1,2).*SMpermute(:,2,3)+SMpermute(:,1,3).*SMpermute(:,3,3);
        Ssquared(:,2,3) = SMpermute(:,2,1).*SMpermute(:,1,3)+SMpermute(:,2,2).*SMpermute(:,2,3)+SMpermute(:,2,3).*SMpermute(:,3,3);
        Ssquared(:,3,3) = SMpermute(:,3,1).*SMpermute(:,1,3)+SMpermute(:,3,2).*SMpermute(:,2,3)+SMpermute(:,3,3).*SMpermute(:,3,3);
        
        %%%degreeOfPolarization(:,ind_a) = sqrt(3/2.*(Ssquared(:,1,1)+Ssquared(:,2,2)+Ssquared(:,3,3)) ...
           %%% ./(SMpermute(:,1,1)+SMpermute(:,2,2)+SMpermute(:,3,3)).^2-1/2); %degree of polarization
        SMpermute=saveSMpermute;
   % degreeOfPolarization(:,ind_a) = sqrt(abs(2*(SMpermute(:,1,1).*SMpermute(:,1,1)+...
    %    SMpermute(:,2,2).*SMpermute(:,2,2)+2*(abs(SMpermute(:,1,2)).*abs(SMpermute(:,1,2))) ...
     %   ./(SMpermute(:,1,1)+SMpermute(:,2,2)).^2)-1));
%     R=real(R);
%     degreeOfPolarization(:,ind_a) = sqrt(2*(real(R(1,1,:)).*real(R(1,1,:))+real(R(2,2,:)).*real(R(2,2,:))+ ...
%         2*abs(real(R(1,2,:)).*real(R(1,2,:))))./(R(1,1,:)+R(2,2,:)).^2-1);
    degreeOfPolarization(:,ind_a) = dop;
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
%degreeOfPolarization=planarity;
%display('note that DOP is the planarity for now');


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

ndata2=size(power2B_SM_plot,1);
% I don't know why these lines aren't working
% if pc12_range || default_range,
%   censur2=floor(.5*a);
% end
% if pc35_range,
%   censur2=floor(.1*a);
% end
censur2=floor(.4*a);

for i=1:length(idx_nan_b)-1,
    if idx_nan_b(i) < idx_nan_b(i+1),
        display('front edge');
        for j=1:length(a),
            censur_index_front=[max(i-censur2(j),1):i];
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
            censur_index_back=[i:min(i+censur2(j),ndata2)];
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

if nargout==4,
	%timeVector = e(:,1);
    timeVector = fix(t1);
	frequencyVector = newfreq;
    BVector = Btot;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
elseif nargout==8,
	%timeVector = e(:,1);
    timeVector = fix(t1);
	frequencyVector = newfreq;
    BVector = Btot;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
    EESum_xxyyzz_ISR2 = power2E_ISR2_plot;
    EE_xxyyzz_FAC(:,:,4) = power2E_plot;
    EE_xxyyzz_FAC(:,:,1) = powerEx_plot
    EE_xxyyzz_FAC(:,:,2) = powerEy_plot
    EE_xxyyzz_FAC(:,:,3) = powerEz_plot
    Poynting_xyz_FAC = S_plot_x;
    Poynting_xyz_FAC(:,:,2) = S_plot_y;
    Poynting_xyz_FAC(:,:,3) = Spar_plot_z;
    Poynting_rThetaPhi_FAC = S_r;
    Poynting_rThetaPhi_FAC(:,:,2) = pi/2-S_elevation;
    Poynting_rThetaPhi_FAC(:,:,3) = S_azimuth;
else
 	%timeVector = e(:,1);
    timeVector = fix(t1);
	frequencyVector = newfreq;
    BVector = Btot;
    BB_xxyyzz_fac = powerBx_SM_plot;
    BB_xxyyzz_fac(:,:,2) = powerBy_SM_plot;
    BB_xxyyzz_fac(:,:,3) = powerBz_SM_plot;
	BB_xxyyzz_fac(:,:,4) = power2B_SM_plot;
    EESum_xxyyzz_ISR2 = power2E_ISR2_plot;
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