

function [Bpsd,degpol,waveangle,elliptict,helict]=irf_wavepolarize_means(varargin)
%irf_wavepolarize_means analysis the polarization of magnetic wave using "means" method
%
%***************************************************************************
%**if you find any bugs in this function, pls contact me (huishan@irfu.se)**
%***************************************************************************
%Input parameter (2 parameters or 3, or 4)
%1. Bwave (1st colomn is time; 2nd to 4th colomn is bx, by, bz)
%2. Bbgd  Background magnetic field (1st colomn is time; 2nd to 4th colomn is Bx, By, Bz)
%3. [minPsd]  -- (optional) thresold for the analysis (for example: 1.0e-7). Below this value, the SVD analysis is meaningless
%   if minPsd is not given, SVD analysis will be done for all waves
%4. [nopfft]  -- (optional) Number of points in FFT (for example: 256).
%   if nopfft is not given, it will be automatically set as 256
%
%Output:
%Bpsd => power spectra density of magnetic filed wave
%degpol => degeree of polarization (form 0 to 1)
%waveangle => (form 0 to 90)
%elliptict => (form -1 to 1)
%helict
%
%Notice: Bwave and Bbgd should be from the same satellite and in the same coordinates
%
%-------------------------------------------------------------------------------------
%Example:
%[Bpsd,degpol,waveangle,elliptict,helict]=irf_wavepolarize_means(Bwave,Bbgd);
%[Bpsd,degpol,waveangle,elliptict,helict]=irf_wavepolarize_means(Bwave,Bbgd,1.0e-7);
%[Bpsd,degpol,waveangle,elliptict,helict]=irf_wavepolarize_means(Bwave,Bbgd,1.0e-7,256);
%-------------------------------------------------------------------------------------
%
%WARNING: If one component is an order of magnitude or more  greater than the other two then the polarisation results saturate and erroneously
%indicate high degrees of polarisation at all times and frequencies. Time series should be eyeballed before running the program.
%For time series containing very rapid changes or spikes the usual problems with Fourier analysis arise.
%Care should be taken in evaluating degree of polarisation results.
%For meaningful results there should be significant wave power at the frequency where the polarisation approaches
%100%. Remembercomparing two straight lines yields 100% polarisation.

%----writen by Huishan Fu at IRFU (2012-08-27)----


[ax,args,nargs]=axescheck(varargin{:});
Bwave=args{1};
Bbgd=args{2};

minPsd=1.0e-25;
nopfft=256;   % number of points in FFT

if nargs==3
  minPsd=args{3};
end
if nargs==4
  minPsd=args{3};
  nopfft=args{4};
end
steplength = nopfft/2;


nopoints=length(Bwave(:,1));
nosteps=(nopoints-nopfft)/steplength;             %total number of FFTs
nosmbins=7;                                       %No. of bins in frequency domain
aa=[0.024,0.093,0.232,0.301,0.232,0.093,0.024];   %smoothing profile based on Hanning


%% change wave to MFA coordinates
Bbgd=irf_resamp(Bbgd,Bwave);
for ii=1:length(Bwave(:,1))
  nb=irf_norm(Bbgd(ii,2:4));
  nperp1=cross(nb,[0 1 0]);
  nperp1=irf_norm(nperp1);
  nperp2=cross(nb,nperp1);

  Bz(ii)=dot(Bwave(ii,2:4),nb);
  Bx(ii)=dot(Bwave(ii,2:4),nperp1);
  By(ii)=dot(Bwave(ii,2:4),nperp2);
end
ct=Bwave(:,1)';


% DEFINE ARRAYS
xs=Bx;  ys=By;  zs=Bz;
sampfreq=1/(ct(2)-ct(1));
endsampfreq=1/(ct(nopoints)-ct(nopoints-1));
if sampfreq ~= endsampfreq
  disp(['Warning: file sampling '  'frequency changes',sampfreq,'Hz to',endsampfreq,'Hz' ]);
else
  disp(['ac '  'file sampling frequency',sampfreq,'Hz']);
end


for j=1:nosteps
  %FFT CALCULATION
  smooth=0.08+0.46*(1-cos(2*pi*(1:nopfft)/nopfft));
  tempx=smooth.*xs(((j-1)*steplength+1) : ((j-1)*steplength+nopfft));
  tempy=smooth.*ys(((j-1)*steplength+1) : ((j-1)*steplength+nopfft));
  tempz=smooth.*zs(((j-1)*steplength+1) : ((j-1)*steplength+nopfft));
  specx(j,:)=(fft(tempx));
  specy(j,:)=(fft(tempy));
  specz(j,:)=(fft(tempz));
  halfspecx(j,:)=specx(j,1:(nopfft/2));
  halfspecy(j,:)=specy(j,1:(nopfft/2));
  halfspecz(j,:)=specz(j,1:(nopfft/2));
  xs=circshift(xs,-steplength);
  ys=circshift(ys,-steplength);
  zs=circshift(zs,-steplength);

  %CALCULATION OF THE SPECTRAL MATRIX
  matspec(j,:,1,1)=halfspecx(j,:).*conj(halfspecx(j,:));
  matspec(j,:,2,1)=halfspecx(j,:).*conj(halfspecy(j,:));
  matspec(j,:,3,1)=halfspecx(j,:).*conj(halfspecz(j,:));
  matspec(j,:,1,2)=halfspecy(j,:).*conj(halfspecx(j,:));
  matspec(j,:,2,2)=halfspecy(j,:).*conj(halfspecy(j,:));
  matspec(j,:,3,2)=halfspecy(j,:).*conj(halfspecz(j,:));
  matspec(j,:,1,3)=halfspecz(j,:).*conj(halfspecx(j,:));
  matspec(j,:,2,3)=halfspecz(j,:).*conj(halfspecy(j,:));
  matspec(j,:,3,3)=halfspecz(j,:).*conj(halfspecz(j,:));


  %CALCULATION OF SMOOTHED SPECTRAL MATRIX
  ematspec(j,:,:,:)=matspec(j,:,:,:)*nan;
  for k=((nosmbins-1)/2+1) : (nopfft/2-(nosmbins-1)/2)
    ematspec(j,k,1,1)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),1,1));
    ematspec(j,k,2,1)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),2,1));
    ematspec(j,k,3,1)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),3,1));
    ematspec(j,k,1,2)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),1,2));
    ematspec(j,k,2,2)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),2,2));
    ematspec(j,k,3,2)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),3,2));
    ematspec(j,k,1,3)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),1,3));
    ematspec(j,k,2,3)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),2,3));
    ematspec(j,k,3,3)=sum(aa(1:nosmbins).*matspec(j,(k-(nosmbins-1)/2):(k+(nosmbins-1)/2),3,3));
  end

  %CALCULATION OF THE MINIMUM VARIANCE DIRECTION AND WAVENORMAL ANGLE
  aaa2(j,:)=sqrt(imag(ematspec(j,:,1,2)).^2+imag(ematspec(j,:,1,3)).^2+imag(ematspec(j,:,2,3)).^2);
  wnx(j,:)=abs(imag(ematspec(j,:,2,3))./aaa2(j,:));
  wny(j,:)=-abs(imag(ematspec(j,:,1,3))./aaa2(j,:));
  wnz(j,:)=imag(ematspec(j,:,1,2))./aaa2(j,:);
  waveangle(j,:)=atan(sqrt(wnx(j,:).^2+wny(j,:).^2)./abs(wnz(j,:)));

  %CALCULATION OF THE DEGREE OF POLARISATION
  %calc of square of smoothed spec matrix
  matsqrd(j,:,1,1)=ematspec(j,:,1,1).*ematspec(j,:,1,1)+ematspec(j,:,1,2).*ematspec(j,:,2,1)+ematspec(j,:,1,3).*ematspec(j,:,3,1);
  matsqrd(j,:,1,2)=ematspec(j,:,1,1).*ematspec(j,:,1,2)+ematspec(j,:,1,2).*ematspec(j,:,2,2)+ematspec(j,:,1,3).*ematspec(j,:,3,2);
  matsqrd(j,:,1,3)=ematspec(j,:,1,1).*ematspec(j,:,1,3)+ematspec(j,:,1,2).*ematspec(j,:,2,3)+ematspec(j,:,1,3).*ematspec(j,:,3,3);
  matsqrd(j,:,2,1)=ematspec(j,:,2,1).*ematspec(j,:,1,1)+ematspec(j,:,2,2).*ematspec(j,:,2,1)+ematspec(j,:,2,3).*ematspec(j,:,3,1);
  matsqrd(j,:,2,2)=ematspec(j,:,2,1).*ematspec(j,:,1,2)+ematspec(j,:,2,2).*ematspec(j,:,2,2)+ematspec(j,:,2,3).*ematspec(j,:,3,2);
  matsqrd(j,:,2,3)=ematspec(j,:,2,1).*ematspec(j,:,1,3)+ematspec(j,:,2,2).*ematspec(j,:,2,3)+ematspec(j,:,2,3).*ematspec(j,:,3,3);
  matsqrd(j,:,3,1)=ematspec(j,:,3,1).*ematspec(j,:,1,1)+ematspec(j,:,3,2).*ematspec(j,:,2,1)+ematspec(j,:,3,3).*ematspec(j,:,3,1);
  matsqrd(j,:,3,2)=ematspec(j,:,3,1).*ematspec(j,:,1,2)+ematspec(j,:,3,2).*ematspec(j,:,2,2)+ematspec(j,:,3,3).*ematspec(j,:,3,2);
  matsqrd(j,:,3,3)=ematspec(j,:,3,1).*ematspec(j,:,1,3)+ematspec(j,:,3,2).*ematspec(j,:,2,3)+ematspec(j,:,3,3).*ematspec(j,:,3,3);

  Trmatsqrd(j,:)=matsqrd(j,:,1,1)+matsqrd(j,:,2,2)+matsqrd(j,:,3,3);
  Trmatspec(j,:)=ematspec(j,:,1,1)+ematspec(j,:,2,2)+ematspec(j,:,3,3);
  degpol(j,:)=Trmatspec(j,:)*nan;
  degpol(j,((nosmbins-1)/2+1):(nopfft/2-(nosmbins-1)/2))=(3*Trmatsqrd(j,((nosmbins-1)/2+1):(nopfft/2-(nosmbins-1)/2))-Trmatspec(j,((nosmbins-1)/2+1):(nopfft/2-(nosmbins-1)/2)).^2)./(2*Trmatspec(j,((nosmbins-1)/2+1):(nopfft/2-(nosmbins-1)/2)).^2);

  %CALCULATION OF HELICITY, ELLIPTICITY AND THE WAVE STATE VECTOR
  alphax(j,:)=sqrt(ematspec(j,:,1,1));
  alphacos2x(j,:)=real(ematspec(j,:,1,2))./sqrt(ematspec(j,:,1,1));
  alphasin2x(j,:)=-imag(ematspec(j,:,1,2))./sqrt(ematspec(j,:,1,1));
  alphacos3x(j,:)=real(ematspec(j,:,1,3))./sqrt(ematspec(j,:,1,1));
  alphasin3x(j,:)=-imag(ematspec(j,:,1,3))./sqrt(ematspec(j,:,1,1));
  lambdau(j,:,1,1)=alphax(j,:);
  lambdau(j,:,1,2)=complex(alphacos2x(j,:),alphasin2x(j,:));
  lambdau(j,:,1,3)=complex(alphacos3x(j,:),alphasin3x(j,:));

  alphay(j,:)=sqrt(ematspec(j,:,2,2));
  alphacos2y(j,:)=real(ematspec(j,:,2,1))./sqrt(ematspec(j,:,2,2));
  alphasin2y(j,:)=-imag(ematspec(j,:,2,1))./sqrt(ematspec(j,:,2,2));
  alphacos3y(j,:)=real(ematspec(j,:,2,3))./sqrt(ematspec(j,:,2,2));
  alphasin3y(j,:)=-imag(ematspec(j,:,2,3))./sqrt(ematspec(j,:,2,2));
  lambdau(j,:,2,1)=alphay(j,:);
  lambdau(j,:,2,2)=complex(alphacos2y(j,:),alphasin2y(j,:));
  lambdau(j,:,2,3)=complex(alphacos3y(j,:),alphasin3y(j,:));

  alphaz(j,:)=sqrt(ematspec(j,:,3,3));
  alphacos2z(j,:)=real(ematspec(j,:,3,1))./sqrt(ematspec(j,:,3,3));
  alphasin2z(j,:)=-imag(ematspec(j,:,3,1))./sqrt(ematspec(j,:,3,3));
  alphacos3z(j,:)=real(ematspec(j,:,3,2))./sqrt(ematspec(j,:,3,3));
  alphasin3z(j,:)=-imag(ematspec(j,:,3,2))./sqrt(ematspec(j,:,3,3));
  lambdau(j,:,3,1)=alphaz(j,:);
  lambdau(j,:,3,2)=complex(alphacos2z(j,:),alphasin2z(j,:));
  lambdau(j,:,3,3)=complex(alphacos3z(j,:),alphasin3z(j,:));

  %HELICITY CALCULATION
  for k=1:nopfft/2
    for xyz=1:3
      upper(j,k)=sum(2*real(lambdau(j,k,xyz,1:3)).*(imag(lambdau(j,k,xyz,1:3))));
      lower(j,k)=sum((real(lambdau(j,k,xyz,1:3))).^2-(imag(lambdau(j,k,xyz,1:3))).^2);
      if upper(j,k)>0
        gamma(j,k)=atan(upper(j,k)/lower(j,k));
      else
        gamma(j,k)=pi+(pi+atan(upper(j,k)/lower(j,k)));
      end

      lambdau(j,k,xyz,:)=exp(complex(0,-0.5*gamma(j,k))).*lambdau(j,k,xyz,:);

      helicity(j,k,xyz)=1/(sqrt(real(lambdau(j,k,xyz,1))^2+real(lambdau(j,k,xyz,2))^2+real(lambdau(j,k,xyz,3))^2)/sqrt(imag(lambdau(j,k,xyz,1))^2+imag(lambdau(j,k,xyz,2))^2+imag(lambdau(j,k,xyz,3))^2));

      %ELLIPTICITY CALCULATION
      uppere=imag(lambdau(j,k,xyz,1))*real(lambdau(j,k,xyz,1))+imag(lambdau(j,k,xyz,2))*real(lambdau(j,k,xyz,2));
      lowere=-imag(lambdau(j,k,xyz,1))^2+real(lambdau(j,k,xyz,1))^2-imag(lambdau(j,k,xyz,2))^2+real(lambdau(j,k,xyz,2))^2;
      if uppere>0
        gammarot(j,k)=atan(uppere/lowere);
      else
        gammarot(j,k)=pi+pi+atan(uppere/lowere);
      end

      lam=lambdau(j,k,xyz,1:2);
      lambdaurot(j,k,:)=exp(complex(0,-0.5*gammarot(j,k)))*lam(:);

      ellip(j,k,xyz)=sqrt(imag(lambdaurot(j,k,1))^2+imag(lambdaurot(j,k,2))^2)/sqrt(real(lambdaurot(j,k,1))^2+real(lambdaurot(j,k,2))^2);
      ellip(j,k,xyz)=-ellip(j,k,xyz)*(imag(ematspec(j,k,1,2))*sin(waveangle(j,k)))/abs(imag(ematspec(j,k,1,2))*sin(waveangle(j,k)));


    end

  end

end  %end of main body


%AVERAGING HELICITY AND ELLIPTICITY RESULTS
elliptict=(ellip(:,:,1)+ellip(:,:,2)+ellip(:,:,3))/3;
helict=(helicity(:,:,1)+helicity(:,:,2)+helicity(:,:,3))/3;



%CREATING OUTPUT PARAMETER
timeline=ct(1)+abs(nopfft/2)/sampfreq+(1:nosteps)*steplength/sampfreq;
binwidth=sampfreq/nopfft;
freqline=binwidth*(1:nopfft/2);
%scaling power results to units with meaning;
W=nopfft*sum(smooth.^2);
powspec(:,2:nopfft/2-1)=1/W*2*Trmatspec(:,2:nopfft/2-1)/binwidth;
powspec(:,1)=1/W*Trmatspec(:,1)/binwidth;
powspec(:,nopfft/2)=1/W*Trmatspec(:,nopfft/2)/binwidth;



%KICK OUT THE ANALYSIS OF THE WEAK SIGNALS
index=find(powspec<minPsd);
waveangle(index)=nan;
degpol(index)=nan;
elliptict(index)=nan;
helict(index)=nan;



%SAVE DATA AS THE STRUCTURE FORMAT
Bpsd=struct('t',timeline,'f',freqline','p',powspec,'f_unit','Hz');
waveangle=struct('t',timeline,'f',freqline','p',waveangle*180/pi,'f_unit','Hz');
degpol=struct('t',timeline,'f',freqline','p',degpol,'f_unit','Hz');
elliptict=struct('t',timeline,'f',freqline','p',elliptict,'f_unit','Hz');
helict=struct('t',timeline,'f',freqline','p',helict,'f_unit','Hz');




end