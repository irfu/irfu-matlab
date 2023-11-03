


function [Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(varargin)
%irf_wavepolarize_magneticSVD analysis the polarization of magnetic wave using "magnetic SVD" method
%
%***************************************************************************
%**if you find any bugs in this function, pls contact me (huishan@irfu.se)**
%***************************************************************************
%Input parameter (2 parameters or 3, or 4)
%1. Bwave (1st colomn is time; 2nd to 4th colomn is bx, by, bz)
%2. Bbgd  Background magnetic field (1st colomn is time; 2nd to 4th colomn is Bx, By, Bz)
%3. [minPsd]  -- (optional) thresold for the analysis (for example: 1.0e-7). Below this value, the SVD analysis is meaningless
%   if minPsd is not given, SVD analysis will be done for all waves
%4. [windwid]  -- (optional) window width for fft analysis (for example: 80).
%   if windwid is not given, it will be automatically set as: sampling rate/200
%
%Output:
%Bpsd => power spectra density of magnetic filed wave
%planarity => planarity of polarization (form 0 to 1)
%waveangle => (form 0 to 90)
%elliptict => (form -1 to 1)
%
%Reference: Santolik et al. 2003 [radio science]
%Notice: Bwave and Bbgd should be from the same satellite and in the same coordinates
%
%--------------------------------------------------------------------------------------
%Example:
%[Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(Bwave,Bbgd);
%[Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(Bwave,Bbgd,1.0e-7);
%[Bpsd,planarity,waveangle,elliptict]=irf_wavepolarize_magneticSVD(Bwave,Bbgd,1.0e-7,80);
%--------------------------------------------------------------------------------------

%----writen by Huishan Fu at IRFU (2012-08-27)----



[ax,args,nargs] = axescheck(varargin{:});
Bwave=args{1};
Bbgd=args{2};

fs=1/(Bwave(2,1)-Bwave(1,1));   % sampling rate in Hz
minPsd=1.0e-25;
windwid=fix(fs/200);   % length of each window in micro second
if windwid<40, windwid=40; end

if nargs==3
  minPsd=args{3};
end
if nargs==4
  minPsd=args{3};
  windwid=args{4};
end

window_overlap=windwid/2; % length of each window overlaps in micro second


%% fft
[Bxcomplex, freq, time]=irf_wavefft(Bwave(:,2), 'hamming', window_overlap, windwid, fs);
[Bycomplex, freq, time]=irf_wavefft(Bwave(:,3), 'hamming', window_overlap, windwid, fs);
[Bzcomplex, freq, time]=irf_wavefft(Bwave(:,4), 'hamming', window_overlap, windwid, fs);
time=linspace(Bwave(1,1),Bwave(end,1),length(time));


%% construct spectra matrix
nt=length(time);
nf=length(freq);
SpecMatrix(:,:,1,1)=Bxcomplex.' .* Bxcomplex';
SpecMatrix(:,:,1,2)=Bxcomplex.' .* Bycomplex';
SpecMatrix(:,:,1,3)=Bxcomplex.' .* Bzcomplex';
SpecMatrix(:,:,2,1)=Bycomplex.' .* Bxcomplex';
SpecMatrix(:,:,2,2)=Bycomplex.' .* Bycomplex';
SpecMatrix(:,:,2,3)=Bycomplex.' .* Bzcomplex';
SpecMatrix(:,:,3,1)=Bzcomplex.' .* Bxcomplex';
SpecMatrix(:,:,3,2)=Bzcomplex.' .* Bycomplex';
SpecMatrix(:,:,3,3)=Bzcomplex.' .* Bzcomplex';


%% average spectra matrix
nosmbins=7;                                       %No. of bins in frequency domain
aa=[0.024,0.093,0.232,0.301,0.232,0.093,0.024];   %smoothing profile based on Hanning
esm=SpecMatrix;
for ii=1:nt
  for jj=((nosmbins-1)/2+1) : (nf-(nosmbins-1)/2)
    esm(ii,jj,1,1)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),1,1));
    esm(ii,jj,2,1)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),2,1));
    esm(ii,jj,3,1)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),3,1));
    esm(ii,jj,1,2)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),1,2));
    esm(ii,jj,2,2)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),2,2));
    esm(ii,jj,3,2)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),3,2));
    esm(ii,jj,1,3)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),1,3));
    esm(ii,jj,2,3)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),2,3));
    esm(ii,jj,3,3)=sum(aa(1:nosmbins).*SpecMatrix(ii,(jj-(nosmbins-1)/2):(jj+(nosmbins-1)/2),3,3));
  end
end


%---------------------------------------------------
Bbgd=irf_resamp(Bbgd,time');  Babstmp=irf_abs(Bbgd);
Bmag=Babstmp(:, [1 5]);
%---------------------------------------------------


%% begin SVD
for ii=1:nt
  for jj=1:nf
    Bsmr=real(esm(ii,jj,:,:));  Bsmr=squeeze(Bsmr);
    Bsmi=imag(esm(ii,jj,:,:));  Bsmi=squeeze(Bsmi);
    X=[Bsmr; Bsmi];

    nb=Bbgd(ii,2:4)/Bmag(ii,2);
    nperp1=cross(nb,[0 1 0]);
    nperp1=nperp1/sqrt(dot(nperp1,nperp1));
    nperp2=cross(nb,nperp1);

    [U,W,V] = svd(X,'econ');

    rowmin=3;  rowmid=2;  rowmax=1;
    w1=W(rowmin,rowmin); w2=W(rowmid,rowmid); w3=W(rowmax,rowmax);

    VT=V';
    k_propagate=VT(rowmin,:);
    k_perp_min=VT(rowmid,:);
    k_perp_max=VT(rowmax,:);

    Bpsd(ii,jj)=Bsmr(1,1)+Bsmr(2,2)+Bsmr(3,3);
    if Bpsd(ii,jj)>minPsd
      waveangle(ii,jj)=acos(dot(k_propagate, nb))*180/pi;
      if waveangle(ii,jj)>90, waveangle(ii,jj)=180-waveangle(ii,jj); end
      planarity(ii,jj)=1-sqrt(w1/w3);
      ellipse(ii,jj)=w2/w3;
      %estimate the sense
      ss=sign(Bsmi(1,2));
      ellipse(ii,jj)=ss*ellipse(ii,jj);
    else
      waveangle(ii,jj)=nan;
      planarity(ii,jj)=nan;
      ellipse(ii,jj)=nan;
    end
  end
end



%% save as structure format
Bpsd=struct('t',time,'f',freq','p',Bpsd,'f_unit','Hz');
planarity=struct('t',time,'f',freq','p',planarity,'f_unit','Hz');
waveangle=struct('t',time,'f',freq','p',waveangle,'f_unit','Hz');
elliptict=struct('t',time,'f',freq','p',ellipse,'f_unit','Hz');






