function [S,Sz,intSz]=irf_Poynting_flux(E,B,Bo)
% function [S,Sz,intSz]=irf_Poynting_flux(E,B,Bo)
% function [S,Sz]=irf_Poynting_flux(E,B,Bo)
% function [Sz]=irf_Poynting_flux(E,B,Bo)
% function [S,intS]=irf_Poynting_flux(E,B)
% function [S]=irf_Poynting_flux(E,B)
%   estimates Poynting flux S and Poynting flux along Bo
%   from electric field E and magnetic field B
%
%   If E and B have different sampling then the lowest sampling
%   is resampled at the highest sampling
%
%   E  units [mV/m]
%   B  units [nT]
%   S  units [mW/m^2]
%intSz units [mW*s/m^2]

% check which Poynting flux to calculate
flag_Sz=0;flag_intSz=0;flag_intS=0;
if nargin == 2 && nargout == 1, flag_Sz=0;flag_intSz=0;flag_intS=0; end
if nargin == 2 && nargout == 2, flag_Sz=0;flag_intSz=0;flag_intS=1; end
if nargin == 3 && nargout == 1, flag_Sz=1;flag_intSz=0;flag_intS=0; end
if nargin == 3 && nargout == 2, flag_Sz=1;flag_intSz=0;flag_intS=0; end
if nargin == 3 && nargout == 3, flag_Sz=1;flag_intSz=1;flag_intS=0; end

% resample if necessary
Fs_E=1/(E(2,1)-E(1,1));Fs_B=1/(B(2,1)-B(1,1));
% interval where both E & B exist
int=[max([min(E(:,1)) min(B(:,1))]) min([max(E(:,1)) max(B(:,1))]) ];
ee=irf_tlim(E,int);
bb=irf_tlim(B,int);
if Fs_E<Fs_B,e=irf_resamp(ee,bb);b=bb;Fs=Fs_B;
elseif Fs_E>Fs_B, b=irf_resamp(bb,ee);e=ee;Fs=Fs_E;
else
  disp('assuming the same sampling. Interpolating B and E to 2x E sampling.');
  t=sort([ee(:,1);ee(:,1)+0.5/Fs_E]);e=irf_resamp(ee,t);b=irf_resamp(bb,t);Fs=2*Fs_E;
end

% Calculate Poynting flux
S=irf_tappl(irf_cross(e,b),'/(4*pi/1e7)*1e-9');

if flag_Sz
  bm=irf_resamp(Bo,e);
  Sz=irf_dot(irf_norm(bm),S);
  if nargout ==1,S=Sz;end
end

% time integral of Poynting flux along ambient magnetic field
if flag_intSz
 ssz=Sz;
 ssz(isnan(Sz))=0; % set to zero points where Sz=NaN
 intSz=ssz;
 intSz(:,2)=cumsum(ssz(:,2))/Fs;
end

if flag_intS% time integral of all Poynting flux components
 ss=S;
 ss(isnan(S))=0; % set to zero points where Sz=NaN
 intS=ss;
 intS(:,2:end)=cumsum(ss(:,2:end))/Fs;
 Sz=intS;
end
