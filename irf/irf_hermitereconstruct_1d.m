function reconstruct = irf_hermitereconstruct_1d(hermitestruct,order)
% IRF_HERMITERECONSTRUCT_1D Reconstruct 1D distributions from hermite
% function coefficients
% 
% Input:
%       hermitestruct - structure of Hermite function coefficients from
%       irf_hermitedecomp_1d
%       order - the highest or Hermite function used in the reconstruction
% 
% Output: 
%       reconstruct - stucture containing distributions constructed from
%       Hermite functions and coefficients 
% 
% Written by D. B. Graham

coeffs = hermitestruct.hermitespectra;
vmat = hermitestruct.vmat*1e3;
vth = hermitestruct.vth*1e3;
vbulk = hermitestruct.Vmoms*1e3;

lengthv = length(vmat(1,:));
%dv = median(diff(vmat(1,:)));

vbulkmat = repmat(vbulk,1,lengthv);
vthmat = repmat(vth,1,lengthv);
vnormmat = (vmat-vbulkmat)./vthmat;
sqrtvthmat = sqrt(vthmat);

Pdistr1D = zeros(size(vmat));

for ii=0:order
  coeffmattemp = repmat(coeffs(:,ii+1),1,lengthv);
  Pdistr1D = Pdistr1D + irf_hermitefunction(vnormmat,ii).*coeffmattemp./sqrtvthmat;
end

reconstruct = struct;
reconstruct.Pdistr1D = Pdistr1D;
reconstruct.orders = 0:order;
reconstruct.times = hermitestruct.times;
reconstruct.vmat = hermitestruct.vmat;
reconstruct.nmoms = hermitestruct.nmoms;
reconstruct.Vmoms = hermitestruct.Vmoms;
reconstruct.Pmoms = hermitestruct.Pmoms;
reconstruct.Tmoms = hermitestruct.Tmoms;
reconstruct.vth = hermitestruct.vth;
reconstruct.species = hermitestruct.species;

end