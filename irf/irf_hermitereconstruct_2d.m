function reconstruct = irf_hermitereconstruct_2d(hermitestruct,order)
% IRF_HERMITERECONSTRUCT_2D Reconstruct 2D distributions from hermite
% function coefficients
% 
% Input:
%       hermitestruct - structure of Hermite function coefficients from
%       irf_hermitedecomp_2d
%       order - the highest or Hermite function used in the reconstruction
% 
% Output: 
%       reconstruct - stucture containing distributions constructed from
%       Hermite functions and coefficients 
% 
% Written by D. B. Graham
tic

coeffs = hermitestruct.hermitespectra;
vxmat = hermitestruct.vxmat*1e3;
vxth = hermitestruct.vxth*1e3;
vxbulk = hermitestruct.Vxmoms*1e3;
vymat = hermitestruct.vymat*1e3;
vyth = hermitestruct.vyth*1e3;
vybulk = hermitestruct.Vymoms*1e3;

lengthvx = length(vxmat(1,:,1));
lengthvy = length(vymat(1,1,:));

vxbulkmat = repmat(vxbulk,1,lengthvx,lengthvy);
vybulkmat = repmat(vybulk,1,lengthvx,lengthvy);

vxthmat = repmat(vxth,1,lengthvx,lengthvy);
vxnormmat = (vxmat-vxbulkmat)./vxthmat;
sqrtvxthmat = sqrt(vxthmat);

vythmat = repmat(vyth,1,lengthvx,lengthvy);
vynormmat = (vymat-vybulkmat)./vythmat;
sqrtvythmat = sqrt(vythmat);

Pdistr2D = zeros(size(vxmat));

c_eval('hermitejj? = irf_hermitefunction(vynormmat,?)./sqrtvythmat;',0:order);

for ii=0:order
  hermiteii = irf_hermitefunction(vxnormmat,ii)./sqrtvxthmat;
  for jj=0:order
    coeffstemp = repmat(coeffs(:,ii+1,jj+1),1,lengthvx,lengthvy);
    c_eval('Pdistr2D = Pdistr2D + coeffstemp.*hermiteii.*hermitejj?;',jj)
  end
end

reconstruct = struct;
reconstruct.Pdistr2D = Pdistr2D;
reconstruct.orders = 0:order;
reconstruct.times = hermitestruct.times;
reconstruct.vxmat = hermitestruct.vxmat;
reconstruct.vymat = hermitestruct.vymat;
reconstruct.nmoms = hermitestruct.nmoms;
reconstruct.Vxmoms = hermitestruct.Vxmoms;
reconstruct.Vymoms = hermitestruct.Vymoms;
reconstruct.Pxxmoms = hermitestruct.Pxxmoms;
reconstruct.Pyymoms = hermitestruct.Pyymoms;
reconstruct.Pxymoms = hermitestruct.Pxymoms;
reconstruct.Txxmoms = hermitestruct.Txxmoms;
reconstruct.Tyymoms = hermitestruct.Tyymoms;
reconstruct.Txymoms = hermitestruct.Txymoms;
reconstruct.vxth = hermitestruct.vxth;
reconstruct.species = hermitestruct.species;

toc

end

