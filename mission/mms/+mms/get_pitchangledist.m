function [paddist,theta] = get_pitchangledist(pdist,B) 
%
% paddist = mms.get_pitchangledist(pdist,B) 
%
%
% Computes the pitch angle distribution from l1b brst particle data. 
%
% Inputs:
%       pdist - electron or ion particle distribution
%       B     - magnetic field
% Output: 
%       paddist - Particle pitch angle distribution
%       theta - Pitch angles
%
% Written by D. B. Graham
%
% To Do: Incorperation angles for v1.0.0 data as an optional input. Also
% kinda slow.

B = B.resample(pdist);
Bmag = B.abs.data;
Bvec = B.data./[Bmag Bmag Bmag];
paddistarr = zeros(length(pdist.time),32,12);

% Define angles
dangle = 180/16;
phi = dangle*[0:31]+dangle/2;
theta = dangle*[0:15]+dangle/2;
[~,energy] = hist([log10(10),log10(30e3)],32);
energy = 10.^energy;

x = -cosd(phi')*sind(theta);
y = -sind(phi')*sind(theta);
z = -ones(length(phi),1)*cosd(theta);
z2 = ones(length(phi),1)*sind(theta);
solida = dangle*dangle*z2;
for jj=1:length(energy);
    allsolide(jj,:,:) = solida;
end
for ii = 1:length(pdist.time);
        allsolide2(ii,:,:,:) = allsolide;
end
pdist.data = pdist.data.*allsolide2;

anglevec = [15 30 45 60 75 90 105 120 135 150 165 180];
pitcha = anglevec-7.5;

for ii = 1:length(B.time);
    dist1 = squeeze(pdist.data(ii,:,:,:));
    thetab = acosd(x*Bvec(ii,1)+y*Bvec(ii,2)+z*Bvec(ii,3));
    c_eval('pos? = ones(length(phi),length(theta));',anglevec); 
    c_eval('sol? = solida;',anglevec);    
    pos15(thetab > 15) = NaN;
    sol15(thetab > 15) = NaN;
    for jj = 2:(length(anglevec)-1)
        c_eval('pos?(thetab < (?-15)) = NaN;',anglevec(jj));
        c_eval('pos?(thetab > ?) = NaN;',anglevec(jj));
        c_eval('sol?(thetab < (?-15)) = NaN;',anglevec(jj));
        c_eval('sol?(thetab > ?) = NaN;',anglevec(jj));
    end
    pos180(thetab < 165) = NaN;
    sol180(thetab < 165) = NaN;
    
    c_eval('dist? = dist1;',anglevec);
    
    for kk = 1:length(energy);
        c_eval('dist?(kk,:,:)  = squeeze(dist1(kk,:,:)).*pos?;',anglevec);
    end 
    
    c_eval('dist? =  squeeze(irf.nanmean(irf.nanmean(dist?,3),2))*1e30/irf.nanmean(irf.nanmean(sol?));',anglevec);
    paddistarr(ii,:,:) = [dist15 dist30 dist45 dist60 dist75 dist90 dist105 dist120 dist135 dist150 dist165 dist180];
end

paddist = TSeries(pdist.time,paddistarr);
theta = pitcha;

end