function moments = hiresmoments(Distpart,SCpot)
% PARTIALMOMENTS Compute moments from partial distribution function
% Written by D. B. Graham
% 
% moments = mms.partialmoments(Distpart,SCpot)
%
% Input:
%   Distpart - partial particle distribution. Must be PDist format. 
%   SCpot - Spacecraft potential
%
% Output:
%   moments - structure containing TSeries' of particle moments
%
% See also: mms.get_hiresdistributions

factor = Distpart.ancillary.pf;

Distpart.data = Distpart.data.*factor;
energyminus = Distpart.ancillary.delta_energy_minus;
energyplus = Distpart.ancillary.delta_energy_plus;

nonnanvalues = ~isnan(Distpart.data);
phi = irf.ts_scalar(Distpart.time,Distpart.depend{1,2});
theta = Distpart.depend{1,3};
stepTable = irf.ts_scalar(Distpart.time,Distpart.ancillary.esteptable);

moments = mms.psd_moments(Distpart,phi,theta,stepTable,Distpart.ancillary.energy0,Distpart.ancillary.energy1,... 
  SCpot,Distpart.species,'partialmoms',nonnanvalues,'energy_minus',energyminus,'energy_plus',energyplus);

end