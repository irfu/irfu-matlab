function moments = hiresmoments(Distpart,SCpot)
% PARTIALMOMENTS Compute moments from partial distribution function
% Written by D. B. Graham
%
% moments = mms.hiresmoments(Distpart,SCpot)
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
nonnanvalues = ~isnan(Distpart.data);
moments = mms.psd_moments(Distpart,SCpot,'partialmoms',nonnanvalues);

end