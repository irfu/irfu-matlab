function moments = hiresmoments(varargin)
% PARTIALMOMENTS Compute moments from partial distribution function
% Written by D. B. Graham
%
% moments = mms.hiresmoments(Distpart,SCpot)
%
% Input:
%   Distpart - partial particle distribution. Must be PDist format.
%   SCpot - Spacecraft potential
%
% Optional Inputs:
%   'energyrange' - set energy range in eV to integrate over [E_min E_max].
%   energy range is applied to energy0 (if applicable) and the same elements are used for energy1 to
%   ensure that the same number of points are integrated over.
%
% Output:
%   moments - structure containing TSeries' of particle moments
%
% See also: mms.get_hiresdistributions

if (nargin < 2)
  nargin
  help hiresmoments;
  return;
end

Distpart = varargin{1};
SCpot = varargin{2};
Erange = [1 40e3];

args=varargin(3:end);

if numel(args)>0
  options=1;
else
  options=0;
end

while options
  l = 2;
  switch(lower(args{1}))
    case 'energyrange'
      if numel(args)>1 && isnumeric(args{2})
        Erange = args{2};
      end
    otherwise
      irf.log('critical',['Unknown flag: ' args{1}]);
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), options=0; end
end

factor = Distpart.ancillary.pf;
Distpart.data = Distpart.data.*factor;
nonnanvalues = ~isnan(Distpart.data);
moments = mms.psd_moments(Distpart,SCpot,'partialmoms',nonnanvalues,'energyrange',Erange);

end