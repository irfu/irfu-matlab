function TS = ts_skymap(time,data,energy,phi,theta,varargin)
%IRF.TS_SKYMAP  Factory for scalar (TSeries)
%
% TsSkymap = irf.ts_skymap(time,data,energy,phi,theta,'optionalArgs1',Val1,...)
%
% Create TSeries object - skymap

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else, epoch = time;
end
args = varargin;
if isempty(energy) % must check that energy0, energy1 and esteptable is given
  energy0_ok = 0;
  energy1_ok = 0;
  esteptable_ok = 0;
  while ~isempty(args)
    x = args{1}; args(1) = [];
    if ischar(x)
      switch lower(x)
        case {'energy0'}; energy0_ok = 1; energy0 = args{1}; args(1) = [];
        case {'energy1'}; energy1_ok = 1; energy1 = args{1}; args(1) = [];
        case {'esteptable'}; esteptable_ok = 1; esteptable = args{1}; args(1) = [];
      end
    end
  end
  if ~(energy0_ok && energy1_ok && esteptable_ok); error('Energy input required'); end
  
  if ~isnumeric(energy0); energy0 = energy.data; end
  if ~isnumeric(energy1); energy1 = energy.data; end
  if ~isnumeric(esteptable); esteptable = esteptable.data; end
  
  energy = repmat(torow(energy0),numel(esteptable),1);
  energy(esteptable==1,:) = repmat(energy1,sum(esteptable),1);
end
if ~isnumeric(phi); phi = phi.data; end
if ~isnumeric(theta); theta = theta.data; end

TS = PDist(epoch,data,'skymap',energy,phi,theta,varargin{:});
TS.ancillary.energy = energy;