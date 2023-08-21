function out = irf_4_v(r1,r2,r3,r4,x)
%IRF_4_V Calculate velocity or time shift of discontinuity.
%
%   v = IRF_4_V(r1,r2,r3,r4,t) Calculates velocity given position data
%   r1...r4 in Cluster format and time vector t.
%
%   dt = IRF_4_V(r1,r2,r3,r4,v) Calculates time shift from SC1 given
%   position and velocity of the discontinuity. v = [t,vx,vy,vz].
%
%   See also: IRF_4_V_GUI


% If the value is more than c, it is considered to be time.
% If EpochTT object then also time
if isa(x,'GenericTimeArray')
  flag='v_from_t';
  x = x.epochUnix;
elseif x(2) > 299792.458
  flag='v_from_t';
else
  flag='dt_from_v';
end

% check for TSeris input
if isa(r1,'TSeries'), r1 = [r1.time.epochUnix double(r1.data)]; end
if isa(r2,'TSeries'), r2 = [r2.time.epochUnix double(r2.data)]; end
if isa(r3,'TSeries'), r3 = [r3.time.epochUnix double(r3.data)]; end
if isa(r4,'TSeries'), r4 = [r4.time.epochUnix double(r4.data)]; end


switch flag
  case 'v_from_t'
    % Time input, velocity output
    t = x;

    dR = get_vol_ten(r1,r2,r3,r4,t);
    dt = [t(2),t(3),t(4)]-t(1);
    m = dR\dt'; % "1/v vector"

    V = m/norm(m)^2;
    out = V';

  case 'dt_from_v'
    % Time and velocity input, time output
    tc = x(1); % Center time
    v = x(2:4); % Input velocity

    m = v/norm(v)^2; % "1/v vector"
    if size(m,2) == 3
      m = m';
    end

    dR = get_vol_ten(r1,r2,r3,r4,tc);

    dt = dR*m;
    out = [0,dt'];
end
end

function dR = get_vol_ten(r1,r2,r3,r4,t)
% Function that returns volumetric tensor with sc1 as center. t can be
% vector or scalar.

if length(t) == 1
  t = [t,t,t,t];
end

% Interpolate position
r1=[t(1) interp1(r1(:,1),r1(:,[2 3 4]),t(1),'spline','extrap')];
r2=[t(2) interp1(r2(:,1),r2(:,[2 3 4]),t(2),'spline','extrap')];
r3=[t(3) interp1(r3(:,1),r3(:,[2 3 4]),t(3),'spline','extrap')];
r4=[t(4) interp1(r4(:,1),r4(:,[2 3 4]),t(4),'spline','extrap')];

r1 = r1(2:4);
r2 = r2(2:4);
r3 = r3(2:4);
r4 = r4(2:4);

% Volumetric tensor with SC1 as center.
dR = [r2;r3;r4]-[r1;r1;r1]; %mabye good?


end