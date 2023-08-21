function [mindist,nvec] = magnetopause_normal(pos_Re_gsm, IMF_Bz_nT, swp_nPa, modelFlag, Alfv_Mach)

% MODEL.MAGNETOPAUSE_NORMAL the distance and normal vector to the magnetopause
% for Shue et al., 1997 or Shue et al., 1998 model. Or bow shock for Farris
% & Russell 1994 model.
%
% [mindist,nvec] = MODEL.MAGNETOPAUSE_NORMAL(pos_Re_gsm, IMF_Bz_nT, swp_nPa, modelFlag, Alfv_Mach)
%
% Input:
%       pos_Re_gsm - GSM position in Re (if more than 3 values assumes that 1st is time)
%       IMF_Bz_nT  - IMF Bz in nT
%       swp_nPa    - Solar wind dynamic pressure in nPa
%       modelflag  - Set to 1 to use the 1998 model, otherwise the 1997
%                    model is used. Alternatively, specify name of model:
%                       'mp_shue1997'    - Shue et al., 1997
%                       'mp_shue1998'    - Shue et al., 1998
%                       'bs' or 'bs97'   - Bow shock, Farris & Russell 1994,
%                                          based on 'mp_shue1997'
%                       'bs98'           - Bow shock, Farris & Russell 1994,
%                                          based on 'mp_shue1998'
%       Alfv_Mach   - Alfvenic Mach number, only needed if bow shock model
%                     is used.
%
% Output:
%       mindist - minimum distance to the magnetopause, in Re. Positive
%       value if spacecraft is inside the magnetopause, negative if
%       outside the magnetopause.
%       nvec    - normal vector to the magnetopause (pointing away from
%       Earth).

% TODO: vectorize, so that input can be vectors

if nargin == 0
  help model.magnetopause_normal;return;
elseif nargin==3
  modelFlag = 'mp_shue1997';
  Alfv_Mach = 4;
elseif nargin==4
  Alfv_Mach = 4;
elseif nargin~=5
  irf.log('critical','Wrong number of input parameters, see help.');
  return;
end

% Keep old syntax
if isnumeric(modelFlag)
  if modelFlag == 1
    modelFlag = 'mp_shue1998';
  else
    modelFlag = 'mp_shue1997';
  end
end


if strcmpi(modelFlag,'mp_shue1998') || strcmpi(modelFlag,'bs98')
  alpha = (0.58 -0.007*IMF_Bz_nT)*(1.0 +0.024*log(swp_nPa));
  r0 = (10.22 + 1.29*tanh(0.184*(IMF_Bz_nT + 8.14)))*swp_nPa^(-1.0/6.6);
  irf.log ('warning','Shue et al., 1998 model used.')
else
  alpha = ( 0.58 -0.01*IMF_Bz_nT)*( 1.0 +0.01*swp_nPa);

  if IMF_Bz_nT>=0, r0 = (11.4 +0.013*IMF_Bz_nT)*swp_nPa^(-1.0/6.6);
  else,            r0 = (11.4 +0.140*IMF_Bz_nT)*swp_nPa^(-1.0/6.6);
  end
  irf.log ('warning','Shue et al., 1997 model used.')
end

%SC pos
x1 = pos_Re_gsm(1);
y1 = pos_Re_gsm(2);
z1 = pos_Re_gsm(3);

x0 = x1;
y0 = sqrt(y1^2+z1^2);

if strcmpi(modelFlag(1:2),'mp')
  % Magnetopause

  d2 = @(theta) (r0)^2*(2./(1+cos(theta))).^(2*alpha)...
    - 2*r0*(2./(1+cos(theta))).^(alpha).*(x0*cos(theta) + y0*sin(theta))...
    + x0^2 + y0^2;

  [thetamin,minval] = fminbnd(d2,-pi/1.2,pi/1.2);
  mindist = sqrt(minval);

  %calculate the direction to the spacecraft normal to the magnetopause
  xn = r0*(2/(1+cos(thetamin)))^alpha * cos(thetamin) - x1;
  phi = atan2(z1,y1);
  yn = cos(phi)*(r0*(2/(1+cos(thetamin)))^alpha * sin(thetamin)) - y1;
  zn = sin(phi)*(r0*(2/(1+cos(thetamin)))^alpha * sin(thetamin)) - z1;

  %disttest = sqrt(xn^2+yn^2+zn^2)

  nvec = [xn yn zn]/mindist;

  %if statement to ensure normal is pointing away from Earth
  if (sqrt(x0^2+y0^2) > r0*(2/(1+cos(thetamin)))^alpha)
    nvec = -nvec;
    mindist = -mindist;
  end

else
  % Bow shock
  irf.log ('warning','Farris & Russell 1994 bow shock model used.')
  gamma = 5/3;
  M = Alfv_Mach;
  % Bow shock standoff distance
  rbs = r0*(1+1.1*((gamma-1)*M^2+2)/((gamma+1)*(M^2-1)));
  % y^2 = -Ax + Bx^2
  A = 45.3;
  B = 0.04;
  x = rbs:-0.001:-100;

  y = sqrt(-A*(x-rbs) + B*(x-rbs).^2);
  x = [fliplr(x),x];
  y = [-fliplr(y),y];

  d2 = (x-x0).^2 + (y-y0).^2;

  [minval,minpos] = min(d2);
  d = [x'-x0,y'-y0];
  dmin = d(minpos,:);
  xn = dmin(1)/norm(dmin);
  mindist = sqrt(minval);

  q = y1/z1;
  zn = sign(z1)*sign(xn)*sqrt((1-xn^2)/(1+q^2));
  yn = zn*q;

  nvec = [xn,yn,zn];

  %if statement to ensure normal is pointing away from Earth
  if nvec(1)<0
    nvec = -nvec;
    mindist = -mindist;
  end

end
end