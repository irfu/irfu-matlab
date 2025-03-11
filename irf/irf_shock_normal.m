function [nd] = irf_shock_normal(spec,leq90)
%IRF_SHOCK_NORMAL Calculates shock normals with different methods.
%
%   Normal vectors are calculated by methods described in ISSI Scientific
%   Report SR-001 Ch 10. (Schwartz 1998), and references therein.
%
%   nst = IRF_SHOCK_NORMAL(spec) returns structure nst which contains data
%   on shock normal vectors given input with plasma parameters spec. The
%   data can be averaged values or values from the time series in matrix
%   format. If series is from time series all parameters are calculated
%   from a random upstream and a random downstream point. This can help set
%   errorbars on shock angle etc. The time series input must have the same
%   size (up- and downstream can be different), so generally the user needs
%   to resample the data first.
%
%   nst = IRF_SHOCK_NORMAL(spec,leq90) for leq90 = 1 angles are forced to
%   be less than 90 (default). For leq90 = 0, angles can be between 0 and
%   180 deg. For time series input and quasi-perp shocks,leq90 = 0 is
%   recommended.
%
%   Input spec contains:
%
%       Bu  -   Upstream magnetic field (nT).
%       Bd  -   Downstream magnetic field.
%       Vu  -   Upstream plasma bulk velocity (km/s).
%       Vd  -   Downstream plasma bulk velocity.
%       nu  -   Upstream number density (cm^-3).
%       nd  -   Downstream number density.
%       Optional:
%       R   -   Spacecraft position in TSeries format, 1x3 vector or as
%               given by R = mms.get_data('R_gse',tint)
%       t   -   Time of shock crossing (EpochTT), for models.
%       d2u -   Down-to-up, is 1 or -1.
%       dTf -   Time duration of shock foot (s).
%       Fcp -   Reflected ion gyrofrequency (Hz).
%       N   -   Number of Monte Carlo particles used in determining
%               errorbars (default 100)
%       n0  -   User defined normal vector from e.g. timing
%
%
%   Output nd contains:
%
%       n   -   Structure containing normal vectors (n always points
%               toward the upstream region):
%           From data:
%               mc  -   Magnetic coplanarity (10.14)
%               vc  -   Velocity coplanarity (10.18)
%               mx1 -   Mixed method 1 (10.15), Abraham-Shrauner, B., 1972
%               mx2 -   Mixed method 2 (10.16), Abraham-Shrauner, B., 1972
%               mx3 -   Mixed method 3 (10.17), Abraham-Shrauner, B., 1972
%           Models (only if R is included in spec):
%               farris  -   Farris, M. H., et al., 1991
%               slho    -   Slavin, J. A. and Holzer, R. E., 1981 (corrected)
%               per     -   Peredo et al., 1995, (z = 0)
%               fa4o    -   Fairfield, D. H., 1971
%               fan4o   -   Fairfield, D. H., 1971
%               foun    -   Formisano, V., 1979
%           User defined (only if n0 is included in spec):
%               n0  -   User defined normal vector from e.g. timing
%
%       thBn -  Angle between normal vector and spec.Bu, same fields as n.
%
%       thVn -  Angle between normal vector and spec.Vu, same fields as n.
%
%       Vsh  -  Structure containing shock velocities:
%            Methods:
%               gt  -   Using shock foot thickness (10.32). Gosling, J. T., and M. F. Thomsen (1985)
%               mf  -   Mass flux conservation (10.29).
%               sb  -   Using jump conditions (10.33). Smith, E. J., and M. E. Burton (1988)
%               mo -    Using shock foot thickness
%
%       info -  Some more info:
%           msh     -   Magnetic shear angle
%           vsh     -   Velocity shear angle
%           cmat    -   Constraints matrix with normalized errors.
%           sig     -   Scaling factor to fit shock models to sc position
%           Calclated from (10.9-10.13) in (Schwartz 1998).
%
%
%   Examples:
%   shp = [];
%   shp.Bu = [.1,-.02,3];
%   shp.Bd = [1,.05,9];
%   shp.Vu = [-300,5,-2];
%   shp.Vd = [-100,5,-2];
%   shp.nu = 5;
%   shp.nd = 15;
%   nd = irf_shock_normal(shp);
%   nd.n
%
%   ans =
%
%       struct with fields:
%
%       mc: [-0.9767 -0.1552 0.1483]
%       vc: [1 0 0]
%       mx1: [0.9889 9.8901e-04 -0.1484]
%       mx2: [0.9889 -8.2406e-04 -0.1483]
%       mx3: [0.9889 -0.0017 -0.1483]
%
%
%   shp = [];
%   shp.Bu = [.1,-.02,3;.3,-.05,2.8];
%   shp.Bd = [1,.05,9;1,-.05,10];
%   shp.Vu = [-300,5,-2;-310,2,2];
%   shp.Vd = [-100,5,-2;-90,10,10];
%   shp.nu = [5;4];
%   shp.nd = [15;16];
%   shp.N = 5;
%   nd = irf_shock_normal(shp);
%   nd.n.mx3
%
% ans =
%
%     0.9904    0.0044   -0.1384
%     0.9917    0.0143   -0.1276
%     0.9924    0.0169   -0.1216
%     0.9903    0.0026   -0.1389
%     0.9922    0.0150   -0.1239
%
%
%   See also: IRF_SHOCK_GUI, IRF_SHOCK_PARAMETERS
%

%   Written by: Andreas Johlander, andreasj@irfu.se
%

%% Check if input are timeseries

% leq90 = 1 if angles should be less or equal to 90 deg. 0 is good if doing
% statistics
if nargin == 1
  leq90 = 1;
end

% checks only Bu and Bd (one can be a scalar)
if size(spec.Bu,1)>1 || size(spec.Bd,1)>1
  Nu = size(spec.Bu,1); % number of points, must be same for all parameters
  Nd = size(spec.Bd,1);

  % randomize points upstream and downstream
  if isfield(spec,'N'); N = spec.N; else; N = 10; end
  % sort of "indices" from 0 to 1
  idtu = rand(1,N); idtd = rand(1,N);

  % copy spec
  tempSpec = spec;


  % loop with up/downstream values randomly picked
  for i = 1:N
    % get value from given point
    tempSpec.Bu = interp1(linspace(0,1,Nu),spec.Bu,idtu(i));
    tempSpec.Vu = interp1(linspace(0,1,Nu),spec.Vu,idtu(i));
    tempSpec.nu = interp1(linspace(0,1,Nu),spec.nu,idtu(i));

    tempSpec.Bd = interp1(linspace(0,1,Nd),spec.Bd,idtd(i));
    tempSpec.Vd = interp1(linspace(0,1,Nd),spec.Vd,idtd(i));
    tempSpec.nd = interp1(linspace(0,1,Nd),spec.nd,idtd(i));

    if i == 1
      % run once to get structure
      nd = irf_shock_normal(tempSpec);
      % replace all content in fields with zeros of appropriate size
      fn1 = fieldnames(nd); % first layer names
      fn2 = [];
      % get field names for Vsh (if Vsh is ever changed, this
      % needs to change as well)
      fnVsh = fieldnames(nd.Vsh.mf); % only field with 3 layers
      % get second layer names
      for k = 1:length(fn1)
        % get second layer names, Vsh does not have two layers
        fn2.(fn1{k}) = fieldnames(nd.(fn1{k}));
      end

      % do a double loop
      for k = 1:length(fn1) % fist layer
        strl1 = fn1{k}; % first layer string
        if ~strcmp(strl1,'Vsh')
          for l = 1:length(fn2.(fn1{k})) % second layer
            strl2 = fn2.(fn1{k}){l}; % second layer string
            % this really messes up cmat
            % replace with zero-arrays
            nd.(strl1).(strl2) = zeros(N,size(nd.(strl1).(strl2),2));
          end
        else % special case for Vsh
          for l = 1:length(fn2.(fn1{k})) % second layer
            strl2 = fn2.(fn1{k}){l}; % second layer string
            for m = 1:length(fnVsh) % third layer
              nd.(strl1).(strl2).(fnVsh{m}) = zeros(N,size(nd.(strl1).(strl2).(fnVsh{m}),2));
            end
          end
        end
      end


    end

    % get shock normal for this point
    tempNd = irf_shock_normal(tempSpec,leq90);

    % do another double loop
    for k = 1:length(fn1) % fist layer
      strl1 = fn1{k}; % first layer string
      if ~strcmp(strl1,'Vsh')
        for l = 1:length(fn2.(fn1{k})) % second layer
          strl2 = fn2.(fn1{k}){l}; % second layer string
          % add value from temporary structure
          % skip difficult ones (cmat and sig)
          if ~strcmp(strl2,'cmat') && ~strcmp(strl2,'sig')
            nd.(strl1).(strl2)(i,:) = tempNd.(strl1).(strl2);
          end
        end
      else % special case for Vsh
        for l = 1:length(fn2.(fn1{k})) % second layer
          strl2 = fn2.(fn1{k}){l}; % second layer string
          for m = 1:length(fnVsh) % third layer
            nd.(strl1).(strl2).(fnVsh{m})(i,:) = tempNd.(strl1).(strl2).(fnVsh{m});
          end
        end
      end
    end


  end
  return; % return results with vector
end


%% Actual program (single input)

% normal vector, according to different models
nd = [];
n = [];

Bu = spec.Bu;
Bd = spec.Bd;
Vu = spec.Vu;
Vd = spec.Vd;

delB = Bd-Bu;
delV = Vd-Vu;
spec.delB = delB;
spec.delV = delV;

% magenetic coplanarity
n.mc = cross(cross(Bd,Bu),delB)/norm(cross(cross(Bd,Bu),delB));

% velocity coplanarity
n.vc = delV/norm(delV);

% Mixed methods
n.mx1 = cross(cross(Bu,delV),delB)/norm(cross(cross(Bu,delV),delB));
n.mx2 = cross(cross(Bd,delV),delB)/norm(cross(cross(Bd,delV),delB));
n.mx3 = cross(cross(delB,delV),delB)/norm(cross(cross(delB,delV),delB));

% user defined normal vector
if isfield(spec,'n0')
  n.n0 = spec.n0;
end

sig = [];
eps = [];
x0 = [];
y0 = [];
alpha = [];
% calculate model normals if R is inputted
if isfield(spec,'R')
  % Farris et al.
  [n.farris,sig.farris,eps.farris,L.farris,alpha.farris,x0.farris,y0.farris] = farris_model(spec);
  % Slavin and Holzer mean
  [n.slho,sig.slho,eps.slho,L.slho,alpha.slho,x0.slho,y0.slho] = slavin_holzer_model(spec);
  % Peredo et al., z = 0
  [n.per,sig.per,eps.per,L.per,alpha.per,x0.per,y0.per] = peredo_model(spec);
  % Fairfield Meridian 4o
  [n.fa4o,sig.fa4o,eps.fa4o,L.fa4o,alpha.fa4o,x0.fa4o,y0.fa4o] = fairfield_meridian_4o_model(spec);
  % Fairfield Meridian No 4o
  [n.fan4o,sig.fan4o,eps.fan4o,L.fan4o,alpha.fan4o,x0.fan4o,y0.fan4o] = fairfield_meridian_no_4o_model(spec);
  % Formisano Unnorm. z = 0
  [n.foun,sig.foun,eps.foun,L.foun,alpha.foun,x0.foun,y0.foun] = formisano_unnorm_model(spec);
end

% make sure all normal vectors are pointing upstream
% based on deltaV, should work for IP shocks also
fnames = fieldnames(n);
for ii = 1:length(fnames)
  if dot(delV,n.(fnames{ii}))<0; n.(fnames{ii}) = -n.(fnames{ii}); end
end

% shock angle
thBn = shock_angle(spec,n,'B',leq90);
thVn = shock_angle(spec,n,'V',leq90);

% Shock velocity from foot time
% Fcp in Hz as given by irf_plasma_calc.
% Shock foot timing (Gosling et al., 1982)
Vsh.gt = shock_speed(spec,n,thBn,'gt');
% Shock velocity from mass flux conservation
Vsh.mf = shock_speed(spec,n,thBn,'mf');
% Shock velocity from (Smith & Burton, 1988)
Vsh.sb = shock_speed(spec,n,thBn,'sb');
% Shock velocity from (Moses et al., 1985)
Vsh.mo = shock_speed(spec,n,thBn,'mo');

% info
info = [];
% magnetic shear angle
info.msh = shear_angle(Bu,Bd);
% velocity shear angle
info.vsh = shear_angle(Vu,Vd);
% compression factors to the models
info.sig = sig;
% other paramters from the models
if isfield(spec,'R')
  info.eps = eps;
  info.L = L;
  info.alpha = alpha;
  info.x0 = x0;
  info.y0 = y0;
end
% constraint matrix
info.cmat = constraint_values(spec,n);

% gather data
nd.info = info;
nd.n = n;
nd.thBn = thBn;
nd.thVn = thVn;
nd.Vsh = Vsh;

if nargout == 0
  for ii = 1:length(fnames)
    fprintf(['n_',fnames{ii},'\t = \t',num2str(n.(fnames{ii})),'\n'])
  end
end


end


function th = shear_angle(Au,Ad)

th = acosd(dot(Au,Ad)/(norm(Au)*norm(Ad)));

end

function th = shock_angle(spec,n,field,leq90)
% field is 'B' or 'V'

switch lower(field)
  case 'b'
    a = spec.Bu;
  case 'v'
    a = spec.Vu;
end

fnames = fieldnames(n);
num = length(fnames);

for i = 1:num
  th.(fnames{i}) = thFn(n.(fnames{i}),a,leq90);
end

end

function th = thFn(nvec,a,leq90)
% Shock normal angle and normal incidence angle
th = acosd(dot(a,nvec)/(norm(a)));
if th>90 && leq90
  th = 180-th;
end
end


function cmat = constraint_values(spec,n)

fnames = fieldnames(n);
num = length(fnames);
cmat = zeros(5,num);

fun1 = @(nvec)dot(spec.delB,nvec)/norm(spec.delB); %#ok<NASGU>
fun2 = @(nvec)dot(cross(spec.Bd,spec.Bu),nvec)/norm(cross(spec.Bd,spec.Bu)); %#ok<NASGU>
fun3 = @(nvec)dot(cross(spec.Bu,spec.delV),nvec)/norm(cross(spec.Bu,spec.delV)); %#ok<NASGU>
fun4 = @(nvec)dot(cross(spec.Bd,spec.delV),nvec)/norm(cross(spec.Bd,spec.delV)); %#ok<NASGU>
fun5 = @(nvec)dot(cross(spec.delB,spec.delV),nvec)/norm(cross(spec.delB,spec.delV)); %#ok<NASGU>

c_eval('cmat(?,:) = structfun(fun?,n);',1:5)

% fields that are by definition zero are set to 0?
idz = sub2ind(size(cmat),[1,1,1,1,2,3,3,4,4,5,5],[1,3,4,5,1,2,3,2,4,2,5]);
cmat(idz) = 0;
end


function Vsp = shock_speed(spec,n,thBn,method)

fn = fieldnames(n);
N = length(fn);


switch lower(method)
  case 'gt' % Gosling & Thomsen
    if ~isfield(spec,'Fcp') || ~isfield(spec,'dTf') || ~isfield(spec,'d2u')
      for k = 1:N
        Vsp.(fn{k}) = 0;
      end
      return;
    else
      Vsp = speed_gosling_thomsen(spec,n,thBn,fn);
    end
  case 'mf' % Mass flux
    Vsp = speed_mass_flux(spec,n,fn);
  case 'sb' % Smith & Burton
    Vsp = speed_smith_burton(spec,n,fn);
  case 'mo' % Moses
    if ~isfield(spec,'Fcp') || ~isfield(spec,'dTf') || ~isfield(spec,'d2u')
      for k = 1:N
        Vsp.(fn{k}) = 0;
      end
      return;
    else
      Vsp = speed_moses(spec,n,thBn,fn);
    end
end


end

function Vsp = speed_gosling_thomsen(spec,n,thBn,fn)
for k = 1:length(fn)
  th = thBn.(fn{k})*pi/180;
  nvec = n.(fn{k});

  % Notation as in (Gosling and Thomsen 1985)
  W = spec.Fcp*2*pi;
  t1 = acos((1-2*cos(th).^2)./(2*sin(th).^2))/W;

  f = @(th)W*t1*(2*cos(th).^2-1)+2*sin(th).^2.*sin(W*t1);
  x0 = f(th)/(W*spec.dTf);

  % the sign of Vsh in this method is ambiguous, assume n points upstream
  Vsp.(fn{k}) = spec.d2u*dot(spec.Vu,nvec)*(x0/(1+spec.d2u*x0));
end
end

function Vsp = speed_mass_flux(spec,n,fn)
% Assume all protons, not very good but no composition available
u = irf_units;

rho_u = spec.nu*u.mp;
rho_d = spec.nd*u.mp;

for k = 1:1:length(fn)
  Vsp.(fn{k}) = (rho_d*dot(spec.Vd,n.(fn{k}))-rho_u*dot(spec.Vu,n.(fn{k})))/(rho_d-rho_u);
end
end

function Vsp = speed_smith_burton(spec,n,fn)
% Assuming n points upstream, there is only one solution to the equation
for k = 1:1:length(fn)
  Vsp.(fn{k}) = dot(spec.Vu,n.(fn{k}))+norm(cross((spec.Vd-spec.Vu),spec.Bd))/norm(spec.Bd-spec.Bu);
end

end

function Vsp = speed_moses(spec,n,thBn,fn)
for k = 1:length(fn)
  th = thBn.(fn{k})*pi/180;
  nvec = n.(fn{k});
  thVn = acos(dot(nvec,spec.Vu)/norm(spec.Vu));

  % Notation as in (Moses et al., 1985)
  W = spec.Fcp*2*pi;
  x = 0.68*sin(th)^2*cos(thVn)/(W*spec.dTf);

  Vsp.(fn{k}) = dot(spec.Vu,nvec)*(x/(1+spec.d2u*x));
end
end


function [n,sig0,eps,L,alpha,x0,y0] = farris_model(spec)
eps = 0.81;
L = 24.8; % in RE
x0 = 0;
y0 = 0;

alpha = 3.8;

[n,sig0] = shock_model(spec,eps,L,x0,y0,alpha);
end

function [n,sig0,eps,L,alpha,x0,y0] = slavin_holzer_model(spec) %
eps = 1.16;
L = 23.3; % in RE
x0 = 3.0;
y0 = 0;
alpha = -atand(spec.Vu(2)/spec.Vu(1));

[n,sig0] = shock_model(spec,eps,L,x0,y0,alpha);
end

function [n,sig0,eps,L,alpha,x0,y0] = peredo_model(spec) %
eps = 0.98;
L = 26.1; % in RE
x0 = 2.0;
y0 = 0.3;
alpha = 3.8-0.6;

[n,sig0] = shock_model(spec,eps,L,x0,y0,alpha);
end

function [n,sig0,eps,L,alpha,x0,y0] = fairfield_meridian_4o_model(spec) %
eps = 1.02;
L = 22.3; % in RE
x0 = 3.4;
y0 = 0.3;
alpha = 4.8;

[n,sig0] = shock_model(spec,eps,L,x0,y0,alpha);
end

function [n,sig0,eps,L,alpha,x0,y0] = fairfield_meridian_no_4o_model(spec) %
eps = 1.05;
L = 20.5; % in RE
x0 = 4.6;
y0 = 0.4;
alpha = 5.2;

[n,sig0] = shock_model(spec,eps,L,x0,y0,alpha);
end

function [n,sig0,eps,L,alpha,x0,y0] = formisano_unnorm_model(spec) %
eps = 0.97;
L = 22.8; % in RE
x0 = 2.6;
y0 = 1.1;
alpha = 3.6;

[n,sig0] = shock_model(spec,eps,L,x0,y0,alpha);
end


function [n,sig0] = shock_model(spec,eps,L,x0,y0,alpha) % Method from ISSI book
u = irf_units;

% rotation matrix
R = [cosd(alpha),-sind(alpha),0;sind(alpha),cosd(alpha),0;0,0,1];

% offset from GSE
r0 = [x0;y0;0];

% sc position in GSE (or GSM or whatever) in Earth radii
if isstruct(spec.R) % MMS specific
  % like mms.get_data returns
  rsc_sum = 0;
  if isfield(spec,'t') % get average position of sc at given time t
    c_eval('rsc_sum = rsc_sum+interp1(spec.R.time.epochUnix,spec.R.gseR?/(u.RE*1e-3),spec.t.epochUnix);')
    rsc = rsc_sum'/4;
  else % get time and spacecraft average
    rsc_sum = 0;
    c_eval('rsc_sum = rsc_sum+mean(spec.R.gseR?)/(u.RE*1e-3);')
    rsc = rsc_sum'/4;
  end
elseif isa(spec.R,'TSeries') % TSeries

  rsc = mean(spec.R.data)'/(u.RE*1e-3);
elseif isnumeric(spec.R) && length(R) == 3 % just a vector
  rsc = spec.R;
end

% Calculate sigma
% sc position in the natural system (cartesian)
rp = @(sig)R*rsc-sig.*r0;
% sc polar angle in the natural system
thp = @(sig)cart2pol([1,0,0]*rp(sig),sqrt(([0,1,0]*rp(sig))^2+([0,0,1]*rp(sig))^2)); % returns angle first
% minimize |LH-RH| in eq 10.22
fval = @(sig)abs(sig.*L./sqrt(sum(rp(sig).^2,1))-1-eps*cos(thp(sig)));

% find the best fit for sigma
sig0 = fminsearch(fval,1);
% to make sure it finds the largest sigma
sig0 = fminsearch(fval,2*sig0);

% calculate normal
xp = [1,0,0]*rp(sig0);
yp = [0,1,0]*rp(sig0);
zp = [0,0,1]*rp(sig0);

% gradient to model surface
gradS = [(xp*(1-eps^2)+eps*sig0*L)*cosd(alpha)+yp*sind(alpha);...
  -(xp*(1-eps^2)+eps*sig0*L)*sind(alpha)+yp*cosd(alpha);...
  zp]/(norm(rp(sig0))*2*sig0*L);
% normal vector (column vector)
n = gradS'/norm(gradS);
end
