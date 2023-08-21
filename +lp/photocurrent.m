function [jPhoto] = photocurrent( iluminatedArea, U, distanceSunAU ,flag)
% LP.PHOTOCURRENT Langmuir probe photocurrent for different materials
%
% j_photo = LP.PHOTOCURRENT( iluminatedArea, U, distanceSunAU )
%
%   Calculates the photo-current emitted by an arbitrary body with
%   cross section area, iluminatedArea [m^2], at a distance distanceSunAU [AU]
%   from the Sun, when the body has a potential, U [V]. Estimates
%   are done for the solar minimum conditions.
%
%   iluminatedArea - scalar
%   U              - scalar or vector or matrix
%
% j_photo = LP.PHOTOCURRENT( iluminatedArea, U, distanceSunAU, surfaceMaterial )
%       returns photo current for the specified material
%
%   some of implemented surface materials
%       'cluster' - same as 'default', based on Pedersen papers
%       'themis'  - based on corrected THEMIS sweeps
%       'cassini' - saturation current based on experimental data
%   for a full list execute LP.PHOTOCURRENT without arguments
%
% j_photo = LP.PHOTOCURRENT( iluminatedArea, U, distanceSunAU, surfacePhotoemission )
%   surfacePhotoemission should be given at 1AU distance in units [A/m^2]
%
% surface_materials=LP.PHOTOCURRENT
%   returns cell array with implemented surface materials
%

% Check # input/output parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narginchk(0,4)

surface_materials={'default','cluster','themis','cassini','aluminium','aquadag','gold','graphite','solar cells','1eV','TiN','elgiloy'};
if nargin==0 && nargout == 0
  for ii=1:numel(surface_materials)
    surf=surface_materials{ii};
    j0=lp.photocurrent(1,0,1,surf);
    disp([surf ': Io= ' num2str(j0*1e6,2) ' uA/m2']);
  end
  return
end
if nargin==0 && nargout == 1
  jPhoto=surface_materials;
  return
end
if nargin==3, flag='default';end
if nargin == 4 && isnumeric(flag) % specified photemission
  photoemission = flag;
  flag = 'photoemission given';
end
switch lower(flag)
  case 'photoemission given'

    jPhoto         = zeros(size(U)); % initialize
    jPhoto(:)      = photoemission*iluminatedArea/distanceSunAU^2; % initialize to current valid for negative potentials
    jPhoto(U >= 0) = photoemission*(iluminatedArea/distanceSunAU^2) .* ...
      ( 5.0e-5/5.6e-5 .* exp( - U(U>=0) ./ 2.74 ) ...
      + 1.2e-5/5.6e-5 .* exp( - (U(U>=0) + 10.0) ./ 14.427 ) );

  case 'default'
    % The Photo-current emitted depends on if the potential of the body is
    % positive or negative. Reference needed TODO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    jPhoto         = zeros(size(U)); % initialize
    jPhoto(:)      = 5.6e-5*iluminatedArea/distanceSunAU^2; % initialize to current valid for negative potentials
    jPhoto(U >= 0) = (iluminatedArea/distanceSunAU^2) .* ...
      ( 5.0e-5 .* exp( - U(U>=0) ./ 2.74 ) ...
      + 1.2e-5 .* exp( - (U(U>=0) + 10.0) ./ 14.427 ) );

  case '1ev'

    jPhoto         = zeros(size(U));
    jPhoto(:)      = 5.0e-5*iluminatedArea/distanceSunAU^2; % initialize to current valid for negative potentials
    jPhoto(U >= 0) = (iluminatedArea/distanceSunAU^2) .* ...
      ( 5.0e-5 .* exp( - U(U>=0) ) );

  case 'themis'
    refU      = [.1 1   5  10  50];
    refJPhoto = [50 27 10   5   .5]*1e-6;
    logU      = log(refU);
    logJ      = log(refJPhoto);

    jPhoto     = ones(size(U)); %
    jPhoto     = jPhoto*refJPhoto(1)*iluminatedArea/distanceSunAU^2; % negative potentials
    ii         = find( U >= refU(1) );
    jPhoto(ii) = exp(interp1(logU,logJ,log(U(ii)),'pchip','extrap'))*iluminatedArea/distanceSunAU^2;

  case {'cassini','tin'}
    jZero       = 25e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case 'cluster'
    % Cluster is like aquadag but closer to alluminium, we take 25 uA/m2 at 1AU
    jZero       = 25e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case 'aluminium'
    %aluminium is 30 uA/m2 at 1AU (roughly from Erik Winkler exjobb)
    jZero       = 30e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case 'aquadag'
    jZero       = 18e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case 'gold'
    jZero       = 29e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case 'graphite'
    jZero       = 7.2e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case {'solar cells','solar cell'}
    jZero       = 20e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  case 'elgiloy'
    jZero       = 30e-6;
    jZeroThemis = lp.photocurrent(1,0,1,'themis');
    jPhoto      = jZero/jZeroThemis ...
      * lp.photocurrent(iluminatedArea,U,distanceSunAU,'themis');

  otherwise
end
