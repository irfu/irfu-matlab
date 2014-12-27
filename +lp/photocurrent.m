function [j_photo] = photocurrent( X_area, U_pot, R_sun ,flag)
% LP.PHOTOCURRENT Langmuir probe photocurrent for different materials
%
% j_photo = LP.PHOTOCURRENT( X_area, U_pot, R_sun )
%
%   Calculates the photo-current emitted by an arbitrary body with 
%   cross section area, X_area [m^2], at a distance R_sun [AU]
%   from the Sun, when the body has a potential, U_pot [V]. Estimates
%   are done for the solar minimum conditions. 
%
%   X_area  - scalar
%   U_pot   - scalar or vector or matrix
%
% j_photo = LP.PHOTOCURRENT( X_area, U_pot, R_sun, surface_material )
%       returns photo current for specified material
% 
%   some of implemented surface materials
%   for full list execute LP_PHOTOCURRENT without arguments
%       'cluster' - same as 'default', based on Pedersen papers
%       'themis'  - based on corrected THEMIS sweeps
%       'cassini' - saturation current based on experimental data
% 
% surface_materials=LP_PHOTOCURRENT
%   returns cell array with implemented surface materials
%

% Check # input/output parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narginchk(0,4)

surface_materials={'default','cluster','themis','cassini','aluminium','aquadag','gold','graphite','solar cells','1eV'};
if nargin==0 && nargout ==0, 
    for ii=1:numel(surface_materials)
        surf=surface_materials{ii};
        j0=lp.photocurrent(1,0,1,surf);
        disp([surf ': Io= ' num2str(j0*1e6,2) ' uA/m2']);
    end
    return
end
if nargin==0 && nargout == 1,
    j_photo=surface_materials;
    return
end
if nargin==3, flag='default';end

% Check size of V_pot, and make it an column vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(flag)
    case 'default'
        % The Photo-current emitted depends on if the potential of the body is
        % positive or negative.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        j_photo = zeros(size(U_pot)); % initialize
        j_photo(:) = 5.6e-5*X_area/R_sun^2; % initialize to current valid for negative potentials
        j_photo(U_pot >= 0) = (X_area/R_sun^2) .* ...
            ( 5.0e-5 .* exp( - U_pot(U_pot>=0) ./ 2.74 ) ...
            + 1.2e-5 .* exp( - (U_pot(U_pot>=0) + 10.0) ./ 14.427 ) );
        
    case '1ev'

        j_photo = zeros(size(U_pot)); % initialize
        j_photo(:) = 5.0e-5*X_area/R_sun^2; % initialize to current valid for negative potentials
        j_photo(U_pot >= 0) = (X_area/R_sun^2) .* ...
            ( 5.0e-5 .* exp( - U_pot(U_pot>=0) ./ 1.0 ) );
        
    case 'themis'
        U_ref  =        [.1 1   5  10  50];
        j_photo_ref   = [50 27 10   5   .5]*1e-6;
       
        logU=log(U_ref);
        logj=log(j_photo_ref);
        
        j_photo = ones(size(U_pot)); %
        j_photo = j_photo*j_photo_ref(1)*X_area/R_sun^2; % negative potentials
        pos_ind = find( U_pot >= U_ref(1) );        
        j_photo(pos_ind) = exp(interp1(logU,logj,log(U_pot(pos_ind)),'pchip','extrap'))*X_area/R_sun^2;
    case 'cassini'        
        % Cassini material is 25 uA/m2 at 1AU
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=25e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');
    case 'cluster'        
        % Cluster is like aquadag but closer to alluminium, we take 25 uA/m2 at 1AU
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=25e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');

    case 'aluminium'        
        %aluminium is 30 uA/m2 at 1AU (roughly from Erik Winkler exjobb)
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=30e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');
        
    case 'aquadag'        
        %aluminium is 40 uA/m2 at 1AU, we scale THEMIS 
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=18e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');
        
    case 'gold'        
        %gold is 29 uA/m2 at 1AU, we scale THEMIS 
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=29e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');
        
    case 'graphite'        
        %graphite is 7.2 uA/m2 at 1AU, we scale THEMIS 
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=7.2e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');
        
    case {'solar cells','solar cell'}        
        %graphite is 20 uA/m2 at 1AU, we scale THEMIS 
        j0=lp.photocurrent(1,0,1,'themis');
        scale_factor=20e-6/j0;
        j_photo=scale_factor*lp.photocurrent(X_area,U_pot,R_sun,'themis');
        
    otherwise
end

