function [j_photo] = lp_photocurrent( X_area, U_pot, R_sun )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [j_photo] = lp_photocurrent( X_area, U_pot, R_sun )
%
%   Matlab function that calculates the photo-current emitted by an
%   arbitrary body with cross area, X_area [m^2], a distance, R_sun [AU], 
%   from the Sun, when the body has a potential, U_pot [V]. 
%
%   Input parameters:    X_area             --- scalar
%                        U_pot              --- scalar or vector
%
%   Created by Jan-Erik Wahlund, Cornell University, October-1994.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check # input/output parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  error(nargchk(3,3,nargin))
%  error(nargchk(1,1,nargout))

% Check size of V_pot, and make it an column vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  U_pot = U_pot(:);
  n_pts = length( U_pot );

% The Photo-current emitted depends on if the potential of the body is
% positive or negative.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pos_ind = find( U_pot >= 0 );
  neg_ind = find( U_pot < 0 );

  j_photo = zeros( n_pts, 1 );

  j_photo(pos_ind) = (X_area/R_sun^2) .* ...
                     ( 5.0e-5 .* exp( - U_pot(pos_ind) ./ 2.74 ) ...
                     + 1.2e-5 .* exp( - (U_pot(pos_ind) + 10.0) ./ 14.427 ) );

  j_photo(neg_ind) = (5.6e-5*X_area/R_sun^2) .* ones(length(neg_ind),1);

  return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

