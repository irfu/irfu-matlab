%-Abstract
%
%   CSPICE_CONICS determines the state (position, velocity) of an orbiting
%   body from a set of elliptic, hyperbolic, or parabolic orbital elements.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   "AS-IS" TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
%   ANY WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR
%   A PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY,
%   OR NASA BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING,
%   BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
%   ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY
%   AND LOST PROFITS, REGARDLESS OF WHETHER CALTECH, JPL, OR
%   NASA BE ADVISED, HAVE REASON TO KNOW, OR, IN FACT, SHALL
%   KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE
%   OF THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO
%   INDEMNIFY CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING
%   FROM THE ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      elts     the array(s) containing the conic osculating elements
%               describing the orbit of a body around a primary.
%
%               [8,n] = size(elts); double = class(elts)
%
%               The elements are, in order:
%
%                  RP      Perifocal distance.
%                  ECC     Eccentricity.
%                  INC     Inclination.
%                  LNODE   Longitude of the ascending node.
%                  ARGP    Argument of periapse.
%                  M0      Mean anomaly at epoch.
%                  T0      Epoch.
%                  MU      Gravitational parameter.
%
%               Units are km, rad, rad/sec, km**3/sec**2.
%
%               The epoch T0 is given in ephemeris seconds past J2000.
%               T0 is the instant at which the state of the body is
%               specified by the elements.
%
%      et       the ephemeris time(s) corresponding one-to-one and onto
%               to each `elts' at which to determine the state of
%               the orbiting body
%
%               [1,n] = size(et); double = class(et)
%
%               Note: The design of cspice_conics supposes the inputs `elts'
%               and `et' originates as the output of another Mice routine
%               and so will have the same vectorization measure.
%
%               Still, in the event the user requires an `elts' constant over
%               a vector of `et', or an `et' constant over an array of
%               `elts', construct the needed variables with the Matlab code:
%
%                  Given a constant `epoch' for an array of `elts', create the
%                  vector `et'.
%
%                     N          = size(elts,2);
%                     et         = zeros(1, N) + epoch;
%
%                  Given a constant element set `elt' for an array of `et',
%                  create the array `elts'.
%
%                     N          = size(et,1);
%                     elts       = zeros(8, N);
%                     elts(1,:)  = elt(1);
%                     elts(2,:)  = elt(2);
%                     elts(3,:)  = elt(3);
%                     elts(4,:)  = elt(4);
%                     elts(5,:)  = elt(5);
%                     elts(6,:)  = elt(6);
%                     elts(7,:)  = elt(7);
%                     elts(8,:)  = elt(8);
%
%   the call:
%
%      [state] = cspice_conics( elts, et )
%
%   returns
%
%      state    the array(s) representing the state (position and velocity) of
%               the body at time `et' in kilometers and kilometers-per-second
%               (the first three components of `state' represent the x-,
%               y-, and z-components of the body's position; the last three
%               components form the corresponding velocity vector)
%
%               [6,n] = size(state); double = class(state)
%
%               `state' returns with the same vectorization measure, N, as
%               `elts' and `et'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Calculate the perturbation between the
%      state elements of the Moon at some time as determined
%      from SPK data and the corresponding state elements
%      determined from propagation of osculating elements.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: conics_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            gm_de431.tpc                  Gravitational constants
%            naif0012.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
%                                'gm_de431.tpc',
%                                'naif0012.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function conics_ex1()
%
%         %
%         % Load the meta kernel listing the needed SPK, PCK, LSK
%         % kernels, and a PCK kernel that contains gravitation constants.
%         %
%         cspice_furnsh( 'conics_ex1.tm' )
%
%         %
%         % Convert the time of interest, provided as a string, to ephemeris
%         % time.
%         %
%         et = cspice_str2et( 'Dec 25, 2007' );
%
%         %
%         % Call cspice_spkezr to retrieve the Moon state
%         % w.r.t. the earth in the 'J2000' frame. Use 'NONE' as aberration
%         % correction since we are comparing geometric states.
%         %
%         [state, lt] = cspice_spkezr( 'Moon',  et,               ...
%                                      'J2000', 'NONE', 'EARTH' );
%
%         %
%         % Read the gravitational parameter for Earth.
%         %
%         mu = cspice_bodvrd( 'EARTH', 'GM', 1 );
%
%         %
%         % Execute the cspice_oscelt call to convert the state 6-vector
%         % to the osculating elements 8-vector, `elts', at `et'. The
%         % osculating elements are relative to the same frame as `state'.
%         %
%         % The elements describe the nominal orbit the Moon would follow
%         % in the absence of all other bodies in the solar system and
%         % and all non-gravitational forces.
%         %
%         % Note: the cspice_bodvrd call returns data as arrays, so
%         % to access the gravitational parameter (the only value in
%         % the array), we use 'mu(1)'.
%         %
%         elts = cspice_oscelt( state, et, mu(1) );
%
%         %
%         % Now, select a time one week from the initial epoch.
%         %
%         later = et + 7. * cspice_spd;
%
%         %
%         % Use the osculating elements to calculate the state vector
%         % of the Moon at the 'later' epoch.
%         %
%         later_state = cspice_conics( elts, later );
%
%         %
%         % Now retrieve the Moon's state at time 'later' from SPK
%         % data.
%         %
%         [state, lt] = cspice_spkezr('Moon',  later,           ...
%                                     'J2000', 'NONE', 'EARTH');
%
%         %
%         % Display the absolute diff between the vector output by
%         % cspice_conics and the state vector returned by cspice_spkezr.
%         %
%         pert = later_state - state;
%
%         txt = sprintf( 'Perturbation in     x: %16.8f', pert(1) );
%         disp( txt )
%
%         txt = sprintf( 'Perturbation in     y: %16.8f', pert(2) );
%         disp( txt )
%
%         txt = sprintf( 'Perturbation in     z: %16.8f', pert(3) );
%         disp( txt )
%
%         txt = sprintf( 'Perturbation in dx/dt: %16.8f', pert(4) );
%         disp( txt )
%
%         txt = sprintf( 'Perturbation in dy/dt: %16.8f', pert(5) );
%         disp( txt )
%
%         txt = sprintf( 'Perturbation in dz/dt: %16.8f', pert(6) );
%         disp( txt )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Perturbation in     x:   -7488.85977321
%      Perturbation in     y:     397.61007948
%      Perturbation in     z:     195.74558097
%      Perturbation in dx/dt:      -0.03615276
%      Perturbation in dy/dt:      -0.00127927
%      Perturbation in dz/dt:      -0.00201459
%
%
%   2) Calculate the magnitude of the perturbation between the
%      position and velocity vectors of the Moon w.r.t. earth as
%      calculated from cspice_conics and as retrieved from an SPK file.
%
%      Use the meta-kernel from the first example.
%
%
%      Example code begins here.
%
%
%      function conics_ex2()
%
%         %
%         % Load the meta kernel listing the needed SPK, PCK, LSK
%         % kernels.
%         %
%         cspice_furnsh( 'conics_ex1.tm' )
%
%         %
%         % Convert the time of interest, provided as a string, to ephemeris
%         % time.
%         %
%         et1 = cspice_str2et( 'Jan 1, 2007' );
%
%         %
%         % Make the cspice_spkezr call to retrieve the state of the
%         % Moon w.r.t. the earth in J2000. Use 'NONE' as aberration
%         % correction since we are comparing geometric states.
%         %
%         [state1, lt] = cspice_spkezr( 'Moon',  et1,              ...
%                                       'J2000', 'NONE', 'EARTH' );
%
%         %
%         % Read the gravitational parameter for Earth.
%         %
%         mu    = cspice_bodvrd( 'EARTH', 'GM', 1 );
%
%         elts1 = cspice_oscelt( state1, et1, mu(1) );
%
%         %
%         % We want to propagate the osculating elements in 'elts1'
%         % by N time steps. Create an array of N copies of 'elts1'
%         %
%         N     = 15;
%         elts  = repmat( elts1, 1, N );
%
%         %
%         % Create an array of N ephemeris times in steps of one day
%         % (measured in seconds) from `et1'.
%         %
%         et             = [1:N]*cspice_spd + et1;
%
%         twobody        = cspice_conics( elts, et );
%         [state, lt] = cspice_spkezr( 'Moon', et, 'J2000', 'NONE', 'EARTH' );
%         utc            = cspice_et2utc( et, 'C', 0 );
%
%         for n=1:N
%            txt = sprintf(                                       ...
%                   '%s perturbation: ||r|| %10.4f, ||v|| %6.4f', ...
%                    utc(n,:)                                   , ...
%                    norm( state(1:3,n) - twobody(1:3,n) )      , ...
%                    norm( state(4:6,n) - twobody(4:6,n) )            );
%            disp( txt )
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      2007 JAN 02 00:00:00 perturbation: ||r||    91.3141, ||v|| 0.0020
%      2007 JAN 03 00:00:00 perturbation: ||r||   312.2194, ||v|| 0.0030
%      2007 JAN 04 00:00:00 perturbation: ||r||   574.8464, ||v|| 0.0030
%      2007 JAN 05 00:00:00 perturbation: ||r||   789.2552, ||v|| 0.0021
%      2007 JAN 06 00:00:00 perturbation: ||r||   880.3755, ||v|| 0.0014
%      2007 JAN 07 00:00:00 perturbation: ||r||   808.2985, ||v|| 0.0033
%      2007 JAN 08 00:00:00 perturbation: ||r||   628.4228, ||v|| 0.0061
%      2007 JAN 09 00:00:00 perturbation: ||r||   760.3389, ||v|| 0.0096
%      2007 JAN 10 00:00:00 perturbation: ||r||  1581.8352, ||v|| 0.0141
%      2007 JAN 11 00:00:00 perturbation: ||r||  2978.3503, ||v|| 0.0202
%      2007 JAN 12 00:00:00 perturbation: ||r||  5011.5684, ||v|| 0.0282
%      2007 JAN 13 00:00:00 perturbation: ||r||  7828.8170, ||v|| 0.0381
%      2007 JAN 14 00:00:00 perturbation: ||r|| 11573.3980, ||v|| 0.0498
%      2007 JAN 15 00:00:00 perturbation: ||r|| 16336.6354, ||v|| 0.0628
%      2007 JAN 16 00:00:00 perturbation: ||r|| 22123.7052, ||v|| 0.0765
%
%
%-Particulars
%
%   None.
%
%-Exceptions
%
%   1)  If the eccentricity supplied is less than 0, the error
%       SPICE(BADECCENTRICITY) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If a non-positive periapse distance is supplied, the error
%       SPICE(BADPERIAPSEVALUE) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If a non-positive value for the attracting mass is supplied,
%       the error SPICE(BADGM) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If `elts' is such that the resulting orbit at periapsis has
%       either its position or velocity equal to zero, or the square
%       of the resulting specific angular momentum's magnitude is
%       zero, an error is signaled by a routine in the call tree of
%       this routine. This is an indication of invalid `elts' elements.
%
%   5)  If `et' is such that the offset in time from periapsis, at which
%       the state is to be determined, is so large that there is a
%       danger of floating point overflow during computation, an error
%       is signaled by a routine in the call tree of this routine.
%
%   6)  If any of the input arguments, `elts' or `et', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   7)  If any of the input arguments, `elts' or `et', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   8)  If the input vectorizable arguments `elts' and `et' do not
%       have the same measure of vectorization (N), an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  R. Bate, D. Mueller, and J. White, "Fundamentals of
%        Astrodynamics," Dover Publications Inc., 1971.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 23-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Reduced
%       number of time steps used in code example #2. Added a call to
%       cspice_kclear in code example #1. Added meta-kernel to example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 30-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Added to -I/O section a description of creating vectorized variables
%       from constant values, i.e. create a vectorized 'et' from a constant
%       (non vectorized) epoch, or create a vectorized 'elts' from a
%       constant (non vectorized) single set of elements.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   state from conic elements
%
%-&

function [state] = cspice_conics( elts, et )

   switch nargin
      case 2

         elts = zzmice_dp(elts);
         et   = zzmice_dp(et);

      otherwise

         error ( 'Usage: [_state(6)_] = cspice_conics( _elts(8)_, _et_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [state] = mice('conics_c', elts, et);
   catch spiceerr
      rethrow(spiceerr)
   end



