%-Abstract
%
%   CSPICE_OSCELT calculates the set of osculating conic
%   orbital elements corresponding to the state 6-vector
%   (position, velocity) of a body at an epoch.
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
%      state    the state(s) (position and velocity) of the body at some
%               epoch.
%
%               [6,n] = size(state); double = class(state)
%
%               Components are x, y, z, dx/dt, dy/dt, dz/dt. `state' must
%               be expressed relative to an inertial reference frame. Units
%               are km and km/sec.
%
%      et       the epoch(s) of the input state, in ephemeris seconds past
%               J2000.
%
%               [1,n] = size(et); double = class(et)
%
%                                                  3    2
%      mu       the gravitational parameter (gm, km /sec ) of the primary
%               body.
%
%               [1,1] = size(mu); double = class(mu)
%
%   the call:
%
%      [elts] = cspice_oscelt( state, et, mu )
%
%   returns:
%
%      elts     equivalent conic elements describing the orbit of the body
%               around its primary.
%
%               [8,n] = size(elts); double = class(elts)
%
%               The elements are, in order:
%
%                  RP      Perifocal distance.
%                  ECC     Eccentricity.
%                  INC     Inclination.
%                  LNODE   Longitude of the ascending node.
%                  ARGP    Argument of periapsis.
%                  M0      Mean anomaly at epoch.
%                  T0      Epoch.
%                  MU      Gravitational parameter.
%
%               The epoch of the elements is the epoch of the input
%               state. Units are km, rad, rad/sec. The same elements
%               are used to describe all three types (elliptic,
%               hyperbolic, and parabolic) of conic orbit.
%
%               `elts' returns with the same vectorization measure, N,
%               as `state' and `et'.
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
%   1) Determine the osculating elements of Phobos with respect to
%      Mars at some arbitrary time in the J2000 inertial reference
%      frame.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: oscelt_ex1.tm
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
%            mar097.bsp                    Mars satellite ephemeris
%            gm_de431.tpc                  Gravitational constants
%            naif0012.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'mar097.bsp',
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
%      function oscelt_ex1()
%
%         %
%         % Load the meta kernel listing the needed SPK, LSK and
%         % PCK with gravitational parameters kernels.
%         %
%         cspice_furnsh( 'oscelt_ex1.tm' )
%
%         %
%         % Convert the time string to ephemeris time
%         %
%         et = cspice_str2et( 'Dec 25, 2007' );
%
%         %
%         % Make the cspice_spkezr call to retrieve the state of
%         % Phobos with respect to Mars in J2000.
%         %
%         [state, lt] = cspice_spkezr( 'PHOBOS', et, 'J2000', ...
%                                      'NONE', 'MARS'         );
%
%         %
%         % Read the gravitational parameter for Mars.
%         %
%         mu = cspice_bodvrd( 'MARS', 'GM', 1 );
%
%         %
%         % make the cspice_oscelt call to convert the state 6-vector
%         % to the elts 8-vector. Note: the  cspice_bodvrd returns
%         % data as arrays, so to access the gravitational parameter
%         % (the only value in the array), we use mu(1).
%         %
%         elts = cspice_oscelt( state, et, mu(1) );
%
%         %
%         % Output the elts vector.
%         %
%         fprintf( 'Perifocal distance          (km): %21.10f\n', elts(1))
%         fprintf( 'Eccentricity                    : %21.10f\n', elts(2))
%         fprintf( 'Inclination                (deg): %21.10f\n',      ...
%                                                   elts(3) * cspice_dpr )
%         fprintf( 'Lon of ascending node      (deg): %21.10f\n',      ...
%                                                   elts(4) * cspice_dpr )
%         fprintf( 'Argument of periapsis      (deg): %21.10f\n',      ...
%                                                   elts(5) * cspice_dpr )
%         fprintf( 'Mean anomaly at epoch      (deg): %21.10f\n',      ...
%                                                   elts(6) * cspice_dpr )
%         fprintf( 'Epoch                        (s): %21.10f\n', elts(7))
%         fprintf( 'Gravitational parameter (km3/s2): %21.10f\n', elts(8))
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
%      Perifocal distance          (km):       9232.5746716211
%      Eccentricity                    :          0.0156113904
%      Inclination                (deg):         38.1225231660
%      Lon of ascending node      (deg):         47.0384055902
%      Argument of periapsis      (deg):        214.1546430017
%      Mean anomaly at epoch      (deg):        340.5048466068
%      Epoch                        (s):  251812865.1837092042
%      Gravitational parameter (km3/s2):      42828.3736206991
%
%
%   2) Calculate the history of Phobos's orbit plane inclination
%      with respect to Mars in the J2000 frame at intervals of six
%      months for a time interval of 10 years.
%
%      Use the meta-kernel from the first example.
%
%
%      Example code begins here.
%
%
%      function oscelt_ex2()
%
%         %
%         % Load the meta kernel listing the needed SPK, LSK and
%         % PCK with gravitational parameters kernels.
%         %
%         cspice_furnsh( 'oscelt_ex1.tm' )
%
%         %
%         % Read the gravitational parameter for Mars.
%         %
%         mu = cspice_bodvrd( 'MARS', 'GM', 1 );
%
%         %
%         % The start epoch.
%         %
%         et0 = cspice_str2et( 'Jan 1, 2000 12:00:00' );
%
%         %
%         % A step of six months - in seconds.
%         %
%         step = 180. * cspice_spd;
%
%         %
%         % Define an array of ephemeris times, covering,
%         % 10 years in steps of six months starting
%         % approximately Jan 1, 2000.
%         %
%         et = [0: 19]*step + et0;
%
%         %
%         % Retrieve the state; convert to osculating elements.
%         %
%         [state,lt] = cspice_spkezr( 'PHOBOS', et, 'J2000', ...
%                                     'NONE', 'MARS'        );
%         elts       = cspice_oscelt( state, et, mu(1) );
%
%         %
%         % Convert the angular measures to degrees.
%         %
%         elts(3,:) = [ elts(3,:) * cspice_dpr ];
%
%         %
%         % Convert the ephemeris time of the state lookup to
%         % calendar UTC, then print the calendar string and the
%         % inclination in degrees of Phobos wrt Mars at the
%         % time.
%         %
%         utcstr = cspice_et2utc( et, 'C', 3 );
%
%         %
%         % Output the epoch and corresponding inclination.
%         %
%         disp( '        UCT Time          Inclination' )
%         disp( '------------------------  -----------' )
%         for n=1:20
%            fprintf( '%s %12.6f\n', utcstr(n,:), elts(3,n) );
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
%              UCT Time          Inclination
%      ------------------------  -----------
%      2000 JAN 01 12:00:00.000    36.055248
%      2000 JUN 29 12:00:00.000    37.112144
%      2000 DEC 26 12:00:00.000    38.152129
%      2001 JUN 24 12:00:00.000    37.552071
%      2001 DEC 21 12:00:00.000    36.242049
%      2002 JUN 19 11:59:59.999    36.330470
%      2002 DEC 16 12:00:00.000    37.674595
%      2003 JUN 14 11:59:59.999    38.121191
%      2003 DEC 11 12:00:00.001    36.973204
%      2004 JUN 08 11:59:59.999    36.033732
%      2004 DEC 05 12:00:00.001    36.844542
%      2005 JUN 03 11:59:59.999    38.077365
%      2005 NOV 30 12:00:00.001    37.786106
%      2006 MAY 29 11:59:58.999    36.413540
%      2006 NOV 25 11:59:59.001    36.171050
%      2007 MAY 24 11:59:58.999    37.448015
%      2007 NOV 20 11:59:59.001    38.189118
%      2008 MAY 18 11:59:58.999    37.223573
%      2008 NOV 14 11:59:59.001    36.084745
%      2009 MAY 13 11:59:57.999    36.608971
%
%
%-Particulars
%
%   The Mice routine cspice_conics is the inverse of this routine:
%   cspice_conics maps a set of osculating elements and a time to a state
%   vector.
%
%-Exceptions
%
%   1)  If `mu' is not positive, the error SPICE(NONPOSITIVEMASS)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If the specific angular momentum vector derived from `state'
%       is the zero vector, the error SPICE(DEGENERATECASE)
%       is signaled by a routine in the call tree of this routine.
%
%   3)  If the position or velocity vectors derived from `state'
%       is the zero vector, the error SPICE(DEGENERATECASE)
%       is signaled by a routine in the call tree of this routine.
%
%   4)  If the inclination is determined to be zero or 180 degrees,
%       the longitude of the ascending node is set to zero.
%
%   5)  If the eccentricity is determined to be zero, the argument of
%       periapse is set to zero.
%
%   6)  If the eccentricity of the orbit is very close to but not
%       equal to zero, the argument of periapse may not be accurately
%       determined.
%
%   7)  For inclinations near but not equal to 0 or 180 degrees,
%       the longitude of the ascending node may not be determined
%       accurately. The argument of periapse and mean anomaly may
%       also be inaccurate.
%
%   8)  For eccentricities very close to but not equal to 1, the
%       results of this routine are unreliable.
%
%   9)  If the specific angular momentum vector is non-zero but
%       "close" to zero, the results of this routine are unreliable.
%
%   10) If `state' is expressed relative to a non-inertial reference
%       frame, the resulting elements are invalid. No error checking
%       is done to detect this problem.
%
%   11) If any of the input arguments, `state', `et' or `mu', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   12) If any of the input arguments, `state', `et' or `mu', is not
%       of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   13) If the input vectorizable arguments `state' and `et' do not
%       have the same measure of vectorization (N), an error is
%       signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The input state vector must be expressed relative to an
%       inertial reference frame.
%
%   2)  Osculating elements are generally not useful for
%       high-accuracy work.
%
%   3)  Accurate osculating elements may be difficult to derive for
%       near-circular or near-equatorial orbits. Osculating elements
%       for such orbits should be used with caution.
%
%   4)  Extracting osculating elements from a state vector is a
%       mathematically simple but numerically challenging task. The
%       mapping from a state vector to equivalent elements is
%       undefined for certain state vectors, and the mapping is
%       difficult to implement with finite precision arithmetic for
%       states near the subsets of R6 where singularities occur.
%
%       In general, the elements found by this routine can have
%       two kinds of problems:
%
%       -  The elements are not accurate but still represent
%          the input state accurately. The can happen in
%          cases where the inclination is near zero or 180
%          degrees, or for near-circular orbits.
%
%       -  The elements are garbage. This can occur when
%          the eccentricity of the orbit is close to but
%          not equal to 1. In general, any inputs that cause
%          great loss of precision in the computation of the
%          specific angular momentum vector or the eccentricity
%          vector will result in invalid outputs.
%
%       For further details, see the -Exceptions section.
%
%       Users of this routine should carefully consider whether
%       it is suitable for their applications. One recommended
%       "sanity check" on the outputs is to supply them to the
%       Mice routine cspice_conics and compare the resulting state
%       vector with the one supplied to this routine.
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
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       examples' problem statement and meta-kernel. Updated code
%       examples to produce formatted output. Reduced the interval time
%       and steps to compute the solution, and fixed a bug in example #2.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 23-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   conic elements from state
%   osculating elements from state
%   convert state to osculating elements
%
%-&

function [elts] = cspice_oscelt( state, et, mu )

   switch nargin
      case 3

         state = zzmice_dp(state);
         et    = zzmice_dp(et);
         mu    = zzmice_dp(mu);

      otherwise

         error( 'Usage: [_elts(8)_] = cspice_oscelt( _state(6)_, _et_, mu )' )

   end

   %
   % Call the MEX library.
   %
   try
      [elts] = mice('oscelt_c', state, et, mu );
   catch spiceerr
      rethrow(spiceerr)
   end


