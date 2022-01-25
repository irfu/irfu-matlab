%-Abstract
%
%   CSPICE_DNEARP computes the state (position and velocity) of an ellipsoid
%   surface point nearest to the position component of a specified state.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      state    a 6-vector giving the position and velocity of some object in
%               the body-fixed coordinates of the ellipsoid.
%
%               [6,1] = size(state); double = class(state)
%
%               In body-fixed coordinates, the semi-axes of the ellipsoid
%               are aligned with the X, Y, and Z-axes of the coordinate
%               system.
%
%      a        the length of the semi-axis of the ellipsoid that is parallel
%               to the X-axis of the body-fixed coordinate system.
%
%               [1,1] = size(a); double = class(a)
%
%      b        the length of the semi-axis of the ellipsoid that is parallel
%               to the Y-axis of the body-fixed coordinate system.
%
%               [1,1] = size(b); double = class(b)
%
%      c        the length of the semi-axis of the ellipsoid that is parallel
%               to the Z-axis of the body-fixed coordinate system.
%
%               [1,1] = size(c); double = class(c)
%
%   the call:
%
%      [dnear, dalt, found] = cspice_dnearp( state, a, b, c )
%
%   returns:
%
%      dnear    the 6-vector giving the position and velocity in body-fixed
%               coordinates of the point on the ellipsoid, closest to the
%               object whose position and velocity are represented by
%               `state'.
%
%               [6,1] = size(dnear); double = class(dnear)
%
%               While the position component of `dnear' is always
%               meaningful, the velocity component of `dnear' will be
%               meaningless if `found' if false (See the discussion of
%               the meaning of `found' below.)
%
%      dalt     an array of two double precision numbers.
%
%               [2,1] = size(dalt); double = class(dalt)
%
%               The first gives the altitude of `state' with respect to the
%               ellipsoid. The second gives the rate of change of the
%               altitude.
%
%               Note that the rate of change of altitude is meaningful if
%               and only if `found' is true (See the discussion of the
%               meaning of `found' below.)
%
%      found    a logical flag indicating whether or not the velocity portion
%               of `dnear' is meaningful.
%
%               [1,1] = size(found); logical = class(found)
%
%               If the velocity portion of `dnear' is meaningful `found'
%               will be returned with a value of true. Under very rare
%               circumstance the velocity of the near point is undefined.
%               Under these circumstances `found' will be returned with the
%               value false.
%
%               `found' can be false only for states whose position
%               components are inside the ellipsoid and then only at
%               points on a special surface contained inside the
%               ellipsoid called the focal set of the ellipsoid.
%
%               A point in the interior is on this special surface only
%               if there are two or more points on the ellipsoid that are
%               closest to it. The origin is such a point and the only
%               such point if the ellipsoid is a sphere. For
%               non-spheroidal ellipsoids the focal set contains small
%               portions of the planes of symmetry of the ellipsoid.
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
%   1) Suppose you wish to compute the velocity of the ground track
%      of a satellite as it passes over a location on Mars and that
%      the moment of passage has been previously determined. (We
%      assume that the spacecraft is close enough to the surface that
%      light time corrections do not matter.)
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: dnearp_ex1.tm
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
%            File name                        Contents
%            ---------                        --------
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0012.tls                     Leapseconds
%            de430.bsp                        Planetary ephemeris
%            mar097.bsp                       Mars satellite ephemeris
%            mro_psp4_ssd_mro95a.bsp          MRO ephemeris
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'pck00010.tpc',
%                                'naif0012.tls',
%                                'de430.bsp',
%                                'mar097.bsp',
%                                'mro_psp4_ssd_mro95a.bsp' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function dnearp_ex1()
%
%         %
%         % Local parameters
%         %
%         BODYNM =   'MARS';
%         META   =   'dnearp_ex1.tm';
%
%         %
%         % Load kernel files via the meta-kernel.
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert the TDB input time string to seconds past
%         % J2000, TDB.
%         %
%         [et] = cspice_str2et( '2007 SEP 30 00:00:00 TDB' );
%
%         %
%         % First get the axes of the body.
%         %
%         [radii] = cspice_bodvrd( BODYNM, 'RADII', 3 );
%         [a, b, c] = cspice_vupack( radii );
%
%         %
%         % Get the geometric state of the spacecraft with
%         % respect to BODYNM in the body-fixed reference frame
%         % at `et' and compute the state of the sub-spacecraft point.
%         %
%         [state, lt] = cspice_spkezr( 'MRO',  et,     'IAU_MARS',         ...
%                                      'NONE', BODYNM              );
%         [dnear, dalt, found] = cspice_dnearp( state, a, b, c );
%
%         if ( found )
%
%            %
%            % `dnear' contains the state of the subspacecraft point.
%            %
%            gtvel = dnear(4:6);
%
%            fprintf( 'Ground-track velocity (km/s): %9.6f %9.6f %9.6f\n', ...
%                                                                  gtvel' )
%            fprintf( 'Ground-track speed    (km/s): %9.6f\n',             ...
%                                        cspice_vnorm( gtvel ) )
%
%         else
%
%            fprintf( 'DNEAR is degenerate.\n' )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Ground-track velocity (km/s):  0.505252  1.986553 -2.475506
%      Ground-track speed    (km/s):  3.214001
%
%
%   2) Suppose you wish to compute the one-way doppler shift of a
%      radar mounted on board a spacecraft as it passes over some
%      region. Moreover, assume that for your purposes it is
%      sufficient to neglect effects of atmosphere, topography and
%      antenna pattern for the sake of this computation.
%
%      Use the meta-kernel from Example 1 above.
%
%
%      Example code begins here.
%
%
%      function dnearp_ex2()
%
%         %
%         % Local parameters
%         %
%         BODYNM =   'MARS';
%         META   =   'dnearp_ex1.tm';
%
%         %
%         % Define the central frequency of the radar,
%         % in megahertz.
%         %
%         RCFRQ =   20.0;
%
%         %
%         % Load kernel files via the meta-kernel.
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert the TDB input time string to seconds past
%         % J2000, TDB.
%         %
%         [et] = cspice_str2et( '2007 SEP 30 00:00:00 TDB' );
%
%         %
%         % First get the axes of the body.
%         %
%         [radii] = cspice_bodvrd( BODYNM, 'RADII', 3 );
%         [a, b, c] = cspice_vupack( radii );
%
%         %
%         % Get the geometric state of the spacecraft with
%         % respect to BODYNM in the body-fixed reference frame
%         % at `et' and compute the state of the sub-spacecraft point.
%         %
%         [state, lt] = cspice_spkezr( 'MRO',  et,     'IAU_MARS',         ...
%                                      'NONE', BODYNM              );
%         [dnear, dalt, found] = cspice_dnearp( state, a, b, c );
%
%         if ( found )
%
%            %
%            % The change in frequency is given by multiplying `shift'
%            % times the carrier frequency
%            %
%            shift = ( dalt(2) / cspice_clight );
%            fprintf( 'Central frequency (MHz): %19.16f\n', RCFRQ )
%            fprintf( 'Doppler shift     (MHz): %19.16f\n', RCFRQ * shift )
%
%         else
%
%            fprintf( 'DNEAR is degenerate.\n' )
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Central frequency (MHz): 20.0000000000000000
%      Doppler shift     (MHz): -0.0000005500991159
%
%
%-Particulars
%
%   If an object is moving relative to some triaxial body along a
%   trajectory c(t) then there is a companion trajectory n(t) that
%   gives the point on the ellipsoid that is closest to c(t) as a
%   function of `t'. The instantaneous position and velocity of c(t),
%   `state', are sufficient to compute the instantaneous position and
%   velocity of n(t), `dnear'.
%
%   This routine computes `dnear' from `state'. In addition it returns the
%   altitude and rate of change of altitude.
%
%   Note that this routine can compute `dnear' for `state' outside, on,
%   or inside the ellipsoid. However, the velocity of `dnear' and
%   derivative of altitude do not exist for a "small" set of `state'
%   in the interior of the ellipsoid. See the discussion of `found'
%   above for a description of this set of points.
%
%-Exceptions
%
%   1)  If the axes are non-positive, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If an object is passing through the interior of an ellipsoid
%       there are points at which there is more than 1 point on the
%       ellipsoid that is closest to the object. At these points the
%       velocity of the near point is undefined. (See the description
%       of the output variable `found').
%
%   3)  If any of the input arguments, `state', `a', `b' or `c', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `state', `a', `b' or `c', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
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
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 01-NOV-2021 (JDR)
%
%-Index_Entries
%
%   Velocity of the nearest point on an ellipsoid
%   Rate of change of the altitude over an ellipsoid
%   Derivative of altitude over an ellipsoid
%   Velocity of a ground track
%
%-&
function [dnear, dalt, found] = cspice_dnearp( state, a, b, c )

   switch nargin
      case 4

         state = zzmice_dp(state);
         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);

      otherwise

         error ( [ 'Usage: [dnear(6), dalt(2), found] = '                   ...
                   'cspice_dnearp( state(6), a, b, c )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [dnear, dalt, found] = mice('dnearp_c', state, a, b, c);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
