%-Abstract
%
%   CSPICE_VDIST returns the distance between two
%   three-dimensional vectors.
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
%      v1       an arbitrary vector(s).
%
%               [3,n] = size(v1); double = class(v1)
%
%      v2       also an arbitrary vector(s).
%
%               [3,n] = size(v2); double = class(v2)
%
%   the call:
%
%      [vdist] = cspice_vdist(v1, v2)
%
%   returns:
%
%      vdist    the value(s) describing the distance(s) between `v1' and `v2',
%               distance defined as:
%
%                   ||  v1 - v2  ||,
%
%                      _                                               _
%               where || x || indicates the Euclidean norm of the vector x.
%
%               [1,n] = size(vdist); double = class(vdist)
%
%               `vdist' returns with the same vectorization measure, N,
%               as `v1' and `v2'.
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
%   1) Define a set of vectors and calculate the distance between them.
%
%      Example code begins here.
%
%
%      function vdist_ex1()
%         %
%         % Define a set of vectors, calculate the distance
%         % between the coordinates.
%         %
%         v1 = [1; 0; 0];
%         v2 = [0; 1; 0];
%
%         vdist = cspice_vdist( v1, v2 );
%         disp( 'Scalar:' )
%         fprintf( '  %12.6f\n', vdist )
%
%         %
%         % Instead of two calls with 3-vectors,
%         % vectorize the input as two 3X2 array.
%         %
%         v1 = [ [1; 0; 0], [1; 0; 0] ];
%         v2 = [ [1; 0; 0], [0; 1; 0] ];
%
%         vdist = cspice_vdist( v1, v2 );
%         disp( 'Vectorized:' )
%         for i=1:2
%            fprintf( '  %12.6f\n', vdist(i))
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Scalar:
%            1.414214
%      Vectorized:
%            0.000000
%            1.414214
%
%
%   2) Given the planetocentric coordinates of a point on the surface of
%      Mars, compute the distance between that point and Phobos.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: vdist_ex2.tm
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
%            de430.bsp                        Planetary ephemeris
%            mar097.bsp                       Mars satellite ephemeris
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0011.tls                     Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de430.bsp',
%                                'mar097.bsp',
%                                'pck00010.tpc',
%                                'naif0011.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function vdist_ex2()
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'vdist_ex2.tm' );
%
%         %
%         % Define the point on the surface of Mars by its planetocentric
%         % coordinates, and the epoch.
%         %
%         epoch  = '2018-07-25 17:14';
%         lon    =    8.544377 * cspice_rpd;
%         lat    =   42.880602 * cspice_rpd;
%         radius = 3380.0;
%
%         %
%         % Convert that point coordinates to rectangular.
%         %
%         [rover] = cspice_latrec( radius, lon, lat );
%
%         %
%         % Convert the UTC epoch to ephemeris time.
%         %
%         [et] = cspice_str2et( epoch );
%
%         %
%         % Compute the position of Phobos with respect to Mars in IAU_MARS
%         % body-fixed reference frame.
%         %
%         [pos, lt] = cspice_spkpos( 'PHOBOS', et,   'IAU_MARS', ...
%                                      'NONE',  'MARS'           );
%
%         %
%         % Compute the distance between Phobos and the point on the surface
%         % of Mars.
%         %
%         dist = cspice_vdist( rover, pos );
%         fprintf( ' Epoch:  %s\n', epoch )
%         fprintf( ' Distance between location and Phobos (km):  %11.5f\n', ...
%                                                                    dist )
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
%       Epoch:  2018-07-25 17:14
%       Distance between location and Phobos (km):   7174.78139
%
%
%-Particulars
%
%   This function is simply shorthand for the code
%
%      diff = v1 - v2;
%
%      dist = cspice_vnorm( diff );
%
%   Using this function saves you the annoyance of declaring local
%   storage for the difference vector `diff'.
%
%
%   The Euclidean norm of a three-dimensional vector (x, y, z) is
%   defined as
%
%                                   1/2
%           2        2        2
%      (   x    +   y    +   z    ).
%
%
%   This number is the distance of the point (x, y, z) from the
%   origin. If `a' and `b' are two vectors whose components are
%
%      ( a(1), a(2), a(3) )    and    ( b(1), b(2), b(3) ),
%
%   then the distance between `a' and `b' is the norm of the difference
%   a - b, which has components
%
%
%      (  a(1) - b(1),  a(2) - b(2),  a(3) - b(3)  ).
%
%-Exceptions
%
%   1)  If any of the input arguments, `v1' or `v2', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `v1' or `v2', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   3)  If the input vectorizable arguments `v1' and `v2' do not have
%       the same measure of vectorization (N), an error is signaled by
%       the Mice interface.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 25-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example #1's problems statement, and a second complete example.
%
%       Changed output argument name "dist" to "vdist" to comply with
%       NAIF standard.
%
%       Added -Parameters, -Particulars -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 18-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   distance between 3-dimensional vectors
%
%-&

function [vdist] = cspice_vdist( v1, v2 )

   switch nargin
      case 2

         v1 = zzmice_dp(v1);
         v2 = zzmice_dp(v2);

      otherwise

         error ( 'Usage: [_vdist_] = cspice_vdist( _v1(3)_, _v2(3)_ )' )

   end

   %
   % Call the MEX library.
   %
   try
      [vdist] = mice('vdist_c',v1, v2);
   catch spiceerr
      rethrow(spiceerr)
   end
