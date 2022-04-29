%-Abstract
%
%   CSPICE_HALFPI returns the double precision value of the constant pi/2.0.
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
%   The call:
%
%      [halfpi] = cspice_halfpi
%
%   returns:
%
%      halfpi   half the value of pi (the ratio of a circle's circumference
%               to its diameter), determined by the acos function.
%
%               [1,1] = size(halfpi); double = class(halfpi)
%
%               That is,
%
%                     halfpi = acos ( -1.0 ) * 0.5;
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
%   1) The following code example returns the double precision value of
%      the constant pi/2.0 and prints it out.
%
%      Example code begins here.
%
%
%      function halfpi_ex1()
%
%         %
%         % Print the double precision value of pi/2.0
%         %
%         fprintf( 'Half pi: %25.22f\n', cspice_halfpi )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Half pi:  1.5707963267948965579990
%
%
%   2) Compute the transformation from inertial to body fixed
%      coordinates, given the directions of the north pole and prime
%      meridian of the body.
%
%      When using the following values for Pluto, extracted from the
%      PCK kernel pck00010.tpc:
%
%         Right ascension (deg): 132.993
%         Declination     (deg):  -6.163
%         Prime meridian  (deg): 302.695
%
%      at ephemeris epoch 2000 Jan 1 12:00:00 TDB, the result should
%      match that obtained using the following call:
%
%         [tipm] = cspice_pxform ( "J2000", "IAU_PLUTO", 0.0 );
%
%      Use the PCK kernel below to load the triaxial ellipsoidal shape
%      model and orientation data for Pluto.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function halfpi_ex2()
%
%         %
%         % Load the PCK.
%         %
%         cspice_furnsh( 'pck00010.tpc' );
%
%         %
%         % Compute the transformation from inertial to body
%         % fixed coordinates, given the directions of the north
%         % pole and prime meridian of the body.
%         %
%
%         %
%         % Assign the values for Pluto, in radians.
%         %
%         ra  = 132.993 * cspice_rpd;
%         dec =  -6.163 * cspice_rpd;
%         w   = 302.695 * cspice_rpd;
%
%         %
%         % The transformation is defined by the compound
%         % rotation
%         %
%         %   [W] [pi/2 - Dec] [RA + pi/2]
%         %      3            1           3
%         %
%         [tipm] = cspice_rotate( ra + cspice_halfpi, 3 );
%         [tipm] = cspice_rotmat( tipm, cspice_halfpi - dec, 1 );
%         [tipm] = cspice_rotmat( tipm, w, 3 );
%
%         %
%         % Print the results
%         %
%         fprintf( 'Rotation matrix, from pole direction and prime\n' )
%         fprintf( 'meridian:\n' )
%         fprintf( '    %11.5f  %11.5f  %11.5f\n', ...
%                                     tipm(1,1), tipm(1,2), tipm(1,3) )
%         fprintf( '    %11.5f  %11.5f  %11.5f\n', ...
%                                     tipm(2,1), tipm(2,2), tipm(2,3) )
%         fprintf( '    %11.5f  %11.5f  %11.5f\n', ...
%                                     tipm(3,1), tipm(3,2), tipm(3,3) )
%
%         %
%         % Use pxform_c to obtain the same transformation.
%         %
%         [tipm] = cspice_pxform( 'J2000', 'IAU_PLUTO', 0.0 );
%         fprintf( '\n' )
%         fprintf( 'Rotation matrix, from pxform_c:\n' )
%         fprintf( '    %11.5f  %11.5f  %11.5f\n', ...
%                                     tipm(1,1), tipm(1,2), tipm(1,3) )
%         fprintf( '    %11.5f  %11.5f  %11.5f\n', ...
%                                     tipm(2,1), tipm(2,2), tipm(2,3) )
%         fprintf( '    %11.5f  %11.5f  %11.5f\n', ...
%                                     tipm(3,1), tipm(3,2), tipm(3,3) )
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
%      Rotation matrix, from pole direction and prime
%      meridian:
%             -0.33349     -0.43443     -0.83669
%             -0.65509     -0.53145      0.53704
%             -0.67797      0.72721     -0.10736
%
%      Rotation matrix, from pxform_c:
%             -0.33349     -0.43443     -0.83669
%             -0.65509     -0.53145      0.53704
%             -0.67797      0.72721     -0.10736
%
%
%-Particulars
%
%   The first time the function is referenced, the value is computed
%   as shown above. The value is saved, and returned directly upon
%   subsequent reference.
%
%-Exceptions
%
%   Error free.
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
%       Edited the header to comply with NAIF standard: created
%       a complete example using the existing code fragments and added
%       example's problem statement. Added example #2.
%
%       Changed output argument name "return_val" to "halfpi" to comply with
%       NAIF standard.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   half the value of pi
%
%-&

function [halfpi] = cspice_halfpi

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: [halfpi] = cspice_halfpi' )

   end

   %
   % Call the MEX library.
   %
   try
      [halfpi] =  mice('halfpi_c');
   catch spiceerr
      rethrow(spiceerr)
   end
