%-Abstract
%
%   CSPICE_XF2RAV determines the rotation matrix and angular velocity of the
%   rotation from a state transformation matrix.
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
%      xform    a state transformation matri(x|ces) from one frame
%               FRAME1 to some other frame FRAME2.
%
%               Either [6,6]   = size(xform); double = class(xform)
%               or     [6,6,n] = size(xform); double = class(xform)
%
%   the call:
%
%      [rot, av] = cspice_xf2rav( xform )
%
%   returns:
%
%      rot      rotation matri(x|ces) that gives the transformation from
%               some frame FRAME1 to another frame FRAME2.
%
%               Either [3,3]   = size(rot); double = class(rot)
%               or     [3,3,n] = size(rot); double = class(rot)
%
%      av       the angular velocity vector(s) of the transformation.
%
%               Either [3,1] = size(av); double = class(av)
%               or     [3,n] = size(av); double = class(av)
%
%               In other words, if `p' is the position of a fixed point in
%               FRAME2, then from the point of view of FRAME1, `p' rotates
%               (in a right handed sense) about an axis parallel to `av'.
%               Moreover the rate of rotation in radians per unit time is
%               given by the length of `av'.
%
%               More formally, the velocity `v' of `p' in FRAME1 is
%               given by
%                                  T
%                   v  = av x ( rot  * p )
%
%               The components of `av' are given relative to FRAME1.
%
%               `rot' and `av' return with the same vectorization
%               measure, N, as `xform'.
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
%   1) Suppose that you wanted to determine the angular velocity
%      of the Earth body-fixed reference frame with respect to
%      J2000 at a particular epoch ET. The following code example
%      illustrates a procedure for computing the angular velocity.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: xf2rav_ex1.tm
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
%            earth_720101_070426.bpc       Earth historical
%                                          binary PCK
%            naif0012.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'earth_720101_070426.bpc',
%                                'naif0012.tls'            )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function xf2rav_ex1()
%
%         %
%         % Local parameters.
%         %
%         META   =   'xf2rav_ex1.tm';
%         UTCSTR =   '2005-OCT-10 16:00:00';
%
%         %
%         % Load SPICE kernels.
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert the input time to seconds past J2000 TDB.
%         %
%         [et] = cspice_str2et( UTCSTR );
%
%         %
%         % Get the transformation matrix from J2000 frame to
%         % ITRF93.
%         %
%         [ftmtrx] = cspice_sxform( 'J2000', 'ITRF93', et );
%
%         %
%         % Now get the angular velocity by calling cspice_xf2rav
%         %
%         [rot, av] = cspice_xf2rav( ftmtrx );
%
%         %
%         % Display the results.
%         %
%         fprintf( 'Rotation matrix:\n' )
%         fprintf( '%16.11f %15.11f %15.11f\n', rot' )
%
%         fprintf( '\n' )
%         fprintf( 'Angular velocity:\n' )
%         fprintf( '%16.11f %15.11f %15.11f\n', av )
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
%      Rotation matrix:
%        -0.18603277688  -0.98254352801   0.00014659080
%         0.98254338275  -0.18603282936  -0.00053610915
%         0.00055402128   0.00004429795   0.99999984555
%
%      Angular velocity:
%         0.00000004025   0.00000000324   0.00007292114
%
%
%   2) Obtain the state transformation matrix from J2000 to IAU_MOON for
%      a set of 10001 ephemeris times based at July 1 2007, convert them
%      to the corresponding rotation matrices and angular velocity
%      vectors and back to state transformation matrices.
%
%      Compare the original state transformation matrices with those
%      computed, and output the maximum absolute difference between any
%      of them.
%
%      Numerical equivalence shall be expected.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: xf2rav_ex2.tm
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
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0012.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'pck00010.tpc',
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
%      function xf2rav_ex2()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'xf2rav_ex2.tm')
%
%         %
%         % Create an array of 10001 ephemeris times based at July 1 2007.
%         %
%         et    = [0: 10000]* cspice_spd + cspice_str2et( 'July 1 2007' );
%
%         %
%         % Calculate the state transformation matrices from J2000 to IAU_MOON
%         % for 'et'.
%         %
%         xform = cspice_sxform( 'J2000', 'IAU_MOON', et );
%
%         %
%         % Convert the set of 'xform' matrices to the corresponding rotation
%         % matrices and angular velocity vectors.
%         %
%         [ rot, av ] = cspice_xf2rav(xform);
%
%         %
%         % Use the converted outputs from cspice_xf2rav to recompute a set
%         % of state transformation matrices.
%         %
%         strans = cspice_rav2xf( rot, av );
%
%         %
%         % Calculate the maximum value of the absolute difference between
%         % 'xform' and 'strans'.
%         %
%         fprintf('Maximum absolute difference: %8.6e\n', ...
%                  max( max( max( abs(strans - xform) ) ) )   )
%
%         %
%         %  It's always good form to unload kernels after use,
%         %  particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Maximum absolute difference: 2.117582e-21
%
%
%-Particulars
%
%   This routine is essentially a macro routine for converting
%   state transformation matrices into the equivalent representation
%   in terms of a rotation and angular velocity.
%
%   This routine is an inverse of the routine cspice_rav2xf.
%
%-Exceptions
%
%   1)  No checks are performed on `xform' to ensure that it is indeed
%       a state transformation matrix.
%
%   2)  If the input argument `xform' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `xform' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
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
%   ROTATION.REQ
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added first example,
%       and second example's problem statement and meta-kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       extended -Particulars section. Added rotation.req to Required
%       readings.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 09-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 11-APR-2007 (EDW)
%
%-Index_Entries
%
%   State transformation to rotation and angular velocity
%
%-&

function [rot, av] = cspice_xf2rav(xform)

   switch nargin
      case 1

         xform = zzmice_dp(xform);

      otherwise

         error ( ['Usage: [_rot(3,3)_, _av(3)_] = ' ...
                  'cspice_xf2rav(_xform(6,6)_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [rot, av] = mice('xf2rav_c', xform);
   catch spiceerr
      rethrow(spiceerr)
   end


