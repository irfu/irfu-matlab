%-Abstract
%
%   CSPICE_M2EUL factors a rotation matrix into a product of
%   three rotations about specified coordinate axes.
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
%      r        the rotation matrix/matrices to factor as a product of three
%               rotations about specified coordinate axes.
%
%               [3,3]   = size(r); double = class(r)
%
%                  or
%
%               [3,3,n] = size(r); double = class(r)
%
%               The angles of these rotations are called "Euler angles."
%
%      axis3,
%      axis2,
%      axis1    the indices defining the rotation axes of the
%               "factor" rotations, whose product is `r'.
%
%               [1,1] = size(axis3); int32 = class(axis3)
%               [1,1] = size(axis2); int32 = class(axis2)
%               [1,1] = size(axis1); int32 = class(axis1)
%
%               `r' is factored as
%
%                  r = [ angle3 ]      [ angle2 ]      [ angle1 ]
%                                axis3           axis2           axis1
%
%               The axis numbers must belong to the set {1, 2, 3}.
%
%               The values of axisX may be 1, 2, or 3, indicating
%               the x, y, and z axes respectively. The second axis number
%               MUST differ from the first and third axis numbers.
%
%   the call:
%
%      [angle3, angle2, angle1] = cspice_m2eul( r, axis3, axis2, axis1 )
%
%   returns:
%
%      angle3,
%      angle2,
%      angle1   the Euler angles measured where the angles satisfy
%
%                  r = [ angle3 ]     [ angle2 ]     [ angle1 ]
%                               axis3          axis2          axis1
%
%               If [3,3] = size(r)
%               then
%
%               [1,1] = size(angle3); double = class(angle3)
%               [1,1] = size(angle2); double = class(angle2)
%               [1,1] = size(angle1); double = class(angle1)
%
%               If [3,3,n] = size(r)
%               then
%
%               [1,n] = size(angle3); double = class(angle3)
%               [1,n] = size(angle2); double = class(angle2)
%               [1,n] = size(angle1); double = class(angle1)
%
%               The range of `angle3' and `angle1' is (-pi, pi].
%
%               The range of `angle2' depends on the exact set of
%               axes used for the factorization. For
%               factorizations in which the first and third axes
%               are the same,
%
%                  r = [ angle3 ]  [ angle2 ]  [ angle1 ]
%                               a           b           a
%
%               the range of `angle2' is [0, pi].
%
%               For factorizations in which the first and third
%               axes are different,
%
%                  r = [ angle3 ]  [ angle2 ]  [ angle1 ],
%                               a           b           c
%
%               the range of angle2 is [-pi/2, pi/2].
%
%               For rotations such that `angle3' and `angle1' are not
%               uniquely determined, `angle3' will always be set to
%               zero; `angle1' is then uniquely determined.
%
%               `angle3', `angle2', and `angle1' return with the same
%               vectorization measure, N, as `r'.
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
%   1) Conversion of instrument pointing from a matrix representation
%      to Euler angles:
%
%      Suppose we want to find camera pointing in alpha, delta, and
%      kappa, given the inertial-to-camera coordinate transformation
%
%
%      .-                                                         -.
%      |  0.491273796781358  0.508726203218642  0.706999085398824  |
%      |                                                           |
%      | -0.508726203218642 -0.491273796781358  0.706999085398824  |
%      |                                                           |
%      |  0.706999085398824 -0.706999085398824  0.017452406437284  |
%      `-                                                         -'
%
%
%      Find angles alpha, delta, kappa such that
%
%         ticam  =  [ kappa ]  [ pi/2 - delta ]  [ pi/2 + alpha ] .
%                            3                 1                 3
%
%      Example code begins here.
%
%
%      function m2eul_ex1()
%
%         %
%         % Scalar example, conversion of instrument pointing from a matrix
%         % representation to Euler angles:
%         %
%         % Suppose we want to find camera pointing in `alpha', `delta', and
%         % `kappa', given the inertial-to-camera coordinate transformation
%         %
%         ticam = [                                                        ...
%             [  0.491273796781358  0.508726203218642  0.706999085398824 ]
%             [ -0.508726203218642 -0.491273796781358  0.706999085398824 ]
%             [  0.706999085398824 -0.706999085398824  0.017452406437284 ] ];
%
%         %
%         % We want to find angles alpha, delta, kappa such that
%         %
%         %
%         %   ticam  =  [ kappa ]  [ pi/2 - delta ]  [ pi/2 + alpha ] .
%         %                      3                 1                 3
%         %
%         %
%         % Factor the matrix to the Euler angles corresponding to a
%         % 3-1-3 rotation.
%         %
%         [ kappa, ang2, ang1  ] = cspice_m2eul( ticam, 3, 1, 3 );
%
%         alpha  =  ang1          - cspice_halfpi;
%         delta  =  cspice_halfpi - ang2;
%
%         %
%         % Calculate the desired angles. If we wish to make sure that
%         % alpha, delta, and kappa are in the ranges [0, 2pi),
%         % [-pi/2, pi/2], and [0, 2pi) respectively, we may add the code
%         %
%
%         if ( alpha < 0. )
%            alpha = alpha + cspice_twopi;
%         end
%
%         if ( kappa < 0. )
%            kappa = kappa + cspice_twopi;
%         end
%
%         %
%         % Output the 3-1-3 Euler rotation angles corresponding to `ticam'.
%         %
%         fprintf( 'Alpha (deg):  %23.14f\n', cspice_dpr * alpha )
%         fprintf( 'Delta (deg):  %23.14f\n', cspice_dpr * delta )
%         fprintf( 'Kappa (deg):  %23.14f\n', cspice_dpr * kappa )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Alpha (deg):       315.00000000000000
%      Delta (deg):         1.00000000000003
%      Kappa (deg):        45.00000000000000
%
%
%   2) Conversion of instrument pointing angles from a non-J2000,
%      not necessarily inertial frame to J2000-relative RA, Dec,
%      and Twist.
%
%      Suppose that we have orientation for the CASSINI Narrow Angle
%      Camera (NAC) frame expressed as
%
%         [ gamma ]  [ beta ]  [ alpha ]
%                  1         2          3
%
%      with respect to the CASSINI spacecraft frame.
%
%      We want to express that orientation with respect to the J2000
%      frame as the sequence of rotations
%
%         [ twist ]  [ dec ]  [ ra ] .
%                  3        1       3
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: m2eul_ex2.tm
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
%            File name                      Contents
%            ---------                      --------
%            naif0010.tls                   Leapseconds
%            cas00145.tsc                   Cassini SCLK
%            cas_v40.tf                     Cassini frames
%            08052_08057ra.bc               Orientation for Cassini
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0010.tls'
%                                'cas00145.tsc'
%                                'cas_v40.tf'
%                                '08052_08057ra.bc')
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function m2eul_ex2()
%
%         %
%         % Local parameters
%         %
%         META  = 'm2eul_ex2.tm';
%         ALPHA =  89.9148;
%         BETA  =  -0.03300;
%         GAMMA = -90.009796;
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( META );
%
%         %
%         % First, we use function cspice_eul2m to form the
%         % transformation from instrument coordinates to
%         % CASSINI spacecraft frame.
%         %
%         [inst2s] = cspice_eul2m( GAMMA * cspice_rpd, BETA * cspice_rpd,  ...
%                                  ALPHA * cspice_rpd, 1,  2,    3       );
%
%         %
%         % Now we compute the transformation from CASSINI
%         % spacecraft frame to J2000, at a given time.
%         %
%         [et]  = cspice_str2et( '2008 Feb 25' );
%         [s2j] = cspice_pxform( 'CASSINI_SC_COORD', 'J2000', et );
%
%         %
%         % Next, we form the transformation from instrument
%         % coordinates to J2000 frame.
%         %
%         inst2j = s2j * inst2s;
%
%         %
%         % Finally, we express `inst2j' using the desired Euler
%         % angles.
%         %
%         [twist, ang1, ang2] = cspice_m2eul( inst2j, 3, 1, 3 );
%
%         ra  =  ang2 - cspice_halfpi;
%         dec =  cspice_halfpi - ang1;
%
%         %
%         % If we wish to make sure that `ra', `dec', and `twist' are in
%         % the ranges [0, 2pi), [-pi/2, pi/2], and [0, 2pi)
%         % respectively, we may add the code
%         %
%         if ( ra < 0.0 )
%            ra    = ra + cspice_twopi;
%         end
%
%         if ( twist < 0.0 )
%            twist = twist + cspice_twopi;
%         end
%
%         %
%         % Now `ra', `dec', and `twist' express the instrument pointing
%         % as RA, Dec, and Twist, relative to the J2000 system.
%         %
%         fprintf( 'RA    (deg):  %23.14f\n', cspice_dpr * ra )
%         fprintf( 'Dec   (deg):  %23.14f\n', cspice_dpr * dec )
%         fprintf( 'Twist (deg):  %23.14f\n', cspice_dpr * twist )
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
%      RA    (deg):        83.77802387778848
%      Dec   (deg):       -14.92925498590898
%      Twist (deg):       294.55732942050986
%
%
%      Note: more than one definition of RA, Dec, and Twist is
%      extant. Before using this example in an application, check
%      that the definition given here is consistent with that used
%      in your application.
%
%
%   3) Calculate the J2000 to IAU_MOON transformation matrices for
%      a set of ephemeris times starting from Jan 1, 2000 at noon,
%      and their corresponding 3-2-1 Euler rotation angles.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: m2eul_ex3.tm
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
%            naif0012.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
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
%      function m2eul_ex3()
%
%         %
%         % Vectorized example, input an array of ephemeris times, calculate
%         % the corresponding J2000 to IAU_MOON transformation matrices.
%         %
%         cspice_furnsh('m2eul_ex3.tm')
%
%         et0 = cspice_str2et( 'Jan 1 2000 12:00:00 TDB' );
%         et1 = cspice_str2et( 'Jan 1 2010 12:00:00 TDB' );
%
%         n     = 10;
%         times = et0 + (1:n)* (et1 - et0)/n;
%         quot   = cspice_pxform( 'J2000', 'IAU_MOON', times );
%
%         %
%         % Factor the matrices to the Euler angles corresponding to a
%         % 3-2-1 rotation set.
%         %
%         [a3,a2,a1] = cspice_m2eul( quot, 1,2,3);
%
%         %
%         % Output the 3-2-1 Euler rotation angles corresponding to `quot'.
%         %
%         fprintf( '%12.5f   %12.5f   %12.5f\n', [a1; a2; a3] * cspice_dpr )
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
%         171.75375       -2.00894      -23.83517
%         -52.93007       18.11962       15.07397
%          77.30266      -22.59555        3.51974
%        -150.68645       12.42680      -18.79120
%         -14.28248        4.91714       21.55874
%         120.06957      -19.09792      -11.00536
%        -109.73801       20.66329       -7.52692
%          23.54335       -8.43440       20.49467
%         160.13917       -9.11890      -20.58629
%         -66.71201       21.70068        7.52880
%
%
%-Particulars
%
%   A word about notation: the symbol
%
%      [ x ]
%           i
%
%   indicates a coordinate system rotation of `x' radians about the
%   ith coordinate axis. To be specific, the symbol
%
%      [ x ]
%           1
%
%   indicates a coordinate system rotation of `x' radians about the
%   first, or x-, axis; the corresponding matrix is
%
%      .-                    -.
%      |  1      0       0    |
%      |                      |
%      |  0    cos(x)  sin(x) |
%      |                      |
%      |  0   -sin(x)  cos(x) |
%      `-                    -'
%
%   Remember, this is a COORDINATE SYSTEM rotation by `x' radians; this
%   matrix, when applied to a vector, rotates the vector by -x
%   radians, not `x' radians. Applying the matrix to a vector yields
%   the vector's representation relative to the rotated coordinate
%   system.
%
%   The analogous rotation about the second, or y-, axis is
%   represented by
%
%      [ x ]
%           2
%
%   which symbolizes the matrix
%
%      .-                    -.
%      | cos(x)   0   -sin(x) |
%      |                      |
%      |  0       1      0    |
%      |                      |
%      | sin(x)   0    cos(x) |
%      `-                    -'
%
%   and the analogous rotation about the third, or z-, axis is
%   represented by
%
%      [ x ]
%           3
%
%   which symbolizes the matrix
%
%      .-                    -.
%      |  cos(x)  sin(x)   0  |
%      |                      |
%      | -sin(x)  cos(x)   0  |
%      |                      |
%      |  0        0       1  |
%      `-                    -'
%
%
%   The input matrix is assumed to be the product of three
%   rotation matrices, each one of the form
%
%      .-                    -.
%      |  1      0       0    |
%      |                      |
%      |  0    cos(r)  sin(r) |     (rotation of `r' radians about the
%      |                      |      x-axis),
%      |  0   -sin(r)  cos(r) |
%      `-                    -'
%
%
%      .-                    -.
%      | cos(s)   0   -sin(s) |
%      |                      |
%      |  0       1      0    |     (rotation of `s' radians about the
%      |                      |      y-axis),
%      | sin(s)   0    cos(s) |
%      `-                    -'
%
%   or
%
%      .-                    -.
%      |  cos(t)  sin(t)   0  |
%      |                      |
%      | -sin(t)  cos(t)   0  |     (rotation of `t' radians about the
%      |                      |      z-axis),
%      |  0        0       1  |
%      `-                    -'
%
%   where the second rotation axis is not equal to the first or
%   third. Any rotation matrix can be factored as a sequence of
%   three such rotations, provided that this last criterion is met.
%
%   This routine is related to the Mice routine cspice_eul2m, which
%   produces a rotation matrix, given a sequence of Euler angles.
%   This routine is a "right inverse" of cspice_eul2m: the sequence of
%   calls
%
%      [angle3, angle2, angle1] = cspice_m2eul( r, axis3, axis2, axis1 )
%      [r] = cspice_eul2m( angle3, angle2, angle1, axis3, axis2, axis1 )
%
%   preserves 'r' to round-off error.
%
%   On the other hand, the sequence of calls
%
%      [r] = cspice_eul2m( angle3, angle2, angle1, axis3, axis2, axis1 )
%      [angle3, angle2, angle1] = cspice_m2eul( r, axis3, axis2, axis1 )
%
%   preserves 'angle3', 'angle2', and 'angle1' only if the initial
%   values of the angle existed within the range of cspice_m2eul's
%   output.
%
%-Exceptions
%
%   1)  If any of `axis3', `axis2', or `axis1' do not have values in
%
%          { 1, 2, 3 }
%
%       the error SPICE(BADAXISNUMBERS) is signaled by a routine in
%       the call tree of this routine.
%
%   2)  If `axis2' is equal to `axis3' or `axis1', the error
%       SPICE(BADAXISNUMBERS) is signaled by a routine in the call
%       tree of this routine. An arbitrary rotation matrix cannot be
%       expressed using a sequence of Euler angles unless the second
%       rotation axis differs from the other two.
%
%   3)  If the input matrix `r' is not a rotation matrix, the error
%       SPICE(NOTAROTATION) is signaled by a routine in the call tree
%       of this routine.
%
%   4)  If `angle3' and `angle1' are not uniquely determined, `angle3'
%       is set to zero, and `angle1' is determined.
%
%   5)  If any of the input arguments, `r', `axis3', `axis2' or
%       `axis1', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   6)  If any of the input arguments, `r', `axis3', `axis2' or
%       `axis1', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
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
%       Edited the header to comply with NAIF standard. Added
%       examples' problem statement, a third example, and the meta-kernel
%       for code example #3.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 09-MAR-2015 (EDW)
%
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.1, 30-DEC-2008 (EDW)
%
%      Corrected misspellings.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   matrix to euler angles
%
%-&

function [angle3, angle2, angle1] = cspice_m2eul(r, axis3, axis2, axis1)

   switch nargin
      case 4

         r     = zzmice_dp(r);
         axis3 = zzmice_int(axis3);
         axis2 = zzmice_int(axis2);
         axis1 = zzmice_int(axis1);

      otherwise

         error ( ['Usage: [_angle3_, _angle2_, _angle1_]  = '  ...
                  'cspice_m2eul( _r(3,3)_, axis3, axis2, axis1)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [angle3, angle2, angle1] = mice('m2eul_c', r, axis3, axis2, axis1);
   catch spiceerr
      rethrow(spiceerr)
   end


