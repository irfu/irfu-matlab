%-Abstract
%
%   CSPICE_RAV2XF determines the state transformation matrix
%   from a rotation matrix and the angular velocity of the
%   rotation.
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
%      rot      a rotation matri(x|ces) that gives the transformation from
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
%                  v  =  av x ( rot  * p )
%
%   the call:
%
%      [xform] = cspice_rav2xf( rot, av )
%
%   returns:
%
%      xform    a state transformation matri(x|ces) associated with `rot'
%               and `av'.
%
%               Either [6,6]   = size(xform); double = class(xform)
%               or     [6,6,n] = size(xform); double = class(xform)
%
%               If `s1' is the state of an object with respect to FRAME1,
%               then the state `s2' of the object with respect to FRAME2 is
%               given by
%
%                  s2  =  xform * s1
%
%               where "*" denotes Matrix-Vector multiplication.
%
%               `xform' returns with the same vectorization measure, N,
%               as `rot' and `av'.
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
%   1) The following example program uses cspice_ckgpav to get C-matrix
%      and associated angular velocity vector for an image whose
%      SCLK count (un-encoded character string version) is known.
%
%      From that matrix and angular velocity vector, the associated
%      state transformation matrix is obtained.
%
%      Note that we need to load a SCLK kernel to convert from clock
%      string to "ticks." Although not required for older spacecraft
%      clocks, most modern spacecraft ones require a leapseconds
%      kernel to be loaded in addition to a SCLK kernel.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: rav2xf_ex1.tm
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
%            File name              Contents
%            --------------------   -----------------------
%            cas00071.tsc           CASSINI SCLK
%            04161_04164ra.bc       CASSINI spacecraft
%                                   reconstructed CK
%
%         \begindata
%
%           KERNELS_TO_LOAD = ( 'cas00071.tsc'
%                               '04161_04164ra.bc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function rav2xf_ex1()
%
%         %
%         % Constants for this program.
%         %
%         % -- The code for the CASSINI spacecraft clock is -82.
%         %
%         % -- The code for CASSINI spacecraft reference frame is
%         %    -82000.
%         %
%         % --  Spacecraft clock tolerance is 1.0 seconds. This may
%         %    not be an acceptable tolerance for some applications.
%         %    It must be converted to "ticks" (units of encoded
%         %    SCLK) for input to cspice_ckgpav.
%         %
%         % -- The reference frame we want is J2000.
%         %
%         META   =   'rav2xf_ex1.tm';
%         REFFRM =   'J2000';
%         SCLKCH =   '1/1465476046.160';
%         SCLTOL =   '1.0';
%         SCID   =   -82;
%         INSTID =   -82000;
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert tolerance from CASSINI formatted character
%         % string SCLK to ticks which are units of encoded SCLK.
%         %
%         [toltik] = cspice_sctiks( SCID, SCLTOL );
%
%         %
%         % cspice_ckgpav requires encoded spacecraft clock.
%         %
%         [sclkdp] = cspice_scencd( SCID, SCLKCH );
%
%         [cmat, av, clkout, found] = cspice_ckgpav( INSTID, sclkdp,       ...
%                                                    toltik, REFFRM  );
%
%         %
%         % Recall that `cmat' and `av' are the rotation and angular
%         % velocity of the transformation from J2000 to the
%         % spacecraft frame.
%         %
%         if ( found )
%
%            %
%            % Display `cmat' and `av'.
%            %
%            fprintf( 'Rotation matrix:\n' )
%            fprintf( '%10.6f %9.6f %9.6f\n', cmat' )
%
%            fprintf( 'Angular velocity:\n' )
%            fprintf( '%20.16f %19.16f %19.16f\n', av )
%
%            %
%            % Get state transformation from J2000 to the spacecraft
%            % frame.
%            %
%            [fxmat] = cspice_rav2xf( cmat, av );
%
%            %
%            % Display the results.
%            %
%            fprintf( '\n' )
%            fprintf( 'State transformation matrix:\n' )
%            fprintf( '%10.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n', fxmat' )
%
%         else
%
%               fprintf( [ 'No rotation matrix/angular velocity',          ...
%                              ' found for %s\n' ], SCLKCH            )
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
%      Rotation matrix:
%       -0.604984  0.796222 -0.005028
%       -0.784160 -0.596891 -0.169748
%       -0.138158 -0.098752  0.985475
%      Angular velocity:
%        0.0000032866819065 -0.0000099372638338  0.0000197597699770
%
%      State transformation matrix:
%       -0.604984  0.796222 -0.005028  0.000000  0.000000  0.000000
%       -0.784160 -0.596891 -0.169748  0.000000  0.000000  0.000000
%       -0.138158 -0.098752  0.985475  0.000000  0.000000  0.000000
%       -0.000016 -0.000012 -0.000003 -0.604984  0.796222 -0.005028
%        0.000013 -0.000015 -0.000010 -0.784160 -0.596891 -0.169748
%       -0.000008 -0.000006 -0.000002 -0.138158 -0.098752  0.985475
%
%
%   2) Compute a state transformation matrix from a rotation matrix
%      for "elementary" frame rotations of 90 degrees about the Z axis
%      and an angular velocity vector, convert that transformation
%      matrix back to a rotation matrix and an angular velocity vector
%      and compute the maximum value of the absolute difference between
%      the rotation matrices and the angular velocity vectors.
%
%      Numerical equivalence shall be expected.
%
%
%      Example code begins here.
%
%
%      function rav2xf_ex2()
%
%         %
%         % Define an angular velocity vector:
%         %
%         e1     =  [ 1.;   0.;  0. ];
%
%         %
%         % Rotation matrix for "elementary" frame rotations:  90 degrees
%         % about the z axis:
%         %
%         rz_90 = [[ 0.,  1.,  0. ];                                       ...
%                  [-1.,  0.,  0. ];                                       ...
%                  [ 0.,  0.,  1. ] ];
%
%         %
%         % The call cspice_rav2xf calculates the state transformation matrix
%         % `strans' associated with the angular velocity vector and the
%         % rotation matrix.
%         %
%         [strans] = cspice_rav2xf( rz_90, e1 );
%
%         %
%         % cspice_xf2rav converts a state transformation to the associated
%         % rotation matrix and angular velocity vector - inverting
%         % the operation of cspice_rav2xf
%         %
%         [rot, av] = cspice_xf2rav( strans );
%
%         %
%         % Calculate the maximum value of the absolute difference between the
%         % output `av' and `rot' vs the inputs `e1' and `rz_90'.
%         %
%         fprintf(['Maximum absolute difference ',                         ...
%                  'between rotation matrices : %15.13f\n'],               ...
%                       max( max( abs(rot - rz_90) ) )  )
%         fprintf(['Maximum absolute difference ',                         ...
%                  'between angular velocities: %15.13f\n'],               ...
%                       max( max(av - e1 ) )            )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Maximum absolute difference between rotation matrices : 0.0000000000000
%      Maximum absolute difference between angular velocities: 0.0000000000000
%
%
%   3) Obtain the state transformation matrix from J2000 to IAU_MOON for
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
%         File name: rav2xf_ex3.tm
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
%      function rav2xf_ex3()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'rav2xf_ex3.tm')
%
%         %
%         % Create an array of 10001 ephemeris times based at July 1 2007.
%         %
%         et    = [0: 10000]* cspice_spd + cspice_str2et( 'July 1 2007' );
%
%         %
%         % Calculate the state transformation matrices from J2000 to IAU_MOON
%         % for `et'.
%         %
%         xform = cspice_sxform( 'J2000', 'IAU_MOON', et );
%
%         %
%         % Convert the set of `xform' matrices to the corresponding rotation
%         % matrices and angular velocity vectors.
%         %
%         [rot, av] = cspice_xf2rav( xform );
%
%         %
%         % Use the converted outputs from cspice_xf2rav to recompute a set
%         % of state transformation matrices.
%         %
%         [strans] = cspice_rav2xf( rot, av );
%
%         %
%         % Calculate the maximum value of the absolute difference between
%         % `xform' and `strans'.
%         %
%         fprintf('Maximum absolute difference: %8.6e\n',                  ...
%                  max( max( max( abs(strans - xform) ) ) ) )
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
%   a rotation and angular velocity of the rotation to the
%   equivalent state transformation matrix.
%
%   This routine is an inverse of cspice_xf2rav.
%
%-Exceptions
%
%   1)  No checks are performed on `rot' to ensure that it is indeed
%       a rotation matrix.
%
%   2)  If any of the input arguments, `rot' or `av', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `rot' or `av', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   4)  If the input vectorizable arguments `rot' and `av' do not have
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added first
%       complete example, examples #2 and #3' problem statement, and
%       meta-kernel for code example #2.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       extended -Particulars section. Added rotation.req to required
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

function [xform] = cspice_rav2xf(rot, av)

   switch nargin
      case 2

         rot = zzmice_dp(rot);
         av  = zzmice_dp(av);

      otherwise

         error ( ['Usage: [_xform(6,6)_] = '                               ...
                  'cspice_rav2xf(_rot(3,3)_, _av(3)_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [xform] = mice('rav2xf_c', rot, av);
   catch spiceerr
      rethrow(spiceerr)
   end
