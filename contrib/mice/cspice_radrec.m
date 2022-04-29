%-Abstract
%
%   CSPICE_RADREC converts the right ascension, declination
%   coordinates of a location to rectangular (Cartesian)
%   coordinates.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      range    the value(s) describing the distance of the position
%               from the origin.
%
%               [1,n] = size(range); double = class(range)
%
%               Input should be in terms of the same units in which the
%               output is desired.
%
%      ra       the value(s) describing the right ascension of the
%               right ascension of the position.
%
%               [1,n] = size(ra); double = class(ra)
%
%               This is the angular distance measured toward the east from
%               the prime meridian to the meridian containing the input
%               point. The direction of increasing right ascension is from
%               the +X axis towards the +Y axis.
%
%               The range (i.e., the set of allowed values) of
%               `ra' is unrestricted. Units are radians.
%
%      dec      the value(s) describing the declination of the position.
%
%               [1,n] = size(dec); double = class(dec)
%
%               This is the angle from the XY plane of the ray from the
%               origin through the point.
%
%               The range (i.e., the set of allowed values) of
%               `dec' is unrestricted. Units are radians.
%
%   the call:
%
%      [rectan] = cspice_radrec( range, ra, dec )
%
%   returns:
%
%      rectan   the array(s) containing the rectangular coordinates of the
%               position(s).
%
%               [3,n] = size(rectan); double = class(rectan)
%
%               `rectan' returns with the same units associated with `range'.
%
%               `rectan' returns with the same vectorization measure, N,
%                as `range', `ra', and `dec'.
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
%   1) Convert to the J2000 frame the right ascension and declination
%      of an object initially expressed with respect to the B1950
%      reference frame.
%
%
%      Example code begins here.
%
%
%      function radrec_ex1()
%
%         %
%         % Set the initial right ascension and declination
%         % coordinates of the object, given with respect
%         % to the B1950 reference frame.
%         %
%         rab  = 135.88680896;
%         decb =  17.50151037;
%
%         %
%         % Convert `rab' and `decb' to a 3-vector expressed in
%         % the B1950 frame.
%         %
%         [v1950] = cspice_radrec( 1.0, rab * cspice_rpd,                  ...
%                                  decb * cspice_rpd     );
%
%         %
%         % We use the Mice routine cspice_pxform to obtain the
%         % transformation  matrix for converting vectors between
%         % the B1950 and J2000 reference frames.  Since
%         % both frames are inertial, the input time value we
%         % supply to cspice_pxform is arbitrary.  We choose zero
%         % seconds past the J2000 epoch.
%         %
%         [mtrans] = cspice_pxform( 'B1950', 'J2000', 0.0 );
%
%         %
%         % Transform the vector to the J2000 frame.
%         %
%         v2000 = mtrans * v1950;
%
%         %
%         % Find the right ascension and declination of the
%         % J2000-relative vector.
%         %
%         [r, raj, decj] = cspice_recrad( v2000 );
%
%         %
%         % Output the results.
%         %
%         fprintf( 'Right ascension (B1950 frame): %f\n', rab )
%         fprintf( 'Declination (B1950 frame)    : %f\n', decb )
%
%         fprintf( 'Right ascension (J2000 frame): %f\n',                  ...
%                                           raj * cspice_dpr )
%         fprintf( 'Declination (J2000 frame)    : %f\n',                  ...
%                                          decj * cspice_dpr )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Right ascension (B1950 frame): 135.886809
%      Declination (B1950 frame)    : 17.501510
%      Right ascension (J2000 frame): 136.587682
%      Declination (J2000 frame)    : 17.300443
%
%
%   2) Define a set of 15 right ascension-declination data pairs for
%      the Earth's pole at different ephemeris epochs, convert them
%      to rectangular coordinates and compute the angular separation
%      between these coordinates and the IAU_EARTH pole given by a
%      PCK kernel.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for the Earth.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function radrec_ex2()
%
%         %
%         % Load a PCK kernel.
%         %
%         cspice_furnsh( 'pck00010.tpc' )
%
%         %
%         % Define a set of 15 right ascension-declination data sets
%         % pairs (in degrees) for the earth's pole and the array of
%         % corresponding ephemeris times J2000 TDB.
%         %
%         right_ascen = [ 180.003739, 180.003205, 180.002671,              ...
%                         180.002137, 180.001602, 180.001068,              ...
%                         180.000534, 360.000000, 359.999466,              ...
%                         359.998932, 359.998397, 359.997863,              ...
%                         359.997329, 359.996795, 359.996261 ]';
%
%          dec        = [  89.996751,  89.997215,  89.997679,              ...
%                          89.998143,  89.998608,  89.999072,              ...
%                          89.999536,  90.000000,  89.999536,              ...
%                          89.999072,  89.998607,  89.998143,              ...
%                          89.997679,  89.997215,  89.996751 ]';
%
%          et         = [ -18408539.52023917, -15778739.49107254,          ...
%                         -13148939.46190590, -10519139.43273926,          ...
%                         -7889339.40357262,   -5259539.37440598,          ...
%                         -2629739.34523934,         60.68392730,          ...
%                          2629860.71309394,    5259660.74226063,          ...
%                          7889460.77142727,   10519260.80059391,          ...
%                         13149060.82976055,   15778860.85892719,          ...
%                         18408660.88809383 ]';
%
%         %
%         % Create a 1xN array of radii, the length of a
%         % unit vector (1) the same size as the above arrays.
%         %
%         n_elements  = size(dec);
%         rad         = ones( 1,  n_elements(1) );
%         z_hat       = [0; 0; 1];
%
%         %
%         % Convert the RA/DEC values to radians.
%         %
%         right_ascen = right_ascen * cspice_rpd;
%         dec         = dec * cspice_rpd;
%
%         %
%         % Convert the angular description of the unit vectors to
%         % Cartesian.
%         %
%         pole        = cspice_radrec( rad, right_ascen', dec');
%
%         %
%         % Retrieve the transformation matrix from the J2000
%         % frame to the IAU_EARTH frame.
%         %
%         mat         = cspice_pxform( 'J2000', 'IAU_EARTH', et');
%
%         %
%         % Rotate the 'pole' vector set into IAU_EARTH. All vectors
%         % should equal (to round-off) the z direction unit vector.
%         %
%
%         disp( '      et            Angular difference' )
%         disp( '------------------  ------------------' )
%
%         for i =1:15
%            z = mat(:,:,i) * pole(:,i);
%
%            %
%            % Output the ephemeris time, the pole vector, and the angular
%            % separation between the calculated and the expected pole
%            % vectors.
%            %
%            txt = sprintf( '%18.8f  %18.16f',                             ...
%                           et(i), cspice_vsep(z,z_hat) * cspice_dpr() );
%            disp(txt)
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
%            et            Angular difference
%      ------------------  ------------------
%      -18408539.52023917  0.0000001559918278
%      -15778739.49107254  0.0000000106799881
%      -13148939.46190590  0.0000001773517911
%      -10519139.43273926  0.0000003440236194
%       -7889339.40357262  0.0000004893045693
%       -5259539.37440598  0.0000003226327536
%       -2629739.34523934  0.0000001559609507
%             60.68392730  0.0000000107108706
%        2629860.71309394  0.0000001773826862
%        5259660.74226063  0.0000003440544891
%        7889460.77142727  0.0000004892736740
%       10519260.80059391  0.0000003226018712
%       13149060.82976055  0.0000001559300556
%       15778860.85892719  0.0000000107417474
%       18408660.88809383  0.0000001774135760
%
%
%-Particulars
%
%   This routine converts the right ascension, declination, and range
%   of a point into the associated rectangular coordinates.
%
%   The input is defined by a distance from a central reference point,
%   an angle from a reference meridian, and an angle above the equator
%   of a sphere centered at the central reference point.
%
%-Exceptions
%
%   1)  If any of the input arguments, `range', `ra' or `dec', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `range', `ra' or `dec', is not
%       of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%   3)  If the input vectorizable arguments `range', `ra' and `dec' do
%       not have the same measure of vectorization (N), an error is
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
%   [1]  L. Taff, "Celestial Mechanics, A Computational Guide for the
%        Practitioner," Wiley, 1985
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
%       Edited the header to comply with NAIF standard. Added example's
%       problem statement, reformatted code example output and added
%       reference to required PCK. Added an additional complete code
%       example.
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
%   -Mice Version 1.0.2, 07-JAN-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   range ra and dec to rectangular coordinates
%   right_ascension and declination to rectangular
%
%-&

function [rectan] = cspice_radrec( range, ra, dec )

   switch nargin
      case 3

         range = zzmice_dp(range);
         ra    = zzmice_dp(ra);
         dec   = zzmice_dp(dec);

      otherwise
         error ( ['Usage: [_rectan(3)_] = cspice_radrec( _range_, ',  ...
                                                    '_ra_, _dec_ )'] )
   end

   %
   % Call the MEX library.
   %
   try
      [rectan] = mice('radrec_c', range, ra, dec);
   catch spiceerr
      rethrow(spiceerr)
   end

