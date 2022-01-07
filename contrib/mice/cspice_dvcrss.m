%-Abstract
%
%   CSPICE_DVCRSS calculates the cross product of the position components of
%   two state vectors and the time derivative of this cross product.
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
%      s1       a SPICE state(s):
%
%                  s1 = ( r1, dr1/dt )
%
%               [6,n] = size(s1); double = class(s1)
%
%               Typically, this might represent the apparent state of a
%               planet or the Sun, which defines the orientation of axes of
%               some coordinate system.
%
%      s2       a second SPICE state(s):
%
%                  s2 = ( r2, dr2/dt )
%
%               [6,n] = size(s2); double = class(s2)
%
%               An implicit assumption exists that `s1' and `s2' are specified
%               in the same reference frame. If this is not the case, the
%               numerical result has no meaning.
%
%   the call:
%
%      [sout] = cspice_dvcrss( s1, s2 )
%
%   returns:
%
%      sout     the cross product(s) associated with the position components
%               of `s1' and `s2' and the derivative of the cross product(s)
%               with respect to time.
%
%               [6,n] = size(sout); double = class(sout)
%
%               In other words, if s1 = (p1,v1) and s2 = (p2,v2) then
%               `sout' is ( p1xp2, d/dt( p1xp2 ) ).
%
%               `sout' returns with the same vectorization measure (N)
%               as `s1' and `s2'.
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
%   1) Compute the cross product of two 3-dimensional vectors
%      and the derivative of this cross product.
%
%
%      Example code begins here.
%
%
%      function dvcrss_ex1()
%
%         %
%         % Set `s1' and `s2' vectors.
%         %
%         s1 = [ [0.0, 1.0, 0.0, 1.0, 0.0, 0.0]',                          ...
%                [5.0, 5.0, 5.0, 1.0, 0.0, 0.0]' ];
%         s2 = [ [ 1.0,  0.0,  0.0, 1.0, 0.0, 0.0]',                       ...
%                [-1.0, -1.0, -1.0, 2.0, 0.0, 0.0]' ];
%
%         %
%         % For each vector `s1' and `s2', compute their cross product
%         % and its derivative.
%         %
%         for i=1:2
%
%            [sout] = cspice_dvcrss( s1(:,i), s2(:,i) );
%
%            fprintf( 'S1  : %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n', s1(:,i) )
%            fprintf( 'S2  : %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n', s2(:,i) )
%            fprintf( 'SOUT: %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n', sout    )
%            fprintf( '\n' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      S1  :    0.0    1.0    0.0    1.0    0.0    0.0
%      S2  :    1.0    0.0    0.0    1.0    0.0    0.0
%      SOUT:    0.0    0.0   -1.0    0.0    0.0   -1.0
%
%      S1  :    5.0    5.0    5.0    1.0    0.0    0.0
%      S2  :   -1.0   -1.0   -1.0    2.0    0.0    0.0
%      SOUT:    0.0    0.0    0.0    0.0   11.0  -11.0
%
%
%   2) One can construct non-inertial coordinate frames from apparent
%      positions of objects or defined directions. However, if one
%      wants to convert states in this non-inertial frame to states
%      in an inertial reference frame, the derivatives of the axes of
%      the non-inertial frame are required.
%
%      Define a reference frame with the apparent direction of the
%      Sun as seen from Earth as the primary axis X. Use the Earth
%      pole vector to define with the primary axis the XY plane of
%      the frame, with the primary axis Y pointing in the direction
%      of the pole.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: dvcrss_ex2.tm
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
%            pck00008.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00008.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function dvcrss_ex2()
%
%         %
%         % Load SPK, PCK, and LSK kernels, use a meta kernel for
%         % convenience.
%         %
%         cspice_furnsh( 'dvcrss_ex2.tm' )
%
%         %
%         % Define the earth body-fixed pole vector (z). The pole
%         % has no velocity in the earth fixed frame "IAU_EARTH."
%         %
%         z_earth = [ 0, 0, 1, 0, 0, 0 ]';
%
%         %
%         % Calculate the state transformation between IAU_EARTH and J2000
%         % at an arbitrary epoch.
%         %
%         utc     = 'Jan 1, 2009';
%         et      = cspice_str2et( utc );
%         trans   = cspice_sxform( 'IAU_EARTH', 'J2000', et );
%
%         %
%         % Transform the earth pole vector from the IAU_EARTH frame to J2000.
%         %
%         z_j2000 = trans * z_earth;
%
%         %
%         % Calculate the apparent state of the sun from earth at the epoch
%         % 'et' in the J2000 frame.
%         %
%         target   = 'Sun';
%         observer = 'Earth';
%
%         [state, lt] = cspice_spkezr( target, et, 'J2000', 'LT+S', ...
%                                                           observer );
%
%         %
%         % Define the X axis of the new frame to aligned with
%         % the computed state. Calculate the state's unit vector
%         % and its derivative to get the X axis and its
%         % derivative.
%         %
%         x_new = cspice_dvhat( state );
%
%         %
%         % Define the Z axis of the new frame as the cross product
%         % between the computed state and the Earth pole.
%         % Calculate the Z direction in the new reference frame,
%         % then calculate the this direction's unit vector and its
%         % derivative to get the Z axis and its derivative.
%         %
%         z_new = cspice_dvcrss( state, z_j2000 );
%         z_new = cspice_dvhat( z_new );
%
%         %
%         % As for `z_new', calculate the y direction in the new
%         % reference frame, then calculate this direction's unit
%         % vector and its derivative to get the Y axis and its
%         % derivative.
%         %
%         y_new = cspice_dvcrss( z_new, state );
%         y_new = cspice_dvhat( y_new );
%
%         fprintf('New X-axis:\n' );
%         fprintf('   position: %15.12f %15.12f %15.12f\n',   x_new(1:3) );
%         fprintf('   velocity: %15.12f %15.12f %15.12f\n\n', x_new(4:6) );
%         fprintf('New Y-axis:\n' );
%         fprintf('   position: %15.12f %15.12f %15.12f\n',   y_new(1:3) );
%         fprintf('   velocity: %15.12f %15.12f %15.12f\n\n', y_new(4:6) );
%         fprintf('New Z-axis:\n' );
%         fprintf('   position: %15.12f %15.12f %15.12f\n',   z_new(1:3) );
%         fprintf('   velocity: %15.12f %15.12f %15.12f\n\n', z_new(4:6) );
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
%      New X-axis:
%         position:  0.183446637633 -0.901919663328 -0.391009273602
%         velocity:  0.000000202450  0.000000034660  0.000000015033
%
%      New Y-axis:
%         position:  0.078846540163 -0.382978080242  0.920386339077
%         velocity:  0.000000082384  0.000000032309  0.000000006387
%
%      New Z-axis:
%         position: -0.979862518033 -0.199671507623  0.000857203851
%         velocity:  0.000000044531 -0.000000218531 -0.000000000036
%
%
%      Note that these vectors define the transformation between the
%      new frame and J2000 at the given `et':
%
%             .-            -.
%             |       :      |
%             |   R   :  0   |
%         M = | ......:......|
%             |       :      |
%             | dRdt  :  R   |
%             |       :      |
%             `-            -'
%
%      with
%
%         R    = [ x_new(1:3); y_new(1:3); z_new(1:3) ]
%
%         dRdt = [ x_new(4:6); y_new(4:6); z_new(4:6) ]
%
%-Particulars
%
%   cspice_dvcrss calculates the three-dimensional cross product of two
%   vectors and the derivative of that cross product according to
%   the definition.
%
%   In this discussion, the notation
%
%      V1 x V2
%
%   indicates the cross product of vectors V1 and V2.
%
%   With s1 = (r1,v1) and s2 = (r2,v2) then
%
%                         d
%      sout = [ r1 x r2 , -- (r1 x r2) ]
%                         dt
%
%-Exceptions
%
%   1)  If `s1' and `s2' are large in magnitude (taken together,
%       their magnitude surpasses the limit allowed by the
%       computer) then it may be possible to generate a
%       floating point overflow from an intermediate
%       computation even though the actual cross product and
%       derivative may be well within the range of double
%       precision numbers.
%
%       cspice_dvcrss does NOT check the magnitude of `s1' or `s2' to
%       insure that overflow will not occur.
%
%   2)  If any of the input arguments, `s1' or `s2', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `s1' or `s2', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   4)  If the input vectorizable arguments `s1' and `s2' do not have
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
%   -Mice Version 1.1.0, 09-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "dvcrss" to "sout".
%
%       Edited the header to comply with NAIF standard. Added example's
%       problem statement and reference to required meta-kernel.
%       Reformatted example's output. Added first example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 09-NOV-2010 (EDW)
%
%-Index_Entries
%
%   Compute the derivative of a cross product
%
%-&

function [sout] = cspice_dvcrss(s1, s2)

   switch nargin
      case 2

         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);

      otherwise

         error ( 'Usage: [_sout(6)_] = cspice_dvcrss(_s1(6)_, _s2(6)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [sout] = mice( 'dvcrss_c', s1, s2);
   catch spiceerr
      rethrow(spiceerr)
   end



