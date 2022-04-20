%-Abstract
%
%   CSPICE_XF2EUL converts a state transformation matrix to Euler angles and
%   their derivatives, given a specified set of axes.
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
%      xform    a state transformation matri(x|ces) from some frame FRAME1
%               to another frame FRAME2.
%
%               Either [6,6]   = size(xform); double = class(xform)
%               or     [6,6,n] = size(xform); double = class(xform)
%
%               Pictorially, `xform' has the structure shown here.
%
%                  .-             -.
%                  |       |       |
%                  |   r   |   0   |
%                  |       |       |
%                  |-------+-------|
%                  |       |       |
%                  | dr/dt |   r   |
%                  |       |       |
%                  `-             -'
%
%               where `r' is a rotation matrix that varies with respect to
%               time and dr/dt is its time derivative.
%
%               More specifically, if `s1' is the state of some object
%               in FRAME1, then `s2', the state of the same object
%               relative to FRAME2 is given by
%
%                  s2 = xform * s1
%
%               where "*" denotes the matrix vector product.
%
%      axisa,
%      axisb,
%      axisc    the axes desired for the factorization of `r'.
%
%               [1,1] = size(axisa); int32 = class(axisa)
%               [1,1] = size(axisb); int32 = class(axisb)
%               [1,1] = size(axisc); int32 = class(axisc)
%
%               All must be in the range from 1 to 3. Moreover
%               it must be the case that `axisa' and `axisb' are distinct
%               and that `axisb' and `axisc' are distinct.
%
%               Every rotation matrix can be represented as a product
%               of three rotation matrices about the principal axes
%               of a reference frame.
%
%                  r =  [ alpha ]      [ beta ]      [ gamma ]
%                                axisa         axisb          axisc
%
%               The value 1 corresponds to the X axis.
%               The value 2 corresponds to the Y axis.
%               The value 3 corresponds to the Z axis.
%
%   the call:
%
%      [eulang, unique] = cspice_xf2eul( xform, axisa, axisb, axisc )
%
%   returns:
%
%      eulang   the vector(s) of Euler angles corresponding to the
%               specified factorization.
%
%               Either [6,1] = size(eulang); double = class(eulang)
%               or     [6,n] = size(eulang); double = class(eulang)
%
%               If we represent `r' as shown here:
%
%                  r =  [ alpha ]      [ beta ]      [ gamma ]
%                                axisa         axisb          axisc
%
%               then
%
%                  eulang(1) = alpha
%                  eulang(2) = beta
%                  eulang(3) = gamma
%                  eulang(4) = dalpha/dt
%                  eulang(5) = dbeta/dt
%                  eulang(6) = dgamma/dt
%
%               or (6xN)
%
%                  eulang(1,N) = alpha_N
%                  eulang(2,N) = beta_N
%                  eulang(3,N) = gamma_N
%                  eulang(4,N) = dalpha_N/dt
%                  eulang(5,N) = dbeta_N/dt
%                  eulang(6,N) = dgamma_N/dt
%
%               The range of `alpha' and `gamma' is (-pi, pi].
%
%               The range of `beta' depends on the exact set of
%               axes used for the factorization. For
%               factorizations in which the first and third axes
%               are the same, the range of `beta' is [0, pi].
%
%               For factorizations in which the first and third
%               axes are different, the range of `beta' is
%               [-pi/2, pi/2].
%
%               For rotations such that `alpha' and `gamma' are not
%               uniquely determined, `alpha' and dalpha/dt will
%               always be set to zero; `gamma' and dgamma/dt are
%               then uniquely determined.
%
%      unique   a logical flag(s) that indicates whether or not the values in
%               `eulang' are uniquely determined.
%
%               Either [1,1] = size(unique); logical = class(unique)
%               or     [1,n] = size(unique); logical = class(unique)
%
%               If the values are unique then `unique' will be set to true.
%               If the values are not unique and some components ( eulang(1)
%               and eulang(4) ) have been set to zero, then `unique' will
%               have the value false.
%
%               `eulang' and `unique' return with the same vectorization
%               measure, N, as `xform'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Determine the rate of change of the right ascension and
%      declination of the pole of the moon, from the state
%      transformation matrix that transforms J2000 states to object
%      fixed states.
%
%      Recall that the rotation component of the state transformation
%      matrix is given by
%
%         [w]  [halfpi_c-dec] [ra+halfpi_c]
%            3               1             3
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: xf2eul_ex1.tm
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
%      function xf2eul_ex1()
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'xf2eul_ex1.tm' )
%
%         %
%         % Define the number of ephemeris times to perform the calculation.
%         %
%         N = 6;
%
%         %
%         % Calculate the separation of each ephemeris time, in seconds,
%         % over an eighteen year span.
%         %
%         STEP = 18 * 365 * cspice_spd/N;
%
%         %
%         % Base the ephemeris time set at May 15, 2007.
%         %
%         et = [0:N]*STEP +  cspice_str2et( 'May 15, 2007' );
%
%         %
%         % Calculate the state transformation matrices corresponding
%         % to 'et', then convert those matrices to Euler angles (3-1-3).
%         %
%         tsipm              = cspice_sxform( 'J2000', 'IAU_MOON', et );
%         [ eulang, unique ] = cspice_xf2eul( tsipm , 3, 1, 3 );
%
%         %
%         % From the Euler angles, calculate right ascension and declination.
%         % Form the UTC time string from 'et' (for output purposes).
%         %
%         ra  = eulang(3,:) - cspice_halfpi;
%         dec = cspice_halfpi - eulang(2,:);
%         utc = cspice_et2utc( et, 'c', 3 );
%
%         %
%         % As a convenience, output in a loop.
%         %
%         for m=1:N+1
%
%            if( unique(m) )
%               fprintf( 'UTC: %s\n', utc(m,:)                 )
%               fprintf( 'w        = %12.6f\n'  , eulang(1,m)  )
%               fprintf( 'dec      = %12.6f\n'  , dec(m)       )
%               fprintf( 'ra       = %12.6f\n'  , ra(m)        )
%               fprintf( 'd w/dt   = %14.9f\n'  , eulang(4,m)  )
%               fprintf( 'd dec/dt = %14.9f\n'  ,-eulang(5,m)  )
%               fprintf( 'd ra/dt  = %14.9f\n\n', eulang(6,m)  )
%            else
%               disp( 'The values in ''eulang'' not uniquely determined.' )
%            end
%
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
%      UTC: 2007 MAY 15 00:00:00.000
%      w        =    -2.649088
%      dec      =     1.186911
%      ra       =    -1.549644
%      d w/dt   =    0.000002658
%      d dec/dt =   -0.000000000
%      d ra/dt  =    0.000000004
%
%      UTC: 2010 MAY 13 23:59:59.000
%      w        =    -2.199870
%      dec      =     1.168228
%      ra       =    -1.504392
%      d w/dt   =    0.000002665
%      d dec/dt =   -0.000000000
%      d ra/dt  =   -0.000000004
%
%      UTC: 2013 MAY 12 23:59:58.000
%      w        =    -1.691639
%      dec      =     1.142338
%      ra       =    -1.523436
%      d w/dt   =    0.000002659
%      d dec/dt =    0.000000000
%      d ra/dt  =    0.000000003
%
%      UTC: 2016 MAY 11 23:59:57.000
%      w        =    -1.145062
%      dec      =     1.135591
%      ra       =    -1.584588
%      d w/dt   =    0.000002663
%      d dec/dt =   -0.000000001
%      d ra/dt  =   -0.000000002
%
%      UTC: 2019 MAY 11 23:59:56.000
%      w        =    -0.611122
%      dec      =     1.151760
%      ra       =    -1.631740
%      d w/dt   =    0.000002662
%      d dec/dt =    0.000000002
%      d ra/dt  =   -0.000000001
%
%      UTC: 2022 MAY 10 23:59:56.000
%      w        =    -0.123811
%      dec      =     1.177463
%      ra       =    -1.627743
%      d w/dt   =    0.000002660
%      d dec/dt =   -0.000000001
%      d ra/dt  =    0.000000002
%
%      UTC: 2025 MAY 09 23:59:56.000
%      w        =     0.307001
%      dec      =     1.188934
%      ra       =    -1.562882
%      d w/dt   =    0.000002663
%      d dec/dt =    0.000000001
%      d ra/dt  =   -0.000000001
%
%
%-Particulars
%
%   A word about notation: the symbol
%
%      [ x ]
%           i
%
%   indicates a coordinate system rotation of x radians about the
%   ith coordinate axis. To be specific, the symbol
%
%      [ x ]
%           1
%
%   indicates a coordinate system rotation of x radians about the
%   first, or x-, axis; the corresponding matrix is
%
%      .-                    -.
%      |  1    0        0     |
%      |                      |
%      |  0    cos(x)  sin(x) |
%      |                      |
%      |  0   -sin(x)  cos(x) |
%      `-                    -'
%
%   Remember, this is a COORDINATE SYSTEM rotation by x radians; this
%   matrix, when applied to a vector, rotates the vector by -x
%   radians, not x radians. Applying the matrix to a vector yields
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
%      |  0       1    0      |
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
%   The input matrix is assumed to be the product of three
%   rotation matrices, each one of the form
%
%      .-                    -.
%      |  1      0       0    |
%      |                      |
%      |  0    cos(r)  sin(r) |     (rotation of r radians about the
%      |                      |      x-axis),
%      |  0   -sin(r)  cos(r) |
%      `-                    -'
%
%
%      .-                    -.
%      | cos(s)   0   -sin(s) |
%      |                      |
%      |  0       1      0    |     (rotation of s radians about the
%      |                      |      y-axis),
%      | sin(s)   0    cos(s) |
%      `-                    -'
%
%   or
%
%      .-                    -.
%      |  cos(t)  sin(t)   0  |
%      |                      |
%      | -sin(t)  cos(t)   0  |     (rotation of t radians about the
%      |                      |      z-axis),
%      |  0        0       1  |
%      `-                    -'
%
%   where the second rotation axis is not equal to the first or
%   third. Any rotation matrix can be factored as a sequence of
%   three such rotations, provided that this last criterion is met.
%
%   This routine is related to the routine cspice_eul2xf which produces
%   a state transformation from an input set of axes, Euler angles
%   and derivatives.
%
%   The two function calls shown here will not change
%   `xform' except for round off errors.
%
%      [eulang, unique] = cspice_xf2eul( xform,  axisa, axisb, axisc );
%      [xform]          = cspice_eul2xf( eulang, axisa, axisb, axisc );
%
%   On the other hand the two calls
%
%      [xform]          = cspice_eul2xf( eulang, axisa, axisb, axisc );
%      [eulang, unique] = cspice_xf2eul( xform,  axisa, axisb, axisc );
%
%   will leave `eulang' unchanged only if the components of `eulang'
%   are in the range produced by cspice_xf2eul and the Euler representation
%   of the rotation component of `xform' is unique within that range.
%
%-Exceptions
%
%   1)  If any of `axisa', `axisb', or `axisc' do not have values in
%
%          { 1, 2, 3 }
%
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If `axisb' is equal to `axisc' or `axisa', an error is signaled by a
%       routine in the call tree of this routine. An arbitrary
%       rotation matrix cannot be expressed using a sequence of Euler
%       angles unless the second rotation axis differs from the other
%       two.
%
%   3)  If the input matrix `xform' is not a rotation matrix, an error
%       is signaled by a routine in the call tree of this routine.
%
%   4)  If eulang(1) and eulang(3) are not uniquely determined,
%       eulang(1) is set to zero, and eulang(3) is determined.
%
%   5)  If any of the input arguments, `xform', `axisa', `axisb' or
%       `axisc', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   6)  If any of the input arguments, `xform', `axisa', `axisb' or
%       `axisc', is not of the expected type, or it does not have the
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
%   PCK.REQ
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
%   -Mice Version 1.3.0, 01-NOV-2021 (EDW) (JDR)
%
%       Edited header to comply with NAIF standard. Reduced the number of
%       ephemeris epochs used as input in code example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section. Added pck.req and rotation.req
%       required readings.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.2.1, 19-SEP-2016 (EDW)
%
%       Corrected usage string typo, 'xform' should show dimension as (6,6).
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.2.0, 10-MAY-2011 (EDW)
%
%       "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 02-APR-2007 (EDW)
%
%-Index_Entries
%
%   Euler angles and derivatives from state transformation
%
%-&

function [eulang, unique] = cspice_xf2eul(xform,  axisa, axisb, axisc)

   switch nargin
      case 4

         xform = zzmice_dp( xform);
         axisa = zzmice_int( axisa);
         axisb = zzmice_int( axisb);
         axisc = zzmice_int( axisc);

      otherwise

         error ( [ 'Usage: [_eulang(6)_, _unique_] = ' ...
                   'cspice_xf2eul(_xform(6,6)_, axisa, axisb, axisc)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [eulang, unique] = mice('xf2eul_c', xform, axisa, axisb, axisc);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      unique = zzmice_logical(unique);
   catch spiceerr
      rethrow(spiceerr)
   end


