%-Abstract
%
%   CSPICE_EUL2XF computes a state transformation from an Euler angle
%   factorization of a rotation and the derivatives of those Euler angles.
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
%      eulang   the array(s) of Euler angles corresponding to the
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
%               then (6x1)
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
%      axisa,
%      axisb,
%      axisc    the indices defining the axes desired for the factorization
%               of `r'.
%
%               [1,1] = size(axis3); int32 = class(axis3)
%               [1,1] = size(axis2); int32 = class(axis2)
%               [1,1] = size(axis1); int32 = class(axis1)
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
%      [xform] = cspice_eul2xf( eulang, axisa, axisb, axisc )
%
%   returns:
%
%      xform    the state transformation matri(x|ces) corresponding to `r'
%               and dr/dt as described above.
%
%               Either [3,3]   = size(xform); double = class(xform)
%               or     [3,3,n] = size(xform); double = class(xform)
%
%               Pictorially,
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
%               `xform' returns with the same vectorization
%               measure, N, as `eulang'.
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
%   1) Suppose you have a set of Euler angles and their derivatives
%      for a 3 1 3 rotation, and that you would like to determine
%      the equivalent angles and derivatives for a 1 2 3 rotation.
%
%         r = [alpha]  [beta]  [gamma]
%                    3       1        3
%
%         r = [roll]  [pitch]  [yaw]
%                   1        2      3
%
%      The following code example will perform the desired
%      computation.
%
%
%      Example code begins here.
%
%
%      function eul2xf_ex1()
%
%         %
%         % Define the initial set of Euler angles.
%         %
%         abgang = [0.01; 0.03; 0.09; -0.001; -0.003; -0.009 ];
%
%         xform            = cspice_eul2xf( abgang, 3, 1, 3 );
%         [rpyang, unique] = cspice_xf2eul( xform , 1, 2, 3 );
%
%         if( unique )
%
%            disp( '1-2-3 equivalent rotation to input (radians):')
%            fprintf( 'Roll   %12.9f, droll/dt   %12.9f\n',                ...
%                                            rpyang(1), rpyang(4) )
%            fprintf( 'Pitch  %12.9f, dPitch/dt  %12.9f\n',                ...
%                                            rpyang(2), rpyang(5) )
%            fprintf( 'Yaw    %12.9f, dYaw/dt    %12.9f\n',                ...
%                                            rpyang(3), rpyang(6) )
%
%         else
%
%            disp( 'The values in ''rpyang'' not uniquely determined.' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      1-2-3 equivalent rotation to input (radians):
%      Roll    0.029998501, droll/dt   -0.002999550
%      Pitch  -0.000299950, dPitch/dt   0.000059980
%      Yaw     0.099995501, dYaw/dt    -0.009998650
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
%   This routine is intended to provide an inverse for cspice_xf2eul.
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
%   2)  If any of the input arguments, `eulang', `axisa', `axisb' or
%       `axisc', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   3)  If any of the input arguments, `eulang', `axisa', `axisb' or
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
%       Edited the header to comply with NAIF standard. Added example's
%       problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section. Added rotation.req required readings.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 06-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.2, 29-FEB-2012 (EDW)
%
%       Edit to "Usage" string. "xform(3,3)" corrected to read
%       "xform(6,6)."
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 02-APR-2007 (EDW)
%
%-Index_Entries
%
%   State transformation from Euler angles and derivatives
%
%-&

function [xform] = cspice_eul2xf(eulang, axisa, axisb, axisc)

   switch nargin
      case 4

         eulang= zzmice_dp(eulang);
         axisa = zzmice_int( axisa);
         axisb = zzmice_int( axisb);
         axisc = zzmice_int( axisc);

      otherwise

         error ( ['Usage: [_xform(6,6)_] = ' ...
                 'cspice_eul2xf(_eulang(6)_, axisa, axisb, axisc)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [xform] = mice('eul2xf_c', eulang, axisa, axisb, axisc);
   catch spiceerr
      rethrow(spiceerr)
   end



