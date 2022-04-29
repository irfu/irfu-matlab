%-Abstract
%
%   CSPICE_NVC2PL constructs a SPICE plane from a normal vector
%   and a constant.
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
%      normal,
%      konst    respectively, a normal vector and constant defining a
%               plane.
%
%               [3,1] = size(normal); double = class(normal)
%               [1,1] = size(konst);  double = class(konst)
%
%               `normal' need not be a unit vector.
%
%               Let the symbol < a, b > indicate the inner product of
%               vectors a and b; then the geometric plane is the set of
%               vectors x in three-dimensional space that satisfy
%
%                    < x,  normal >  =  konst.
%
%   the call:
%
%      plane = cspice_nvc2pl( normal, konst )
%
%   returns:
%
%      plane    a structure describing a SPICE plane defined by
%               `normal' and `konst'
%
%               [1,1] = size(plane); struct = class(plane)
%
%               The structure has the fields:
%
%                  normal:   [3,1] = size(normal);   double = class(normal)
%                  constant: [1,1] = size(constant); double = class(constant)
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
%   1) Construct a SPICE plane from a normal vector and a constant.
%
%      Example code begins here.
%
%
%      function nvc2pl_ex1()
%
%         %
%         % Define an arbitrary normal and constant...
%         %
%         normal    = [ 1.; 1.; 1. ];
%         konst  = 23.;
%         fprintf( 'Inputs:\n' );
%         fprintf( '  Normal vector: %15.12f %15.12f %15.12f\n', ...
%                                                  normal       )
%         fprintf( '  Constant     : %15.12f\n\n', konst        )
%
%         %
%         % ...then construct the plane.
%         %
%         plane = cspice_nvc2pl( normal, konst );
%
%         fprintf( 'Generated plane:\n' )
%         fprintf( '  Normal vector: %15.12f %15.12f %15.12f\n', ...
%                                                       plane.normal   )
%         fprintf( '  Constant     : %15.12f\n\n',      plane.constant )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Inputs:
%        Normal vector:  1.000000000000  1.000000000000  1.000000000000
%        Constant     : 23.000000000000
%
%      Generated plane:
%        Normal vector:  0.577350269190  0.577350269190  0.577350269190
%        Constant     : 13.279056191361
%
%
%   2) Apply a linear transformation represented by a matrix to
%      a plane represented by a normal vector and a constant.
%
%      Find a normal vector and constant for the transformed plane.
%
%
%      Example code begins here.
%
%
%      function nvc2pl_ex2()
%
%         %
%         % Set the normal vector and the constant defining the
%         % initial plane.
%         %
%         normal = [-0.1616904, 0.8084521, -0.5659165]';
%         konst  =   4.8102899;
%
%         %
%         % Define a transformation matrix to the right-handed
%         % reference frame having the +i unit vector as primary
%         % axis, aligned to the original frame's +X axis, and
%         % the -j unit vector as second axis, aligned to the +Y
%         % axis.
%         %
%         axdef  = [1.0,  0.0,  0.0]';
%         plndef = [0.0, -1.0,  0.0]';
%
%         [m]    = cspice_twovec( axdef, 1, plndef, 2 );
%
%         %
%         % Make a SPICE plane from `normal' and `konst', and then
%         % find a point in the plane and spanning vectors for the
%         % plane.  `normal' need not be a unit vector.
%         %
%         [plane]               = cspice_nvc2pl( normal, konst );
%         [point, span1, span2] = cspice_pl2psv( plane );
%
%         %
%         % Apply the linear transformation to the point and
%         % spanning vectors.  All we need to do is multiply
%         % these vectors by `m', since for any linear
%         % transformation T,
%         %
%         %       T ( point  +  t1 * span1     +  t2 * span2 )
%         %
%         %    =  T (point)  +  t1 * T(span1)  +  t2 * T(span2),
%         %
%         % which means that T(point), T(span1), and T(span2)
%         % are a point and spanning vectors for the transformed
%         % plane.
%         %
%         tpoint = m * point;
%         tspan1 = m * span1;
%         tspan2 = m * span2;
%
%         %
%         % Make a new SPICE plane `tplane' from the
%         % transformed point and spanning vectors, and find a
%         % unit normal and constant for this new plane.
%         %
%         [tplane]         = cspice_psv2pl( tpoint, tspan1, tspan2 );
%         [tnorml, tkonst] = cspice_pl2nvc( tplane );
%
%         %
%         % Print the results.
%         %
%         fprintf( 'Unit normal vector: %11.7f %11.7f %11.7f\n',           ...
%                                tnorml(1), tnorml(2), tnorml(3) )
%         fprintf( 'Constant          : %11.7f\n', tkonst )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Unit normal vector:  -0.1616904  -0.8084521   0.5659165
%      Constant          :   4.8102897
%
%
%-Particulars
%
%   Mice geometry routines that deal with planes use the `plane'
%   data type to represent input and output planes. This data type
%   makes the subroutine interfaces simpler and more uniform.
%
%   The Mice routines that produce SPICE planes from data that
%   define a plane are:
%
%      cspice_nvc2pl ( Normal vector and constant to plane )
%      cspice_nvp2pl ( Normal vector and point to plane    )
%      cspice_psv2pl ( Point and spanning vectors to plane )
%
%   The Mice routines that convert SPICE planes to data that
%   define a plane are:
%
%      cspice_pl2nvc ( Plane to normal vector and constant )
%      cspice_pl2nvp ( Plane to normal vector and point    )
%      cspice_pl2psv ( Plane to point and spanning vectors )
%
%   Any of these last three routines may be used to convert this
%   routine's output, 'plane', to another representation of a
%   geometric plane.
%
%-Exceptions
%
%   1)  If the input vector `normal' is the zero vector, the error
%       SPICE(ZEROVECTOR) is signaled by a routine in the call tree of
%       this routine.
%
%   2)  If any of the input arguments, `normal' or `konst', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `normal' or `konst', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  No checking is done to prevent arithmetic overflow.
%
%-Required_Reading
%
%   MICE.REQ
%   PLANES.REQ
%
%-Literature_References
%
%   [1]  G. Thomas and R. Finney, "Calculus and Analytic Geometry,"
%        7th Edition, Addison Wesley, 1988.
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
%       Changed the argument name "constant" to "konst" for consistency
%       with other routines.
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement, modified code example to produce
%       formatted output and added second example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 27-AUG-2012 (EDW)
%
%      Edited -I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.0, 30-DEC-2008 (EDW)
%
%-Index_Entries
%
%   normal vector and constant to plane
%
%-&

function [plane] = cspice_nvc2pl( normal, konst )

   switch nargin

      case 2

         normal = zzmice_dp(normal);
         konst  = zzmice_dp(konst);

      otherwise

         error ( ['Usage: [plane] = cspice_nvc2pl( normal(3), konst )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [plane] = mice('nvc2pl_c', normal, konst );
   catch spiceerr
      rethrow(spiceerr)
   end
