%-Abstract
%
%   CSPICE_PL2NVC returns a unit normal vector and constant that define a
%   specified plane.
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
%      plane    a SPICE plane.
%
%               [1,1] = size(plane); struct = class(plane)
%
%               The structure has the fields:
%
%                  normal:   [3,1] = size(normal);   double = class(normal)
%                  constant: [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [normal, konst] = cspice_pl2nvc( plane )
%
%   returns:
%
%      normal,
%      konst    respectively, a unit normal vector and constant that
%               define the geometric plane represented by `plane'.
%
%               [3,1] = size(normal); double = class(normal)
%               [1,1] = size(konst);  double = class(konst)
%
%               Let the symbol < a, b > indicate the inner product of
%               vectors `a' and `b'; then the geometric plane is the set of
%               vectors `x' in three-dimensional space that satisfy
%
%                  < x,  normal >  =  konst.
%
%               `normal' is a unit vector. `konst' is the distance of
%               the plane from the origin;
%
%                  konst * normal
%
%               is the closest point in the plane to the origin.
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
%   1) Determine the distance of a plane from the origin, and
%      confirm the result by calculating the dot product (inner
%      product) of a vector from the origin to the plane and a
%      vector in that plane.
%
%      The dot product between these two vectors should be zero,
%      to double precision round-off, so orthogonal to that
%      precision.
%
%
%      Example code begins here.
%
%
%      function pl2nvc_ex1()
%
%         %
%         % A simple task, determine the distance of a plane
%         % from the origin.
%         %
%         % Define the plane with a vector normal to the plane
%         % and a point in the plane.
%         %
%         normal = [ -1.;  5.;    -3.5 ];
%         point  = [  9.; -0.65;  -12. ];
%
%         %
%         % create the SPICE plane from the normal and point.
%         %
%         plane = cspice_nvp2pl( normal, point );
%
%         %
%         % Calculate the normal vector and constant defining
%         % the plane. The constant value is the distance from
%         % the origin to the plane.
%         %
%         [normal, konst] = cspice_pl2nvc( plane );
%         fprintf( 'Distance to the plane: %12.7f\n', konst );
%
%         %
%         % Confirm the results. Calculate a vector
%         % from the origin to the plane.
%         %
%         vec = konst * normal;
%         fprintf( 'Vector from origin   : %12.7f %12.7f %12.7f\n\n', ...
%                                                               vec );
%
%         %
%         % Now calculate a vector in the plane from the
%         % location in the plane defined by 'vec'.
%         %
%         plane_vec = vec - point;
%
%         %
%         % These vectors should be orthogonal.
%         %
%         fprintf('dot product          : %12.7f\n', ...
%                            dot( plane_vec, vec ) );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Distance to the plane:    4.8102899
%      Vector from origin   :   -0.7777778    3.8888889   -2.7222222
%
%      dot product          :   -0.0000000
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
%      function pl2nvc_ex2()
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
%   makes the routine interfaces simpler and more uniform.
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
%-Exceptions
%
%   1)  The input plane MUST have been created by one of the Mice
%       routines
%
%          cspice_nvc2pl ( Normal vector and constant to plane )
%          cspice_nvp2pl ( Normal vector and point to plane    )
%          cspice_psv2pl ( Point and spanning vectors to plane )
%
%       Otherwise, the results of this routine are unpredictable.
%
%   2)  If the input argument `plane' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `plane' is not of the expected type, or
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
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Changed the argument name "constant" to "konst" for consistency
%       with other routines.
%
%       Edited the -Examples section to comply with NAIF standard.
%       Reformatted example's output, added problem statement and a
%       second example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 27-AUG-2012 (EDW)
%
%-Index_Entries
%
%   plane to normal vector and constant
%
%-&

function [normal, konst] = cspice_pl2nvc( plane )

   switch nargin

      case 1

         plane = zzmice_pln( plane );

      otherwise

         error( ['Usage: [normal(3), konst] = cspice_pl2nvc( plane )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [normal, konst] = mice('pl2nvc_c', plane );
   catch spiceerr
      rethrow(spiceerr)
   end
