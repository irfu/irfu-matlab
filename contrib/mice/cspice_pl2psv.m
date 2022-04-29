%-Abstract
%
%   CSPICE_PL2PSV returns a point and two orthogonal spanning vectors that
%   generate a specified plane.
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
%      plane    a SPICE plane.
%
%               [1,1] = size(plane); struct = class(plane)
%
%               The structure has the fields:
%
%                  normal:   [3,1] = size(normal); double = class(normal)
%                  constant: [1,1] = size(constant); double = class(constant)
%
%   the call:
%
%      [point, span1, span2] = cspice_pl2psv( plane )
%
%   returns:
%
%      point,
%      span1,
%      span2    respectively, a point and two orthogonal
%               spanning vectors that generate the geometric plane
%               represented by `plane'.
%
%               [3,1] = size(point); double = class(point)
%               [3,1] = size(span1); double = class(span1)
%               [3,1] = size(span2); double = class(span2)
%
%               The geometric plane is the set of vectors
%
%                  point   +   s * span1   +   t * span2
%
%               where `s' and `t' are real numbers. `point' is the closest
%               point in the plane to the origin; this point is always a
%               multiple of the plane's normal vector. `span1' and `span2'
%               are an orthonormal pair of vectors. `point', `span1', and
%               `span2' are mutually orthogonal.
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
%   1) Construct a SPICE plane from a normal vector and a point on
%      that plane, and calculate a point and two orthogonal spanning
%      vectors that generate the specified plane.
%
%      Example code begins here.
%
%
%      function pl2psv_ex1()
%         %
%         % Define a normal vector from a plane and a
%         % point in a plane.
%         %
%         normal = [ -1.;  5.;   -3.5 ];
%         point  = [  9.; -0.65; -12. ];
%
%         %
%         % Create a plane from the vectors.
%         %
%         plane = cspice_nvp2pl( normal, point );
%         fprintf( 'Input plane:\n' )
%         fprintf( '  Normal vector  : %15.12f %15.12f %15.12f\n', ...
%                                                        plane.normal   )
%         fprintf( '  Constant       : %15.12f\n\n',      plane.constant )
%
%         %
%         % Calculate a point in the plane, and
%         % two spanning vectors in the plane such that
%         % the point and spanning are mutually orthogonal.
%         %
%         [point, span1, span2] = cspice_pl2psv( plane );
%
%         fprintf( 'Point            : %15.12f %15.12f %15.12f\n',   point )
%         fprintf( 'Spanning vector 1: %15.12f %15.12f %15.12f\n',   span1 )
%         fprintf( 'Spanning vector 2: %15.12f %15.12f %15.12f\n\n', span2 )
%
%         %
%         % Test `point', `span1', and `span2' orthogonality. The dot
%         % products of any two vectors should equal zero to
%         % within round-off.
%         %
%         fprintf( 'dot(point,span1) : %20.17f\n', dot( point, span1) )
%         fprintf( 'dot(point,span2) : %20.17f\n', dot( point, span2) )
%         fprintf( 'dot(span1,span2) : %20.17f\n', dot( span1, span2) )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Input plane:
%        Normal vector  : -0.161690416691  0.808452083454 -0.565916458418
%        Constant       :  4.810289896554
%
%      Point            : -0.777777777778  3.888888888889 -2.722222222222
%      Spanning vector 1:  0.000000000000  0.573462344363  0.819231920519
%      Spanning vector 2:  0.986841531934  0.132461950595 -0.092723365417
%
%      dot(point,span1) :  0.00000000000000000
%      dot(point,span2) :  0.00000000000000006
%      dot(span1,span2) :  0.00000000000000000
%
%
%      Note that, as expected, the dot products of any two vectors equal
%      zero to within round-off.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement and modified code example to produce
%       formatted output.
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
%   plane to point and spanning vectors
%
%-&

function [point, span1, span2] = cspice_pl2psv( plane )

   switch nargin

      case 1

         plane = zzmice_pln( plane );

      otherwise

         error ( ['Usage: [point(3), span1(3), span2(3)] ' ...
                  '= cspice_pl2psv( plane )'] )

   end

   %
   % Call the MEX library.
   %
   % The developer decided to not complicate the interface call and so
   % use the individual fields of the 'plane' structure as arguments.
   %
   try
      [point, span1, span2] = mice('pl2psv_c', plane );
   catch spiceerr
      rethrow(spiceerr)
   end
