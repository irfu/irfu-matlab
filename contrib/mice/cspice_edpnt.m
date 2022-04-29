%-Abstract
%
%   CSPICE_EDPNT scales a point so that it lies on the surface of a specified
%   triaxial ellipsoid that is centered at the origin and aligned
%   with the Cartesian coordinate axes.
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
%      p        a non-zero point in three-dimensional space.
%
%               [3,1] = size(p); double = class(p)
%
%      a,
%      b,
%      c        respectively, the semi-axis lengths of a triaxial
%               ellipsoid in the X, Y, and Z directions.
%
%               [1,1] = size(a); double = class(a)
%               [1,1] = size(b); double = class(b)
%               [1,1] = size(c); double = class(c)
%
%               The axes of the ellipsoid are aligned with the axes of the
%               Cartesian coordinate system.
%
%   the call:
%
%      [ep] = cspice_edpnt( p, a, b, c )
%
%   returns:
%
%      ep       the result of scaling the input point `p' so that it lies on
%               the surface of the triaxial ellipsoid defined by the input
%               semi-axis lengths.
%
%               [3,1] = size(ep); double = class(ep)
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
%   1) Find the surface intercept point on an ellipsoid having radii
%
%          ( 3, 2, 1 )
%
%      of the ray emanating from the origin and having direction
%      vector
%
%          ( 1, 1, 1 )
%
%
%      Example code begins here.
%
%
%      function edpnt_ex1()
%
%         a    = 3.0;
%         b    = 2.0;
%         c    = 1.0;
%
%         v    = [ 1.0, 1.0, 1.0 ]';
%
%         [ep] = cspice_edpnt( v, a, b, c );
%
%         fprintf( 'EP    =  %17.14f %17.14f %17.14f\n', ...
%                                    ep(1), ep(2), ep(3) )
%
%         %
%         % Verify that `ep' is on the ellipsoid.
%         %
%         level =   (ep(1)/a) ^ 2 + (ep(2)/b) ^ 2 + (ep(3)/c) ^ 2;
%
%         fprintf( 'LEVEL =  %17.14f\n', level )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      EP    =   0.85714285714286  0.85714285714286  0.85714285714286
%      LEVEL =   1.00000000000000
%
%
%-Particulars
%
%   This routine efficiently computes the ellipsoid surface point
%   corresponding to a specified ray emanating from the origin.
%   Practical examples of this computation occur in the Mice
%   routines cspice_latsrf and cspice_srfrec.
%
%-Exceptions
%
%   1)  If any of the target ellipsoid's semi-axis lengths is
%       non-positive, the error SPICE(INVALIDAXES) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If `p' is the zero vector, the error SPICE(ZEROVECTOR) is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If the level surface parameter of the input point underflows,
%       the error SPICE(POINTTOOSMALL) is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If any of the input arguments, `p', `a', `b' or `c', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `p', `a', `b' or `c', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   scale point to lie on ellipsoid
%
%-&
function [ep] = cspice_edpnt( p, a, b, c )

   switch nargin
      case 4

         p = zzmice_dp(p);
         a = zzmice_dp(a);
         b = zzmice_dp(b);
         c = zzmice_dp(c);

      otherwise

         error ( 'Usage: [ep(3)] = cspice_edpnt( p(3), a, b, c )' )

   end

   %
   % Call the MEX library.
   %
   try
      [ep] = mice('edpnt_c', p, a, b, c);
   catch spiceerr
      rethrow(spiceerr)
   end
