%-Abstract
%
%   CSPICE_SURFPV finds the state (position and velocity) of the surface
%   intercept defined by a specified ray, ray velocity, and ellipsoid.
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
%      stvrtx   the state of a ray's vertex.
%
%               [6,1] = size(stvrtx); double = class(stvrtx)
%
%               The first three components of `stvrtx' are the vertex's x,
%               y, and z position components; the vertex's x, y, and z
%               velocity components follow.
%
%               The reference frame relative to which `stvrtx' is
%               specified has axes aligned with with those of a
%               triaxial ellipsoid. See the description below of
%               the arguments `a', `b', and `c'.
%
%               The vertex may be inside or outside of this
%               ellipsoid, but not on it, since the surface
%               intercept is a discontinuous function at
%               vertices on the ellipsoid's surface.
%
%               No assumption is made about the units of length
%               and time, but these units must be consistent with
%               those of the other inputs.
%
%
%      stdir    the state of the input ray's direction vector.
%
%               [6,1] = size(stdir); double = class(stdir)
%
%               The first three components of `stdir' are a non-zero vector
%               giving the x, y, and z components of the ray's direction; the
%               direction vector's x, y, and z velocity components follow.
%
%               `stdir' is specified relative to the same reference
%               frame as is `stvrtx'.
%
%      a,
%      b,
%      c        respectively, the lengths of a triaxial ellipsoid's semi-axes
%               lying along the x, y, and z axes of the reference frame
%               relative to which `stvrtx' and `stdir' are specified.
%
%               [1,1] = size(a); double = class(a)
%               [1,1] = size(b); double = class(b)
%               [1,1] = size(c); double = class(c)
%
%   the call:
%
%      [stx, found] = cspice_surfpv( stvrtx, stdir, a, b, c )
%
%   returns:
%
%      stx      the state of the intercept of the input ray on the surface of
%               the input ellipsoid.
%
%               [6,1] = size(stx); double = class(stx)
%
%               The first three components of `stx' are the intercept's x,
%               y, and z position components; the intercept's x, y, and z
%               velocity components follow.
%
%               `stx' is specified relative to the same reference
%               frame as are `stvrtx' and `stdir'.
%
%               `stx' is defined if and only if both the intercept
%               and its velocity are computable, as indicated by
%               the output argument `found'.
%
%               The position units of `stx' are the same as those of
%               `stvrtx', `stdir', and `a', `b', and `c'. The time units are
%               the same as those of `stvrtx' and `stdir'.
%
%
%      found    a logical flag indicating whether `stx' is defined.
%
%               [1,1] = size(found); logical = class(found)
%
%               `found' is true if and only if both the intercept and its
%               velocity are computable. Note that in some cases the
%               intercept may computable while the velocity is not; this can
%               happen for near-tangency cases.
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
%   1) Illustrate the role of the ray vertex velocity and
%      ray direction vector velocity via several simple cases. Also
%      show the results of a near-tangency computation.
%
%
%      Example code begins here.
%
%
%      function surfpv_ex1()
%
%         a        = 1.0;
%         b        = 2.0;
%         c        = 3.0;
%
%         fprintf( '\n' )
%         fprintf( 'Ellipsoid radii:\n' )
%         fprintf( '     A = %f\n', a )
%         fprintf( '     B = %f\n', b )
%         fprintf( '     C = %f\n', c )
%
%         fprintf( '\n' )
%         fprintf( 'Case 1: Vertex varies, direction is constant\n' )
%         fprintf( '\n' )
%
%         stvrtx   = [ 2.0, 0.0, 0.0, 0.0, 0.0, 3.0 ]';
%
%         stdir    = [ -1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]';
%
%         fprintf( 'Vertex:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stvrtx(1:3) )
%         fprintf( 'Vertex velocity:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stvrtx(4:6) )
%         fprintf( 'Direction:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stdir(1:3) )
%         fprintf( 'Direction velocity:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stdir(4:6) )
%
%         [stx, found] = cspice_surfpv( stvrtx, stdir, a, b, c );
%
%         if ( ~ found )
%
%            fprintf( ' No intercept state found.\n' )
%
%         else
%
%            fprintf( 'Intercept:\n' )
%            fprintf( '  %19.12e %19.12e %19.12e\n', stx(1:3) )
%            fprintf( 'Intercept velocity:\n' )
%            fprintf( '  %19.12e %19.12e %19.12e\n', stx(4:6) )
%            fprintf( '\n' )
%
%         end
%
%         fprintf( '\n' )
%         fprintf( 'Case 2: Vertex and direction both vary\n' )
%         fprintf( '\n' )
%
%         stdir(6) =  4.0;
%
%         fprintf( 'Vertex:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stvrtx(1:3) )
%         fprintf( 'Vertex velocity:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stvrtx(4:6) )
%         fprintf( 'Direction:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stdir(1:3) )
%         fprintf( 'Direction velocity:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stdir(4:6) )
%
%         [stx, found] = cspice_surfpv( stvrtx, stdir, a, b, c );
%
%         if ( ~ found )
%
%            fprintf( ' No intercept state found.\n' )
%
%         else
%
%            fprintf( 'Intercept:\n' )
%            fprintf( '  %19.12e %19.12e %19.12e\n', stx(1:3) )
%            fprintf( 'Intercept velocity:\n' )
%            fprintf( '  %19.12e %19.12e %19.12e\n', stx(4:6) )
%            fprintf( '\n' )
%
%         end
%
%         fprintf( '\n' )
%         fprintf( 'Case 3: Vertex and direction both vary;\n' )
%         fprintf( '        near-tangent case.\n' )
%         fprintf( '\n' )
%
%         stvrtx(3) = c - 1.e-15;
%         stvrtx(6) =  1.e299;
%         stdir(6)  =  1.e299;
%
%         fprintf( 'Vertex:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stvrtx(1:3) )
%         fprintf( 'Vertex velocity:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stvrtx(4:6) )
%         fprintf( 'Direction:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stdir(1:3) )
%         fprintf( 'Direction velocity:\n' )
%         fprintf( '  %19.12e %19.12e %19.12e\n', stdir(4:6) )
%
%         [stx, found] = cspice_surfpv( stvrtx, stdir, a, b, c );
%
%         if ( ~ found )
%
%            fprintf( ' No intercept state found.\n' )
%
%         else
%
%            fprintf( 'Intercept:\n' )
%            fprintf( '  %19.12e %19.12e %19.12e\n', stx(1:3) )
%            fprintf( 'Intercept velocity:\n' )
%            fprintf( '  %19.12e %19.12e %19.12e\n', stx(4:6) )
%            fprintf( '\n' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Ellipsoid radii:
%           A = 1.000000
%           B = 2.000000
%           C = 3.000000
%
%      Case 1: Vertex varies, direction is constant
%
%      Vertex:
%         2.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Vertex velocity:
%         0.000000000000e+00  0.000000000000e+00  3.000000000000e+00
%      Direction:
%        -1.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Direction velocity:
%         0.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Intercept:
%         1.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Intercept velocity:
%         0.000000000000e+00  0.000000000000e+00  3.000000000000e+00
%
%
%      Case 2: Vertex and direction both vary
%
%      Vertex:
%         2.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Vertex velocity:
%         0.000000000000e+00  0.000000000000e+00  3.000000000000e+00
%      Direction:
%        -1.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Direction velocity:
%         0.000000000000e+00  0.000000000000e+00  4.000000000000e+00
%      Intercept:
%         1.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Intercept velocity:
%         0.000000000000e+00  0.000000000000e+00  7.000000000000e+00
%
%
%      Case 3: Vertex and direction both vary;
%              near-tangent case.
%
%      Vertex:
%         2.000000000000e+00  0.000000000000e+00  3.000000000000e+00
%      Vertex velocity:
%         0.000000000000e+00  0.000000000000e+00 1.000000000000e+299
%      Direction:
%        -1.000000000000e+00  0.000000000000e+00  0.000000000000e+00
%      Direction velocity:
%         0.000000000000e+00  0.000000000000e+00 1.000000000000e+299
%      Intercept:
%         2.580956827952e-08  0.000000000000e+00  3.000000000000e+00
%      Intercept velocity:
%        -3.874532036208e+306  0.000000000000e+00 2.999999974190e+299
%
%
%-Particulars
%
%   The position and velocity of the ray's vertex as well as the
%   ray's direction vector and velocity vary with time. The
%   inputs to cspice_surfpv may be considered the values of these
%   vector functions at a particular time, say t0. Thus
%
%      State of vertex:            stvrtx = ( v(t0), v'(t0) )
%
%      State of direction vector: stdir  = ( d(t0), d'(t0) )
%
%   To determine the intercept point, w(t0), we simply compute the
%   intersection of the ray originating at v(t0) in the direction of
%   d(t0) with the ellipsoid
%
%         2        2        2
%        x        y        z
%      ----- +  ----- +  -----  =  1
%         2        2        2
%        a        b        c
%
%   w(t) is the path of the intercept point along the surface of
%   the ellipsoid. To determine the velocity of the intercept point,
%   we need to take the time derivative of w(t), and evaluate it at
%   t0. Unfortunately w(t) is a complicated expression, and its
%   derivative is even more complicated.
%
%   However, we know that the derivative of w(t) at t0, w'(t0), is
%   tangent to w(t) at t0. Thus w'(t0) lies in the plane that is
%   tangent to the ellipsoid at t0. Let x(t) be the curve in the
%   tangent plane that represents the intersection of the ray
%   emanating from v(t0) with direction d(t0) with that tangent
%   plane.
%
%      x'(t0) = w'(t0)
%
%   The expression for x'(t) is much simpler than that of w'(t);
%   cspice_surfpv evaluates x'(t) at t0.
%
%
%   Derivation of x(t) and x'(t)
%   ----------------------------
%
%   w(t0) is the intercept point. Let `n' be a surface normal at i(t0).
%   Then the tangent plane at w(t0) is the set of points x(t) such
%   that
%
%      < x(t) - i(t0), n > = 0
%
%   x(t) can be expressed as the vector sum of the vertex
%   and some scalar multiple of the direction vector,
%
%      x(t) = v(t) + s(t) * d(t)
%
%   where s(t) is a scalar function of time. The derivative of
%   x(t) is given by
%
%      x'(t) = v'(t)  +  s(t) * d'(t)  +  s'(t) * d(t)
%
%   We have v(t0), V'(t0), d(t0), D'(t0), w(t0), and `n', but to
%   evaluate X'(t0), we need s(t0) and s'(t0). We derive an
%   expression for s(t) as follows.
%
%   Because x(t) is in the tangent plane, it must satisfy
%
%      < x(t) - w(t0), n > = 0.
%
%   Substituting the expression for x(t) into the equation above
%   gives
%
%      < v(t) + s(t) * d(t) - w(t0), n > = 0.
%
%   Thus
%
%      < v(t) - w(t0), n >  +  s(t) * < d(t), n > = 0,
%
%   and
%                  < v(t) - w(t0), n >
%      s(t)  =  -  -------------------
%                      < d(t), n >
%
%   The derivative of s(t) is given by
%
%      s'(t) =
%
%          < d(t),n > * < v'(t),n >  -  < v(t)-w(t0),n > * < d'(t),n >
%      -   -----------------------------------------------------------
%                                           2
%                                < d(t), n >
%
%-Exceptions
%
%   1)  If the input ray's direction vector is the zero vector, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If any of the ellipsoid's axis lengths is nonpositive, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   3)  If the vertex of the ray is on the ellipsoid, the error
%       SPICE(INVALIDVERTEX) is signaled by a routine in the call tree
%       of this routine.
%
%   4)  If any of the input arguments, `stvrtx', `stdir', `a', `b' or
%       `c', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   5)  If any of the input arguments, `stvrtx', `stdir', `a', `b' or
%       `c', is not of the expected type, or it does not have the
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
%   ellipsoid surface point and velocity
%
%-&
function [stx, found] = cspice_surfpv( stvrtx, stdir, a, b, c )

   switch nargin
      case 5

         stvrtx = zzmice_dp(stvrtx);
         stdir  = zzmice_dp(stdir);
         a      = zzmice_dp(a);
         b      = zzmice_dp(b);
         c      = zzmice_dp(c);

      otherwise

         error ( [ 'Usage: [stx(6), found] = '                              ...
                   'cspice_surfpv( stvrtx(6), stdir(6), a, b, c )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [stx, found] = mice('surfpv_c', stvrtx, stdir, a, b, c);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
