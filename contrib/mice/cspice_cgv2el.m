%-Abstract
%
%   CSPICE_CGV2EL forms a SPICE ellipse from a center vector and two
%   generating vectors.
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
%      center   the location for an ellipse center.
%
%               [3,1] = size(center); double = class(center)
%
%      vec1,
%      vec2     the two vectors defining the ellipse (the generating vectors)
%               with the `center' in three-dimensional space. The ellipse is
%               the set of points
%
%                  center  +  cos(theta) vec1  +  sin(theta) vec2
%
%               where theta ranges over the interval (-pi, pi].
%
%               `vec1' and `vec2' need not be linearly independent.
%
%               [3,1] = size(vec1); double = class(vec1)
%               [3,1] = size(vec2); double = class(vec2)
%
%   the call:
%
%      [ellips] = cspice_cgv2el( center, vec1, vec2 )
%
%   returns:
%
%      ellips   a structure describing a SPICE ellipse defined by the input
%               vectors.
%
%               [1,1] = size(ellips); struct = class(ellips)
%
%               The structure has the fields:
%
%               center:    [3,1] = size(center); double = class(center)
%               semiMinor: [3,1] = size(semiMinor); double = class(semiMinor)
%               semiMajor: [3,1] = size(semiMajor); double = class(semiMajor)
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
%   1) Create a SPICE ellipse given its center and two linearly independent
%      generating vectors of the ellipse.
%
%      Example code begins here.
%
%
%      function cgv2el_ex1()
%
%         %
%         % Define the center and two linearly independent
%         % generating vectors of an ellipse (the vectors need not
%         % be linearly independent).
%         %
%         center = [ -1.;  1.; -1. ];
%         vec1   = [  1.;  1.;  1. ];
%         vec2   = [  1.; -1.;  1. ];
%
%         %
%         % Create the CSPICE_ELLIPSE structure.
%         %
%         ellips = cspice_cgv2el( center, vec1, vec2 );
%
%         fprintf( 'SPICE ellipse:\n' );
%         fprintf( '  Semi-minor axis: %10.6f %10.6f %10.6f\n',            ...
%                                              ellips.semiMinor);
%         fprintf( '  Semi-major axis: %10.6f %10.6f %10.6f\n',            ...
%                                              ellips.semiMajor);
%         fprintf( '  Center         : %10.6f %10.6f %10.6f\n',            ...
%                                              ellips.center   );
%         fprintf( '\n' );
%
%         %
%         % Obtain the center and generating vectors from the
%         % `ellips'.
%         %
%         [ecentr, smajor, sminor] = cspice_el2cgv( ellips );
%         fprintf( 'SPICE ellipse (using cspice_el2cgv):\n' );
%         fprintf( '  Semi-minor axis: %10.6f %10.6f %10.6f\n', sminor);
%         fprintf( '  Semi-major axis: %10.6f %10.6f %10.6f\n', smajor);
%         fprintf( '  Center         : %10.6f %10.6f %10.6f\n', ecentr);
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      SPICE ellipse:
%        Semi-minor axis:   0.000000   1.414214   0.000000
%        Semi-major axis:   1.414214  -0.000000   1.414214
%        Center         :  -1.000000   1.000000  -1.000000
%
%      SPICE ellipse (using cspice_el2cgv):
%        Semi-minor axis:   0.000000   1.414214   0.000000
%        Semi-major axis:   1.414214  -0.000000   1.414214
%        Center         :  -1.000000   1.000000  -1.000000
%
%
%   2) Find the intersection of an ellipse with a plane.
%
%
%      Example code begins here.
%
%
%      function cgv2el_ex2()
%
%         %
%         % Local variables.
%         %
%         xpts = zeros(2,3);
%
%         %
%         % The ellipse is defined by the vectors `center', `vec1', and
%         % `vec2'. The plane is defined by the normal vector `normal'
%         % and the `center'.
%         %
%         center = [ 0.0,  0.0,  0.0]';
%         vec1   = [ 1.0,  7.0,  2.0]';
%         vec2   = [-1.0,  1.0,  3.0]';
%
%         normal = [ 0.0,  1.0,  0.0]';
%
%         %
%         % Make a SPICE ellipse and a plane.
%         %
%         [ellips] = cspice_cgv2el( center, vec1, vec2 );
%         [plane]  = cspice_nvp2pl( normal, center     );
%
%         %
%         % Find the intersection of the ellipse and plane.
%         % `nxpts' is the number of intersection points; `xpts'
%         % are the points themselves.
%         %
%         [nxpts, xpts(1,:), xpts(2,:)] = cspice_inelpl( ellips, plane );
%
%         fprintf( 'Number of intercept points: %2d\n', nxpts )
%
%         for i=1:nxpts
%            fprintf( ' Point %1d : %9.6f %9.6f %9.6f\n', i, xpts(i,:) )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Number of intercept points:  2
%       Point 1 :  1.131371  0.000000 -2.687006
%       Point 2 : -1.131371 -0.000000  2.687006
%
%
%-Particulars
%
%   SPICE ellipses serve to simplify calling sequences and reduce
%   the chance for error in declaring and describing argument lists
%   involving ellipses.
%
%   The set of ellipse conversion routines is
%
%      cspice_cgv2el( Center and generating vectors to ellipse )
%      cspice_el2cgv( Ellipse to center and generating vectors )
%
%-Exceptions
%
%   1)  If `vec1' and `vec2' are linearly dependent, `ellips' will be
%       degenerate. SPICE ellipses are allowed to represent
%       degenerate geometric ellipses.
%
%   2)  If any of the input arguments, `center', `vec1' or `vec2', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `center', `vec1' or `vec2', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
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
%   ELLIPSES.REQ
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
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "ellipse" to "ellips".
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement, reformatted example's output and added
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
%   -Mice Version 1.0.1, 09-NOV-2012 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 30-DEC-2008 (EDW)
%
%-Index_Entries
%
%   center and generating vectors to ellipse
%
%-&
function [ellips] = cspice_cgv2el( center, vec1, vec2 )

   switch nargin

      case 3

         center = zzmice_dp(center);
         vec1   = zzmice_dp(vec1);
         vec2   = zzmice_dp(vec2);

      otherwise

         error ( ['Usage: [ellips] = ' ...
                  'cspice_cgv2el( center(3), vec1(3), vec2(3) )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [ellips] = mice('cgv2el_c', center, vec1, vec2 );
   catch spiceerr
      rethrow(spiceerr)
   end
