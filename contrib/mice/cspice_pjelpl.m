%-Abstract
%
%   CSPICE_PJELPL projects orthogonally an ellipse onto a plane.
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
%      elin    a structure describing a SPICE ellipse.
%
%              [1,1] = size(elin); struct = class(elin)
%
%              The structure has the fields:
%
%                 center:    [3x1 double]
%                 semiMajor: [3x1 double]
%                 semiMinor: [3x1 double]
%
%      plane   a structure describing a SPICE plane.
%
%              [1,1] = size(plane); struct = class(plane)
%
%              The structure has the fields:
%
%                  normal:     [3x1 double]
%                  constant:   [1x1 double]
%
%              are, respectively, a SPICE ellipse and a SPICE plane. The
%              geometric ellipse represented by 'elin' is to be orthogonally
%              projected onto the geometric plane represented by 'plane'.
%
%   the call:
%
%      elout = cspice_pjelpl( elin, plane )
%
%   returns:
%
%      elout   the SPICE ellipse that represents the geometric
%              ellipse resulting from orthogonally projecting the ellipse
%              represented by 'elin' onto the plane represented by 'plane'.
%
%              [1,1] = size(elout); struct = class(elout)
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
%   1) Given an ellipse and a plane, compute the projection of the
%      ellipse on the plane.
%
%      Example code begins here.
%
%
%      function pjelpl_ex1()
%
%         %
%         % Assign the values for plane/ellipse definition
%         % vectors.
%         %
%         center  = [ 1,  1,  1 ]';
%         vect1   = [ 2,  0,  0 ]';
%         vect2   = [ 0,  1,  1 ]';
%         normal  = [ 0,  0,  1 ]';
%
%         %
%         % Create a plane using a constant value of 0...
%         %
%         plane = cspice_nvc2pl( normal, 0 );
%
%         %
%         % ...and an ellipse.
%         %
%         elin = cspice_cgv2el( center, vect1, vect2 );
%
%         %
%         % Project the ellipse onto the plane.
%         %
%         elout = cspice_pjelpl( elin, plane );
%
%         %
%         % Output the ellipse in the plane.
%         %
%         fprintf( 'Center    :  %f  %f  %f\n', elout.center    )
%         fprintf( 'Semi-minor:  %f  %f  %f\n', elout.semiMinor )
%         fprintf( 'Semi-major:  %f  %f  %f\n', elout.semiMajor )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Center    :  1.000000  1.000000  0.000000
%      Semi-minor:  0.000000  1.000000  0.000000
%      Semi-major:  2.000000  0.000000  0.000000
%
%
%-Particulars
%
%   Projecting an ellipse orthogonally onto a plane can be thought of
%   finding the points on the plane that are `under' or `over' the
%   ellipse, with the `up' direction considered to be perpendicular
%   to the plane. More mathematically, the orthogonal projection is
%   the set of points Y in the plane such that for some point X in
%   the ellipse, the vector Y - X is perpendicular to the plane.
%   The orthogonal projection of an ellipse onto a plane yields
%   another ellipse.
%
%-Exceptions
%
%   1)  If the input plane is invalid, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  The input ellipse may be degenerate--its semi-axes may be
%       linearly dependent. Such ellipses are allowed as inputs.
%
%   3)  The ellipse resulting from orthogonally projecting the input
%       ellipse onto a plane may be degenerate, even if the input
%       ellipse is not.
%
%   4)  If any of the input arguments, `elin' or `plane', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `elin' or `plane', is not of
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
%   ELLIPSES.REQ
%   PLANES.REQ
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
%   -Mice Version 1.1.0, 27-AUG-2021 (EDW) (JDR)
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 11-JUN-2013 (EDW)
%
%-Index_Entries
%
%   project ellipse onto plane
%
%-&

function [elout] = cspice_pjelpl( elin, plane )

   switch nargin
      case 2

         elin  = zzmice_ell(elin);
         plane = zzmice_pln(plane);

      otherwise

         error ( 'Usage: [elout] = cspice_pjelpl( elin, plane )' )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [elout] = mice( 'pjelpl_s', elin, plane );
   catch spiceerr
      rethrow(spiceerr)
   end
