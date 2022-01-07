%-Abstract
%
%   CSPICE_PLTVOL computes the volume of a three-dimensional region bounded by
%   a collection of triangular plates.
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
%      vrtces   an array containing the plate model's vertices.
%
%               [3,nv] = size(vrtces); double = class(vrtces)
%
%               Elements
%
%                  vrtces(1,i)
%                  vrtces(2,i)
%                  vrtces(3,i)
%
%               are, respectively, the X, Y, and Z components of
%               the ith vertex, where `i' ranges from 1 to nv.
%
%               This routine doesn't associate units with the
%               vertices.
%
%      plates   an array containing 3-tuples of integers
%               representing the model's plates. The elements of
%               `plates' are vertex indices. The vertex indices are
%               1-based: vertices have indices ranging from 1 to
%               nv.
%
%               [3,np] = size(plates); int32 = class(plates)
%
%               The elements
%
%                  plates(1,i)
%                  plates(2,i)
%                  plates(3,i)
%
%               are, respectively, the indices of the vertices
%               comprising the ith plate.
%
%               Note that the order of the vertices of a plate is
%               significant: the vertices must be ordered in the
%               positive (counterclockwise) sense with respect to
%               the outward normal direction associated with the
%               plate. In other words, if v1, v2, v3 are the
%               vertices of a plate, then
%
%                 ( v2 - v1 )  x  ( v3 - v2 )
%
%               points in the outward normal direction. Here
%               "x" denotes the vector cross product operator.
%
%   the call:
%
%      [pltvol] = cspice_pltvol( vrtces, plates )
%
%   returns:
%
%      pltvol   the volume of the spatial region bounded
%               by the plates.
%
%               [1,1] = size(pltvol); double = class(pltvol)
%
%               If the components of the vertex array have distance unit L,
%               then the output volume has units
%
%                   3
%                  L
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
%   1) Compute the volume of the pyramid defined by the four
%      triangular plates whose vertices are the 3-element
%      subsets of the set of vectors:
%
%         ( 0, 0, 0 )
%         ( 1, 0, 0 )
%         ( 0, 1, 0 )
%         ( 0, 0, 1 )
%
%
%      Example code begins here.
%
%
%      function pltvol_ex1()
%
%         %
%         % Let the notation
%         %
%         %    < A, B >
%         %
%         % denote the dot product of vectors A and B.
%         %
%         % The plates defined below lie in the following planes,
%         % respectively:
%         %
%         %    Plate 1:    { P :  < P, (-1,  0,  0) > = 0 }
%         %    Plate 2:    { P :  < P, ( 0, -1,  0) > = 0 }
%         %    Plate 3:    { P :  < P, ( 0,  0, -1) > = 0 }
%         %    Plate 4:    { P :  < P, ( 1,  1,  1) > = 1 }
%         %
%         vrtces =[  [ 0.0, 0.0, 0.0 ]', ...
%                    [ 1.0, 0.0, 0.0 ]', ...
%                    [ 0.0, 1.0, 0.0 ]', ...
%                    [ 0.0, 0.0, 1.0 ]'  ];
%
%         plates =[ [ 1, 4, 3 ]', ...
%                   [ 1, 2, 4 ]', ...
%                   [ 1, 3, 2 ]', ...
%                   [ 2, 3, 4 ]'  ];
%
%           vol = cspice_pltvol( vrtces, plates );
%
%           fprintf ( 'Expected volume  =      1/6\n'       )
%           fprintf ( 'Computed volume  =   %24.17e\n', vol )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Expected volume  =      1/6
%      Computed volume  =    1.66666666666666657e-01
%
%
%-Particulars
%
%   This routine computes the volume of a spatial region bounded by
%   a set of triangular plates. If the plate set does not actually
%   form the boundary of a spatial region, the result of this routine
%   is invalid.
%
%   Examples:
%
%      Valid inputs
%      ------------
%      Tetrahedron
%      Box
%      Tiled ellipsoid
%      Two disjoint boxes
%
%      Invalid inputs
%      --------------
%      Single plate
%      Tiled ellipsoid with one plate removed
%      Two boxes with intersection having positive volume
%
%-Exceptions
%
%   1)  The input plate model must define a spatial region with
%       a boundary. This routine does not check the inputs to
%       verify this condition. See the -Restrictions section below.
%
%   2)  If the number of vertices is less than 4, the error
%       SPICE(TOOFEWVERTICES) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the number of plates is less than 4, the error
%       SPICE(TOOFEWPLATES) is signaled by a routine in the call tree
%       of this routine.
%
%   4)  If any plate contains a vertex index outside of the range
%
%          [1, nv]
%
%       where `nv' is the number of vertices, the error
%       SPICE(INDEXOUTOFRANGE) is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If any of the input arguments, `vrtces' or `plates', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   6)  If any of the input arguments, `vrtces' or `plates', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The plate collection must describe a surface and enclose a
%       volume such that the divergence theorem (see [1]) is
%       applicable.
%
%-Required_Reading
%
%   DSK.REQ
%   MICE.REQ
%
%-Literature_References
%
%   [1]  T. Apostol, "Calculus, Vol. II," John Wiley & Sons, 1969.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 07-AUG-2020 (EDW) (JDR)
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Fixed
%       minor typos in header.
%
%       Edited the header to comply with NAIF standard.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 16-MAR-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   compute plate model volume
%
%-&

function [pltvol] = cspice_pltvol( vrtces, plates )

   switch nargin
      case 2

         vrtces  = zzmice_dp(vrtces);
         plates  = zzmice_int(plates);

      otherwise

         error ( 'Usage: [pltvol] = cspice_pltvol( vrtces, plates )' )

   end

   %
   % Call the MEX library.
   %
   try
      [pltvol] = mice('pltvol_c', vrtces, plates);
   catch spiceerr
      rethrow(spiceerr)
   end
