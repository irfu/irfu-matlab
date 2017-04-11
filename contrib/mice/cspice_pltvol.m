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
%      vrtces   is an array containing the plate model's vertices.
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
%      plates   is an array containing 3-tuples of integers
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
%      pltvol = cspice_pltvol( vrtces, plates )
%
%   returns:
%
%      pltvol   The function returns the volume of the spatial region bounded
%               by the plates.
%
%               [1,1] = size(pltvol); double = class(pltvol)
%
%               If the components of the vertex array have distance unit L,
%               then the output volume has units
%
%                3
%               L
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Compute the volume of the pyramid defined by the four
%      triangular plates whose vertices are the 3-element
%      subsets of the set of vectors:
%
%         ( 0, 0, 0 )
%         ( 1, 0, 0 )
%         ( 0, 1, 0 )
%         ( 0, 0, 1 )
%
%      function pltvol_t
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
%           fprintf ( 'Expected volume =      1/6\n'        )
%           fprintf ( 'Computed volume  =   %24.17e\n', vol )
%
%   Matlab outputs:
%
%      Expected volume =      1/6
%      Computed volume  =    1.66666666666666657e-01
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
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pltvol_c.
%
%   MICE.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 16-MAR-2016, EDW (JPL), NJB (JPL)
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
   catch
      rethrow(lasterror)
   end
