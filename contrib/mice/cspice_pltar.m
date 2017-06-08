%-Abstract
%
%   CSPICE_PLTAR computes the total area of a collection of triangular plates.
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
%               [3,m] = size(vrtces); double = class(vrtces)
%
%               Elements
%
%                  vrtces(1,i)
%                  vrtces(2,i)
%                  vrtces(3,i)
%
%               are, respectively, the X, Y, and Z components of
%               the ith vertex, where `i' ranges from 1 to m.
%
%               This routine doesn't associate units with the
%               vertices.
%
%      plates   is an array containing 3-tuples of integers
%               representing the model's plates. The elements of
%               `plates' are vertex indices. The vertex indices are
%               1-based: vertices have indices ranging from 1 to
%               n.
%
%               [3,n] = size(plates); int32 = class(plates)
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
%      pltar = cspice_pltar( vrtces, plates )
%
%   returns:
%
%      pltar   The function returns the total area of the input set of
%              plates. Each plate contributes the area of the triangle
%              defined by the plate's vertices.
%
%              [1,1] = size(pltar); double = class(pltar)
%
%              If the components of the vertex array have length unit L, then
%              the output area has units
%
%               2
%              L
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Compute the area of the pyramid defined by the four
%      triangular plates whose vertices are the 3-element
%      subsets of the set of vectors:
%
%         ( 0, 0, 0 )
%         ( 1, 0, 0 )
%         ( 0, 1, 0 )
%         ( 0, 0, 1 )
%
%      function pltar_t
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
%           area = cspice_pltar( vrtces, plates );
%
%           fprintf ( ['Expected area   =    (3 + sqrt(3))/2\n' ...
%                      '                =    0.23660254037844384e+01\n'] )
%           fprintf (  'Computed volume =   %24.17e\n', area )
%
%   Matlab outputs:
%
%      Expected area   =    (3 + sqrt(3))/2
%                      =    0.23660254037844384e+01
%      Computed volume =    2.36602540378443837e+00
%
%-Particulars
%
%   This routine computes the total area of a set of triangular
%   plates. The plates need not define a closed surface.
%
%   Examples of valid plate sets:
%
%      Tetrahedron
%      Box
%      Tiled ellipsoid
%      Tiled ellipsoid with one plate removed
%      Two disjoint boxes
%      Two boxes with intersection having positive volume
%      Single plate
%      Empty plate set
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pltar_c.
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
%   compute plate model area
%
%-&

function [pltar] = cspice_pltar( vrtces, plates )

   switch nargin
      case 2

         vrtces  = zzmice_dp(vrtces);
         plates  = zzmice_int(plates);

      otherwise

         error ( 'Usage: [pltar] = cspice_pltar( vrtces, plates )' )

   end

   %
   % Call the MEX library.
   %
   try
      [pltar] = mice('pltar_c', vrtces, plates);
   catch
      rethrow(lasterror)
   end
