%-Abstract
%
%   CSPICE_PLTEXP expands a triangular plate by a specified amount.
%   The expanded plate is co-planar with, and has the same orientation as,
%   the  original. The centroids of the two plates coincide.
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
%      iverts   an array containing three vertices of a triangular
%               plate. Each vertex is a three-dimensional vector. The
%               elements
%
%               [3,3]   = size(iverts); double = class(iverts)
%
%      delta    a fraction by which the plate is to be scaled.
%               Scaling is done so that the scaled plate has the
%               following properties:
%
%                  -  it is co-planar with the input plate
%
%                  -  its centroid coincides with that of the input
%                     plate
%
%                  -  its sides remain parallel to the corresponding
%                     sides of the input plate
%
%                  -  the distance of each vertex from the centroid is
%                     (1+delta) times the corresponding distance for
%                     the input plate
%
%                 [1,1]   = size(delta); double = class(delta)
%
%   the call:
%
%      overts = cspice_pltexp( iverts, delta)
%
%   returns:
%
%      overts   an array containing three vertices of the triangular
%               plate resulting from scaling the input plate.
%
%               If `ctroid' is the centroid (the average of the vertices)
%               of the input plate, then the jth vertex of `overts'
%
%                  overts(i,j), i = 1 ... 3
%
%               is equal to
%
%                  ctroid(i) + (1+delta)*( iverts(i,j) - ctroid(i) ),
%
%                  i = 1 ... 3
%
%
%               [3,3]   = size(overts); double = class(overts)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      function pltexp_t
%
%         s     = sqrt( 3.0 ) / 2.0;
%
%         iverts = [ [ s; -0.5; 7.0] [ 0.0; 1.0; 7.0] [ -s; -0.5; 7.0] ];
%
%
%         delta = 1.0;
%
%         overts = cspice_pltexp ( iverts, delta );
%
%         fprintf ( '\nVertices of input plate: \n')
%         fprintf ( ' I1 = %20.12f %20.12f %20.12f\n', iverts(:,1) )
%         fprintf ( ' I2 = %20.12f %20.12f %20.12f\n', iverts(:,2) )
%         fprintf ( ' I3 = %20.12f %20.12f %20.12f\n', iverts(:,3) )
%
%
%
%         fprintf ( '\nVertices of output plate: \n')
%         fprintf ( ' O1 = %20.12f %20.12f %20.12f\n', overts(:,1) )
%         fprintf ( ' O2 = %20.12f %20.12f %20.12f\n', overts(:,2) )
%         fprintf ( ' O3 = %20.12f %20.12f %20.12f\n', overts(:,3) )
%
%   MATLAB outputs:
%
%      Vertices of input plate:
%       I1 =       0.866025403784      -0.500000000000       7.000000000000
%       I2 =       0.000000000000       1.000000000000       7.000000000000
%       I3 =      -0.866025403784      -0.500000000000       7.000000000000
%
%      Vertices of output plate:
%       O1 =       1.732050807569      -1.000000000000       7.000000000000
%       O2 =       0.000000000000       2.000000000000       7.000000000000
%       O3 =      -1.732050807569      -1.000000000000       7.000000000000
%
%-Particulars
%
%   This routine supports "greedy" ray-plate intercept algorithms.
%   Such algorithms attempt to ensure that false negatives---in which
%   an intersection is not found due to round-off error---do not
%   occur. In such an algorithm, the plate of interest is expanded
%   slightly before the intersection test is performed.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pltexp_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 28-NOV-2016, EDW (JPL)
%
%-Index_Entries
%
%   expand triangular plate
%
%-&

function [overts] = cspice_pltexp( iverts, delta)

   switch nargin
      case 2

         iverts = zzmice_dp(iverts);
         delta  = zzmice_dp(delta);

      otherwise

         error ( [ 'Usage: [overts(3,3)] = ' ...
                   'cspice_pltexp( iverts(3,3), delta)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [overts] = mice('pltexp_c', iverts, delta );
   catch
      rethrow(lasterror)
   end


