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
%               plate.
%
%               [3,3]   = size(iverts); double = class(iverts)
%
%               Each vertex is a three-dimensional vector. The elements
%
%                 iverts(j,i), j = 1 ... 3
%
%               are, respectively, the X, Y, and Z components of the
%               ith vertex.
%
%
%      delta    a fraction by which the plate is to be scaled.
%
%               [1,1]   = size(delta); double = class(delta)
%
%               Scaling is done so that the scaled plate has the following
%               properties:
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
%   the call:
%
%      [overts] = cspice_pltexp( iverts, delta )
%
%   returns:
%
%      overts   an array containing three vertices of the triangular
%               plate resulting from scaling the input plate.
%
%               [3,3]   = size(overts); double = class(overts)
%
%               If `ctroid' is the centroid (the average of the vertices)
%               of the input plate, then the jth vertex of `overts'
%
%                  overts(j,i), j = 1 ... 3
%
%               is equal to
%
%                  ctroid(j) + (1+delta)*( iverts(j,i) - ctroid(j) ),
%
%                  j = 1 ... 3
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
%   1) Expand an equilateral triangle that lies in the plane
%
%         { (x,y,z) : z = 7 }
%
%      Use an expansion fraction of 1.0; this doubles the size of
%      the plate.
%
%      Example code begins here.
%
%
%      function pltexp_ex1()
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
%         fprintf ( '\nVertices of output plate: \n')
%         fprintf ( ' O1 = %20.12f %20.12f %20.12f\n', overts(:,1) )
%         fprintf ( ' O2 = %20.12f %20.12f %20.12f\n', overts(:,2) )
%         fprintf ( ' O3 = %20.12f %20.12f %20.12f\n', overts(:,3) )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
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
%
%-Particulars
%
%   This routine supports "greedy" ray-plate intercept algorithms.
%   Such algorithms attempt to ensure that false negatives---in which
%   an intersection is not found due to round-off error---do not
%   occur. In such an algorithm, the plate of interest is expanded
%   slightly before the intersection test is performed.
%
%-Exceptions
%
%   1)  If any of the input arguments, `iverts' or `delta', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   2)  If any of the input arguments, `iverts' or `delta', is not of
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
%   DSK.REQ
%   MICE.REQ
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
%   -Mice Version 1.1.0, 07-AUG-2020 (EDW) (JDR)
%
%       Updated description of input argument "iverts".
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Fixed
%       minor typos in header.
%
%       Edited the header to comply with NAIF standard. Added
%       example task statement.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 28-NOV-2016 (EDW)
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
   catch spiceerr
      rethrow(spiceerr)
   end
