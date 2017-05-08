%-Abstract
%
%   CSPICE_PLTNRM computes an outward normal vector of a triangular plate.
%   The vector does not necessarily have unit length.
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
%      v1,
%      v2,
%      v3         are 3-vectors constituting the vertices of
%                 a triangular plate.
%
%                 [3,1] = size(v1); double = class(v1)
%                 [3,1] = size(v2); double = class(v2)
%                 [3,1] = size(v3); double = class(v3)
%
%   the call:
%
%      [normal] = cspice_pltnrm(v1, v2, v3,)
%
%   returns:
%
%      normal     is an outward normal vector of the plate defined by
%                 the input vertices. The order of the vertices is
%                 used to determine the choice of normal direction:
%                 the normal vector is
%
%                    ( V2 - V1 ) x ( V3 - V2 )
%
%                 [3,1] = size(normal); double = class(normal)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Compute an outward normal of an equilateral triangle lying
%      in the X-Y plane and centered at the origin.
%
%      function pltnrm_t
%
%         s = sqrt(3.0)/2;
%
%         v1 = [  s,  -0.5,  0.0]';
%         v2 = [ 0.0,  1.0,  0.0]';
%         v3 = [ -s,  -0.5,  0.0]';
%
%
%         normal = cspice_pltnrm( v1, v2, v3 );
%
%         fprintf ( 'NORMAL = %18.11e %18.11e %18.11e\n', ...
%                   normal(1), normal(2), normal(3)      );
%
%   Matlab outputs:
%
%         NORMAL =  0.00000000000e+00  0.00000000000e+00  2.59807621135e+00
%
%-Particulars
%
%   This routine saves computation time by not scaling the output
%   vector to unit length. The caller can scale the vector using
%   the routine cspice_vhat.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine pltnrm_c.
%
%   MICE.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 17-MAR-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   compute normal vector of triangular plate from vertices
%
%-&

function [normal] = cspice_pltnrm( v1, v2, v3)

   switch nargin
      case 3

         v1    = zzmice_dp(v1);
         v2    = zzmice_dp(v2);
         v3    = zzmice_dp(v3);

      otherwise

         error ( 'Usage: [normal] = cspice_pltnrm( v1, v2, v3)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [normal] = mice( 'pltnrm_c', v1, v2, v3);
   catch
      rethrow(lasterror)
   end



