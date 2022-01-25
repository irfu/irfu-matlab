%-Abstract
%
%   CSPICE_UNORM normalizes a double precision 3-vector and returns its
%   magnitude.
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
%      v1       arbitrary 3-vector(s), including the zero vector.
%
%               [3,n] = size(v1); double = class(v1)
%
%   the call:
%
%      [vout, vmag] = cspice_unorm( v1 )
%
%   returns:
%
%      vout     the unit vector(s) in the direction of `v1'.
%
%               [3,n] = size(vout); double = class(vout)
%
%               If `v1' is the zero vector, then `vout' will also be the
%               zero vector.
%
%      vmag     the magnitude(s) of `v1'.
%
%               [1,n] = size(vmag); double = class(vmag)
%
%               `vout' and `vmag' return with the same vectorization measure,
%               N, as `v1'.
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
%   1) Define a set of vectors and compute their corresponding
%      unit vectors and magnitudes.
%
%
%      Example code begins here.
%
%
%      function unorm_ex1()
%
%         %
%         % Local parameters.
%         %
%         SETSIZ = 2;
%
%         %
%         % Define the vector set.
%         %
%         seta = [ [5.0,  12.0,  0.0]', [1.e-7,  2.e-7, 2.e-7]' ];
%
%         %
%         % Calculate the unit vectors and magnitudes.
%         %
%         fprintf('Scalar case:\n')
%
%         for i=1:SETSIZ
%
%            [vout, vmag] = cspice_unorm( seta(:,i) );
%
%            fprintf( 'Vector     :  %12.8f %12.8f %12.8f\n',              ...
%                             seta(1,i), seta(2,i), seta(3,i) )
%            fprintf( 'Unit vector:  %12.8f %12.8f %12.8f\n',              ...
%                                   vout(1), vout(2), vout(3) )
%            fprintf( 'Magnitude  :  %12.8f\n', vmag )
%            fprintf( ' \n' )
%
%         end
%
%         %
%         % Repeat the operation with one single call to cspice_vhat.
%         %
%         [vout, vmag] = cspice_unorm( seta );
%
%         fprintf('Vectorized case:\n')
%
%         for i=1:SETSIZ
%
%            fprintf( 'Vector     :  %12.8f %12.8f %12.8f\n',              ...
%                             seta(1,i), seta(2,i), seta(3,i) )
%            fprintf( 'Unit vector:  %12.8f %12.8f %12.8f\n',              ...
%                             vout(1,i), vout(2,i), vout(3,i) )
%            fprintf( 'Magnitude  :  %12.8f\n', vmag(i) )
%            fprintf( ' \n' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Scalar case:
%      Vector     :    5.00000000  12.00000000   0.00000000
%      Unit vector:    0.38461538   0.92307692   0.00000000
%      Magnitude  :   13.00000000
%
%      Vector     :    0.00000010   0.00000020   0.00000020
%      Unit vector:    0.33333333   0.66666667   0.66666667
%      Magnitude  :    0.00000030
%
%      Vectorized case:
%      Vector     :    5.00000000  12.00000000   0.00000000
%      Unit vector:    0.38461538   0.92307692   0.00000000
%      Magnitude  :   13.00000000
%
%      Vector     :    0.00000010   0.00000020   0.00000020
%      Unit vector:    0.33333333   0.66666667   0.66666667
%      Magnitude  :    0.00000030
%
%
%-Particulars
%
%   cspice_unorm references a function called cspice_vnorm (which itself is
%   numerically stable) to calculate the norm of the input vector `v1'.
%   If the norm is equal to zero, then each component of the output
%   vector `vout' is set to zero. Otherwise, `vout' is calculated by
%   dividing `v1' by the norm.
%
%-Exceptions
%
%   1)  If the input argument `v1' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   2)  If the input argument `v1' is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 27-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       complete code example to -Examples section.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.2, 09-NOV-2012 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.1, 21-APR-2010 (EDW)
%
%       Corrected description of call example to list all output arguments,
%       from:
%
%          vout = cspice_unorm(v1)
%
%       to
%
%          [vout, vmag] = cspice_unorm(v1)
%
%   -Mice Version 1.0.0, 25-APR-2006 (EDW)
%
%-Index_Entries
%
%   3-dimensional unit vector and norm
%
%-&

function [vout, vmag] = cspice_unorm( v1 )

   switch nargin
      case 1

         v1 = zzmice_dp(v1);

      otherwise

         error ( 'Usage: [_vout(3)_, _vmag_] = cspice_unorm(_v1(3)_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [vout, vmag] = mice('unorm_c',v1);
   catch spiceerr
      rethrow(spiceerr)
   end



