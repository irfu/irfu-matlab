%-Abstract
%
%   CSPICE_DVDOT returns the time derivative of the dot product of
%   two position vectors.
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
%      s1       any state vector(s).
%
%               [6,n] = size(s1); double = class(s1)
%
%                            dr1
%                  s1 = (r1, --- )
%                            dt
%
%               The components are in order (x, y, z, dx/dt, dy/dt, dz/dt)
%
%      s2       a second state vector(s);
%
%               [6,n] = size(s2); double = class(s2)
%
%                            dr2
%                  s2 = (r2, --- ).
%                            dt
%
%               An implicit assumption exists that `s1' and `s2' are specified
%               in the same reference frame. If this is not the case, the
%               numerical result has no meaning.
%
%   the call:
%
%      [dvdot] = cspice_dvdot( s1, s2 )
%
%   returns:
%
%      dvdot    the time derivative(s) of the dot product between the position
%               components of `s1' and `s2'.
%
%               [1,n] = size(dvdot); double = class(dvdot)
%
%               `dvdot' returns with the same vectorization measure (N)
%               as `s1' and `s2'.
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
%   1) Suppose that given two state vectors whose position components
%      are unit vectors, and that we need to compute the rate of
%      change of the angle between the two vectors.
%
%      Example code begins here.
%
%
%      function dvdot_ex1()
%
%         %
%         % Define the two state vectors whose position
%         % components are unit vectors.
%         %
%         s1 = [  7.2459e-01,  6.6274e-01,  1.8910e-01,                    ...
%                -1.5990e-06,  1.6551e-06,  7.4873e-07 ]';
%         s2 = [  8.4841e-01, -4.7790e-01, -2.2764e-01,                    ...
%                 1.0951e-07,  1.0695e-07,  4.8468e-08 ]';
%
%         %
%         % We know that the Cosine of the angle `theta' between them
%         % is given by
%         %
%         %    cos(theta) = dot( p1, p2 )
%         %
%         % where `p1' and `p2' are the position components of the
%         % `s1' and `s2' state vectors, respectively.
%         %
%         % Thus by the chain rule, the derivative of the angle is
%         % given by:
%         %
%         %    sin(theta) dtheta/dt = cspice_dvdot(s1,s2)
%         %
%         % Thus for values of `theta' away from zero we can compute
%         % dtheta/dt as:
%         %
%         dtheta = cspice_dvdot(s1,s2) /                                   ...
%                  sqrt( 1 - dot( s1(1:3), s2(1:3) )^2 );
%
%         fprintf( 'Rate of change of angle between S1 and S2: %17.12f\n', ...
%                                                                   dtheta )
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Rate of change of angle between S1 and S2:   -0.000002232415
%
%
%      Note that if the position components of `s1' and `s2' are parallel,
%      the derivative of the  angle between the positions does not
%      exist. Any code that computes the derivative of the angle
%      between two position vectors should account for the case
%      when the position components are parallel.
%
%-Particulars
%
%   In this discussion, the notation
%
%      < v1, v2 >
%
%   indicates the dot product of vectors `v1' and `v2'.
%
%   Given two state vectors `s1' and `s2' made up of position and velocity
%   components (r1,v1) and (r2,v2) respectively, cspice_dvdot calculates
%   the derivative of the dot product of `p1' and `p2', i.e. the time
%   derivative
%
%         d
%         -- < r1, r2 > = < v1, r2 > + < r1, v2 >
%         dt
%
%-Exceptions
%
%   1)  If any of the input arguments, `s1' or `s2', is undefined, an
%       error is signaled by the Matlab error handling system.
%
%   2)  If any of the input arguments, `s1' or `s2', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%   3)  If the input vectorizable arguments `s1' and `s2' do not have
%       the same measure of vectorization (N), an error is signaled by
%       the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The user is responsible for determining that the states `s1' and
%       `s2' are not so large as to cause numeric overflow. In most
%       cases this won't present a problem.
%
%   2)  An implicit assumption exists that `s1' and `s2' are specified in
%       the same reference frame. If this is not the case, the
%       numerical result has no meaning.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standards. Added complete
%       example code.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 20-APR-2010 (EDW)
%
%-Index_Entries
%
%   time derivative of a dot product
%
%-&

function [dvdot] = cspice_dvdot(s1, s2)

   switch nargin
      case 2

         s1 = zzmice_dp(s1);
         s2 = zzmice_dp(s2);

      otherwise

         error ( 'Usage: [_dvdot_] = cspice_dvdot(_s1(6)_, _s2(6)_)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [dvdot] = mice('dvdot_c', s1, s2);
   catch spiceerr
      rethrow(spiceerr)
   end

