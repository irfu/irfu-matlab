%-Abstract
%
%   CSPICE_DVHAT calculates the unit vector corresponding to a state or states
%   and the derivative of the unit vector.
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
%      s1       the array(s) defining a state or states;
%
%                            dr1
%                  s1 = (r1, --- ).
%                             dt
%
%               [6,n] = size(s1); double = class(s1)
%
%   the call:
%
%      [sout] = cspice_dvhat( s1 )
%
%   returns:
%
%      sout     the array(s) containing the unit vector(s) pointing in the
%               direction of the position component(s) of `s1' and the
%               derivative of the unit vector with respect to time;
%
%                             du               r1
%                  sout = [u, -- ] where u = ------
%                             dt             ||r1||
%
%               [6,n] = size(sout); double = class(sout)
%
%               `sout' returns with the same vectorization measure (N)
%               as `s1'.
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
%   1) Suppose that `state' gives the apparent state of a body with
%      respect to an observer. This routine can be used to compute the
%      instantaneous angular rate of the object across the sky as seen
%      from the observers vantage.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: dvhat_ex1.tm
%
%         This meta-kernel is intended to support operation of SPICE
%         example programs. The kernels shown here should not be
%         assumed to contain adequate or correct versions of data
%         required by SPICE-based user applications.
%
%         In order for an application to use this meta-kernel, the
%         kernels referenced here must be present in the user's
%         current working directory.
%
%         The names and contents of the kernels referenced
%         by this meta-kernel are as follows:
%
%            File name                     Contents
%            ---------                     --------
%            de421.bsp                     Planetary ephemeris
%            pck00008.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00008.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function dvhat_ex1()
%
%         target   = 'MOON';
%         frame    = 'J2000';
%         abcorr   = 'LT+S';
%         observer = 'EARTH BARYCENTER';
%
%         %
%         % Load SPK, PCK, and LSK kernels, use a meta kernel for
%         % convenience.
%         %
%         cspice_furnsh( 'dvhat_ex1.tm' );
%
%         %
%         % Define an arbitrary epoch, convert the epoch to ephemeris time.
%         %
%         EPOCH = 'Jan 1 2009';
%         et    = cspice_str2et( EPOCH );
%
%         %
%         % Calculate the state of the moon with respect to the earth-moon
%         % barycenter in J2000, corrected for light time and stellar
%         % aberration at `et'.
%         %
%         [ state, lt ] = cspice_spkezr( target, et, frame, ...
%                                           abcorr, observer       );
%
%         %
%         % Calculate the unit vector of `state' and the derivative of the
%         % unit vector.
%         %
%         ustate = cspice_dvhat( state );
%
%         %
%         % Calculate the instantaneous angular velocity from the magnitude
%         % of the derivative of the unit vector.
%         %
%         %   v = r x omega
%         %
%         %   ||omega|| = ||v||  for  r . v = 0
%         %               -----
%         %               ||r||
%         %
%         %   ||omega|| = ||v||  for  ||r|| = 1
%         %
%         omega = cspice_vnorm( ustate(4:6) );
%
%         fprintf( 'Instantaneous angular velocity (rad/sec): %18.12e\n', ...
%                                                                      omega )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Instantaneous angular velocity (rad/sec): 2.481066592694e-06
%
%
%-Particulars
%
%   Let `s1' be a state vector with position and velocity components P
%   and V respectively. From these components one can compute the
%   unit vector parallel to P, call it U and the derivative of U
%   with respect to time, DU. This pair (U,DU) is the state returned
%   by this routine in `sout'.
%
%-Exceptions
%
%   1)  If `s1' represents the zero vector, then the position
%       component of `sout' will also be the zero vector. The
%       velocity component will be the velocity component
%       of `s1'.
%
%   2)  If the input argument `s1' is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   3)  If the input argument `s1' is not of the expected type, or it
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "dvhat" to "sout".
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section. Corrected minor typos in header.
%
%       Edited the header to comply with NAIF standard. Added
%       meta-kernel to the example.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 03-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 04-MAY-2010 (EDW)
%
%-Index_Entries
%
%   state of a unit vector parallel to a state vector
%
%-&

function [sout] = cspice_dvhat( s1 )

   switch nargin
      case 1

         s1 = zzmice_dp(s1);

      otherwise

         error ( 'Usage: [_sout(6)_] = cspice_dvhat( _s1(6)_ )' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      [sout] = mice('dvhat_c',s1);
   catch spiceerr
      rethrow(spiceerr)
   end
