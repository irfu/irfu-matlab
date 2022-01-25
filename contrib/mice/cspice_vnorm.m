%-Abstract
%
%   CSPICE_VNORM computes the magnitude of a double precision 3-dimensional
%   vector.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      v1       any double precision 3-dimensional vector(s).
%
%               [3,n] = size(v1); double = class(v1)
%
%   the call:
%
%      [vnorm] = cspice_vnorm( v1 )
%
%   returns:
%
%      vnorm    the magnitude(s) of `v1' calculated in a numerically stable
%               way.
%
%               [1,n] = size(vnorm); double = class(vnorm)
%
%               `vnorm' returns with the same vectorization measure, N,
%               as `v1'.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Define a set of 3-dimensional vectors and compute the
%      magnitude of each vector within.
%
%
%      Example code begins here.
%
%
%      function vnorm_ex1()
%
%         %
%         % Local parameters.
%         %
%         SETSIZ =   3;
%
%         %
%         % Define a set of 3-dimensional vectors.
%         %
%         v1 = [ [ 1.0,    2.0,  2.0   ]',                                 ...
%                [ 5.0,   12.0,  0.0   ]',                                 ...
%                [-5.e-17, 0.0, 12.e-17]' ];
%
%         %
%         % Calculate the magnitude of each vector
%         %
%         for i=1:SETSIZ
%
%            fprintf( 'Input vector:  %9.2e %9.2e %9.2e\n', v1(:,i) )
%            fprintf( 'Magnitude   :  %23.20f\n', cspice_vnorm( v1(:,i) ) )
%            fprintf( '\n' )
%
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Input vector:   1.00e+00  2.00e+00  2.00e+00
%      Magnitude   :   3.00000000000000000000
%
%      Input vector:   5.00e+00  1.20e+01  0.00e+00
%      Magnitude   :  13.00000000000000000000
%
%      Input vector:  -5.00e-17  0.00e+00  1.20e-16
%      Magnitude   :   0.00000000000000013000
%
%
%   2) Calculate the distance between Mercury and the Earth for a
%      period of 5 days from Jan 1, 2010 at 12:34:56.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: vnorm_ex2.tm
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
%      function vnorm_ex2()
%
%         %
%         % Set a size for the Nx1 array of ephemeris times.
%         %
%         N = 4;
%
%         %
%         %  Load a set of kernels: an SPK file, a PCK
%         %  file and a leapseconds file. Use a meta
%         %  kernel for convenience.
%         %
%         cspice_furnsh( 'vnorm_ex2.tm' )
%
%         %
%         % Set a reference epoch, convert the string representation
%         % to ET.
%         %
%         utc = 'Jan 1 2010 12:34:56';
%         et_0 = cspice_str2et( utc );
%
%         %
%         % Create an array of N elements off the reference epoch in
%         % steps of one day in ET seconds.
%         %
%         et = [1:N]*cspice_spd() + et_0;
%
%         %
%         % Calculate the geometric position of Mercury with respect to
%         % the earth, without aberration correction, at time 'et'.
%         %
%         target   = 'Mercury';
%         frame    = 'J2000';
%         abcorr   = 'none';
%         observer = 'Earth';
%
%         [pos, lt] = cspice_spkpos( target, et_0, frame,  ...
%                                       abcorr, observer );
%
%         %
%         % Calculate the  magnitude of the position vector returned
%         % from cspice_spkpos.
%         %
%         vmag = cspice_vnorm( pos );
%         utcstr = cspice_et2utc( et_0, 'C', 0 );
%
%         disp('      UTC time          Distance (km)'  )
%         disp('--------------------  ----------------' )
%         disp('Scalar:')
%         fprintf( '%s  %16.6f\n', utcstr, vmag )
%
%         [pos, lt] = cspice_spkpos( target, et, frame, abcorr, observer );
%
%         %
%         % Calculate the 1xN array of magnitudes of the N position vectors
%         % returned from cspice_spkpos.
%         %
%         vmag   = cspice_vnorm ( pos );
%         utcstr = cspice_et2utc( et, 'C', 0 );
%
%         disp('Vectorized:')
%         for i=1:4
%            fprintf( '%s  %16.6f\n', utcstr(i,:), vmag(i) )
%         end
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
%            UTC time          Distance (km)
%      --------------------  ----------------
%      Scalar:
%      2010 JAN 01 12:34:56  104117901.898226
%      Vectorized:
%      2010 JAN 02 12:34:56  102527780.780346
%      2010 JAN 03 12:34:56  101381073.443063
%      2010 JAN 04 12:34:56  100691733.096109
%      2010 JAN 05 12:34:56  100461190.556787
%
%
%-Particulars
%
%   cspice_vnorm takes care to avoid overflow while computing the norm of the
%   input vector `v1'. cspice_vnorm finds the component of `v1' whose
%   magnitude is the largest. Calling this magnitude `v1max', the norm is
%   computed using the formula:
%
%                             ||    1         ||
%      cspice_vnorm = v1max * || ------- * v1 ||
%                             ||  v1max       ||
%
%   where the notation ||x|| indicates the norm of the vector `x'.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added first complete
%       code example and the second example's problem statement and
%       meta-kernel. Updated the second code example to produce formatted
%       output.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 18-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 24-APR-2010 (EDW)
%
%-Index_Entries
%
%   norm of 3-dimensional vector
%
%-&

function [vnorm] = cspice_vnorm(v1)

   switch nargin
      case 1

         v1 = zzmice_dp(v1);

      otherwise

         error ( 'Usage: [_vnorm_] = cspice_vnorm(_v1(3)_)' )

   end

   %
   % Call the MEX library.
   %
   try
      [vnorm] = mice('vnorm_c', v1);
   catch spiceerr
      rethrow(spiceerr)
   end
