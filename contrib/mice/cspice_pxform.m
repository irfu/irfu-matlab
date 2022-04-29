%-Abstract
%
%   CSPICE_PXFORM returns the matrix that transforms position
%   vectors from one specified frame to another at a specified epoch.
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
%      from     the name of a reference frame in which a position vector is
%               known.
%
%               [1,c1] = size(from); char = class(from)
%
%                  or
%
%               [1,1] = size(from); cell = class(from)
%
%      to       the name of a reference frame in which it is desired to
%               represent a position vector.
%
%               [1,c2] = size(to); char = class(to)
%
%                  or
%
%               [1,1] = size(to); cell = class(to)
%
%      et       the epoch(s) in ephemeris seconds past the epoch of J2000
%               (TDB) at which the position transformation matrix `rotate'
%               should be evaluated.
%
%               [1,n] = size(et); double = class(et)
%
%   the call:
%
%      [rotate] = cspice_pxform( from, to, et )
%
%   returns:
%
%      rotate   the matri(x|ces) that transforms position vectors from the
%               reference frame `from' to the frame `to' at epoch `et'.
%
%               If [1,1] = size(et) then [3,3]   = size(rotate)
%               If [1,n] = size(et) then [3,3,n] = size(rotate)
%                                         double = class(rotate)
%
%               If (x, y, z) is a position relative to the frame `from'
%               then the vector ( x', y', z') is the same position relative
%               to the frame `to' at epoch `et'. Here the vector (x', y', z')
%               is defined by the equation:
%
%                  .-   -.     .-        -.   .-  -.
%                  | x'  |     |          |   | x  |
%                  | y'  |  =  |  rotate  |   | y  |
%                  | z'  |     |          |   | z  |
%                  `-   -'     `-        -'   `-  -'
%
%               `rotate' returns with the same vectorization measure, N,
%               as `et'.
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
%   1) Output the right ascension and declination of the earth's pole
%      in the J2000 frame approximately every month for the time
%      interval January 1, 1990 to January 1, 1991 (UTC).
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: pxform_ex1.tm
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
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0009.tls'
%                                'pck00009.tpc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function pxform_ex1()
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'pxform_ex1.tm' )
%
%         %
%         % Define the time bounds for the time interval,
%         % 20 years, convert to ephemeris time J2000.
%         %
%         utc_bounds = strvcat( '1 Jan 1990', '1 Jan 1991' );
%         et_bounds  = cspice_str2et( utc_bounds );
%
%         %
%         % Step in units of a month.
%         %
%         step = (et_bounds(2) - et_bounds(1) ) / 12.;
%
%         %
%         % Create an array of 12 ephemeris times ending at
%         % ~et_bound(2) in intervals of 'step'.
%         %
%         et = [1:12]*step + et_bounds(1);
%
%         %
%         % Set the conversion constant "radians to degrees."
%         %
%         r2d = cspice_dpr;
%
%         %
%         % Convert the 12-vector of 'et' to an array of corresponding
%         % transformation matrices (dimensions (3,3,12) ).
%         %
%         mat = cspice_pxform( 'IAU_EARTH', 'J2000', et );
%
%         %
%         % Extract the pole vector from the transformation matrix,
%         % convert to RA and DEC expressed in degrees.
%         %
%
%         %
%         % The last column in each matrix is the pole vector (z = (0,0,1))
%         % of the earth in IAU expressed in J2000.
%         %
%         % Recall, MATLAB uses 1 based indexing, so (:,3,:) represents.
%         % the third column of the matrices.
%         %
%         pole = mat(:,3,:);
%
%         %
%         % 'pole' ready for use in cspice_radrec.
%         %
%         [radius, ra, dec] = cspice_recrad( pole );
%
%         %
%         % Output the ephemeris time and the corresponding
%         % angular values (in degrees). 'ra' and 'dec' return
%         % as double precision 12-vectors.
%         %
%         ra  = ra  * r2d;
%         dec = dec * r2d;
%
%         %
%         % Create an array of values for output.
%         %
%         output = [  et; ra; dec ];
%
%         fprintf( '%24.8f %16.8f %16.8f\n', output );
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%           -312947942.73273718     180.06356619      89.94476386
%           -310319942.64940447     180.06303239      89.94522771
%           -307691942.56607175     180.06249859      89.94569156
%           -305063942.48273909     180.06196478      89.94615541
%           -302435942.39940637     180.06143098      89.94661925
%           -299807942.31607366     180.06089718      89.94708310
%           -297179942.23274094     180.06036338      89.94754695
%           -294551942.14940822     180.05982958      89.94801080
%           -291923942.06607556     180.05929578      89.94847465
%           -289295941.98274285     180.05876198      89.94893850
%           -286667941.89941013     180.05822818      89.94940235
%           -284039941.81607741     180.05769438      89.94986620
%
%
%-Particulars
%
%   This routine provides the user level interface to computing
%   position transformations from one reference frame to another.
%
%   Note that the reference frames may be inertial or non-inertial.
%   However, the user must take care that sufficient SPICE kernel
%   information is loaded to provide a complete position
%   transformation path from the `from' frame to the `to' frame.
%
%   A common type of reference frame transformation is one from one
%   time-dependent frame to another, where the orientations of the
%   frames are computed at different epochs. For example, a remote
%   sensing application may compute the transformation from a target
%   body-fixed frame, with its orientation evaluated at the epoch of
%   photon emission, to a spacecraft instrument frame, with its
%   orientation evaluated at the epoch of photon reception. The
%   Mice routine cspice_pxfrm2 computes this type of frame
%   transformation.
%
%-Exceptions
%
%   1)  If sufficient information has not been supplied via loaded
%       SPICE kernels to compute the transformation between the
%       two frames, an error is signaled by a routine
%       in the call tree of this routine.
%
%   2)  If either frame `from' or `to' is not recognized, the error
%       SPICE(UNKNOWNFRAME) is signaled by a routine in the call tree
%       of this routine.
%
%   3)  If any of the input arguments, `from', `to' or `et', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `from', `to' or `et', is not of
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
%   MICE.REQ
%   ROTATION.REQ
%   FRAMES.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR) (NJB)
%
%       Edited the header to comply with NAIF standard. Reduced
%       the time interval and ephemeris epochs used as input in code
%       example.
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
%       Updated -I/O to extend "rotate" description and -Particulars
%       to mention cspice_pxfrm2.
%
%   -Mice Version 1.0.2, 05-FEB-2015 (EDW)
%
%       Minor edits to -I/O header to match corresponding header
%       in cspice_sxform.
%
%   -Mice Version 1.0.1, 09-NOV-2012 (EDW) (SCK)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   Find a position transformation matrix
%
%-&

function [rotate] = cspice_pxform(from, to, et)

   switch nargin
      case 3

         from = zzmice_str(from);
         to   = zzmice_str(to);
         et   = zzmice_dp(et);

      otherwise

         error ( [ 'Usage: [_rotate(3,3)_] = ' ...
                   'cspice_pxform( `from`, `to`, _et_ )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [rotate] = mice('pxform_c',from,to,et);
   catch spiceerr
      rethrow(spiceerr)
   end




