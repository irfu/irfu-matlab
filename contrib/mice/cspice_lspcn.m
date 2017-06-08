%-Abstract
%
%   CSPICE_LSPCN computes L_s, the planetocentric longitude of the sun,
%   as seen from a specified body.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      body     is the name of the central body, typically a planet.
%
%               [1,c1] = size(body); char = class(body)
%
%                  or
%
%               [1,1] = size(body); cell = class(body)
%
%      et       is the epoch at which the longitude of the sun (L_s) is
%               to be computed. 'et' is expressed as seconds past J2000
%               TDB (Barycentric Dynamical Time).
%
%               [1,n] = size(et); double = class(et)
%
%      abcorr   indicates the aberration corrections to be applied
%               when computing the longitude of the sun.
%
%               [1,c2] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               'abcorr' may be any of the following.
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       Correct the position of the sun,
%                             relative to the central body, for
%                             planetary (light time) aberration.
%
%                  'LT+S'     Correct the position of the sun,
%                             relative to the central body, for
%                             planetary and stellar aberrations.
%
%   the call:
%
%      lspcn = cspice_lspcn( body, et, abcorr )
%
%   returns:
%
%      lspcn   the planetocentric longitude of the sun, often called
%              "L_s," for the specified body at the specified time.
%
%              [1,n] = size(et); double = class(et)
%
%              The longitude is defined in a right-handed frame whose
%              basis vectors are defined as follows:
%
%              - The positive Z direction is given by the instantaneous
%                angular velocity vector of the orbit of the body about
%                the sun.
%
%              - The positive X direction is that of the cross product
%                of the instantaneous north spin axis of the body with
%                the positive Z direction.
%
%              - The positive Y direction is Z x X.
%
%              Units are radians; the range is 0 to 2*pi. Longitudes are
%              positive to the east.
%
%              'lspcn' returns with the same vectorization measure (N) as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         File name: standard.tm
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
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0011.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00010.tpc',
%                                'naif0011.tls'  )
%
%         \begintext
%
%   Example(1):
%
%      Scalar
%
%      cspice_furnsh( 'standard.tm' )
%
%      et = cspice_str2et('21 march 2006')
%      lspcn = cspice_lspcn( 'earth', et, 'none' ) * cspice_dpr
%
%   MATLAB outputs:
%
%      lspcn =
%
%          0.2365
%
%   Example(2):
%
%      Vector
%
%      et = cspice_str2et('21 march 2005') + [0:1000000.:10000000.];
%      lspcn = cspice_lspcn( 'earth', et, 'none' ) * cspice_dpr
%
%   MATLAB outputs:
%
%      lspcn =
%
%          0.4815   11.9353   23.3193   34.6294   45.8611   57.0476   68.1626
%         79.2509   90.3058  101.3366  112.3833
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine lspcn_c.c.
%
%   MICE.REQ
%   ABCORR.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 10-NOV-2015, EDW (JPL)
%
%      Script rewritten to call CSPICE interface.
%
%   -Mice Version 0.9.0, 23-OCT-2015, EDW (JPL)
%
%      Pure Matlab script.
%
%-Index_Entries
%
%   planetocentric longitude of sun
%   compute L_s
%   compute Ls
%   compute L_sub_s
%
%-&

function lspcn = cspice_lspcn( body, et, abcorr )

   switch nargin
      case 3

         body   = zzmice_str(body);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);

      otherwise

         error( 'Usage: [_lspcn_] = cspice_lspcn( `body`, _et_, `abcorr`)' )

   end

   %
   % Call the MEX library, catch any error then rethrow the error from
   % this script.
   %
   try
      lspcn = mice('lspcn_c', body, et, abcorr );
   catch
      rethrow(lasterror)
   end
