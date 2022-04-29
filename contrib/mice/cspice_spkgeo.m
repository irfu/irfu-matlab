%-Abstract
%
%   CSPICE_SPKGEO computes the geometric state (position and velocity) of a
%   target body relative to an observing body.
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
%      targ     the standard NAIF ID code for a target body.
%
%               [1,1] = size(targ); int32 = class(targ)
%
%      et       the epoch (ephemeris time) at which the state of the target
%               body is to be computed.
%
%               [1,1] = size(et); double = class(et)
%
%      ref      the name of the reference frame to which the state vector
%               returned by the routine should be rotated.
%
%               [1,c1] = size(ref); char = class(ref)
%
%                  or
%
%               [1,1] = size(ref); cell = class(ref)
%
%               This may be any frame supported by the SPICELIB function
%               FRMCHG. See also the Frames Required Reading for a list of
%               supported frames.
%
%      obs      the standard NAIF ID code for an observing body.
%
%               [1,1] = size(obs); int32 = class(obs)
%
%   the call:
%
%      [state, lt] = cspice_spkgeo( targ, et, ref, obs )
%
%   returns:
%
%      state    contains the geometric position and velocity of the target
%               body, relative to the observing body, at epoch `et'.
%
%               [6,1] = size(state); double = class(state)
%
%               `state' has six elements: the first three contain the
%               target's position; the last three contain the target's
%               velocity. These vectors are transformed into the specified
%               reference frame.
%
%               Units are always km and km/sec.
%
%      lt       the one-way light time from the observing body to the
%               geometric position of the target body in seconds at the
%               specified epoch.
%
%               [1,1] = size(lt); double = class(lt)
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
%   1) Return the geometric state vector of Mars (499) as seen from
%      Earth (399) in the J2000 frame and the one-way light time
%      between them at the epoch July 4, 2003 11:00 AM PST.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: spkgeo_ex1.tm
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
%            File name                        Contents
%            ---------                        --------
%            de430.bsp                        Planetary ephemeris
%            mar097.bsp                       Mars satellite ephemeris
%            naif0011.tls                     Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de430.bsp',
%                                'mar097.bsp',
%                                'naif0011.tls' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function spkgeo_ex1()
%
%         %
%         % Load a set of kernels. Use a meta
%         % kernel for convenience.
%         %
%         cspice_furnsh( 'spkgeo_ex1.tm' );
%
%         %
%         % Define parameters for a state lookup.
%         %
%         target = 499;
%         epoch  = 'July 4, 2003 11:00 AM PST';
%         reffrm = 'J2000';
%         obsrvr = 399;
%
%         %
%         % Convert the epoch to ephemeris time.
%         %
%         [et] = cspice_str2et( epoch );
%
%         %
%         % Look-up the state for the defined parameters.
%         %
%         [state, lt] = cspice_spkgeo( target, et, reffrm, obsrvr );
%
%         %
%         % Output...
%         %
%         fprintf( 'The position of    : %2d\n', target )
%         fprintf( 'As observed from   : %2d\n', obsrvr )
%         fprintf( 'In reference frame :  %s\n', reffrm )
%         fprintf( 'At epoch           :  %s\n', epoch )
%         fprintf( ' \n' )
%
%         %
%         % The first three entries of state contain the
%         % X, Y, Z position components. The final three contain
%         % the Vx, Vy, Vz velocity components.
%         %
%         fprintf( 'R   (km): %17.5f %17.5f %17.5f\n', ...
%                         state(1), state(2), state(3) )
%         fprintf( 'V (km/s): %17.5f %17.5f %17.5f\n', ...
%                         state(4), state(5), state(6) )
%         fprintf( ' \n' )
%         fprintf( [ 'Light time (s) between observer and target: ', ...
%                    ' %18.13f\n' ], lt                              )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      The position of    : 499
%      As observed from   : 399
%      In reference frame :  J2000
%      At epoch           :  July 4, 2003 11:00 AM PST
%
%      R   (km):    73826216.43529   -27128030.73241   -18741973.86829
%      V (km/s):          -6.80950           7.51381           3.00129
%
%      Light time (s) between observer and target:   269.7026477631753
%
%
%-Particulars
%
%   cspice_spkgeo computes the geometric state, targ(t), of the target
%   body and the geometric state, obs(t), of the observing body
%   relative to the first common center of motion. Subtracting
%   obs(t) from targ(t) gives the geometric state of the target
%   body relative to the observer.
%
%
%      center ----- obs(t)
%          |      /
%          |     /
%          |    /
%          |   /  targ(t) - obs(t)
%          |  /
%        targ(t)
%
%
%   The one-way light time, tau, is given by
%
%
%             | targ(t) - obs(t) |
%      tau = ----------------------
%                      C
%
%
%   For example, if the observing body is -94, the Mars Observer
%   spacecraft, and the target body is 401, Phobos, then the
%   first common center is probably 4, the Mars Barycenter.
%   obs(t) is the state of -94 relative to 4 and targ(t) is the
%   state of 401 relative to 4.
%
%   The center could also be the Solar System Barycenter, body 0.
%   For example, if the observer is 399, Earth, and the target
%   is 299, Venus, then obs(t) would be the state of 399 relative
%   to 0 and targ(t) would be the state of 299 relative to 0.
%
%   Ephemeris data from more than one segment may be required
%   to determine the states of the target body and observer
%   relative to a common center. cspice_spkgeo reads as many segments
%   as necessary, from as many files as necessary, using files
%   that have been loaded by previous calls to cspice_furnsh or
%   cspice_spklef (load ephemeris file).
%
%   cspice_spkgeo is similar to cspice_spkez but returns geometric states
%   only, with no option to make planetary (light-time) nor
%   stellar aberration corrections. The geometric states
%   returned by cspice_spkez and cspice_spkgeo are the same.
%
%-Exceptions
%
%   1)  If insufficient ephemeris data has been loaded to compute
%       the necessary states, the error SPICE(SPKINSUFFDATA) is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If any of the input arguments, `targ', `et', `ref' or `obs',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   3)  If any of the input arguments, `targ', `et', `ref' or `obs',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See -Restrictions.
%
%-Restrictions
%
%   1)  The ephemeris files to be used by cspice_spkgeo must be loaded
%       by cspice_furnsh or cspice_spklef before cspice_spkgeo is called.
%
%-Required_Reading
%
%   MICE.REQ
%   SPK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   geometric state of one body relative to another
%
%-&
function [state, lt] = cspice_spkgeo( targ, et, ref, obs )

   switch nargin
      case 4

         targ = zzmice_int(targ);
         et = zzmice_dp(et);
         ref = zzmice_str(ref);
         obs = zzmice_int(obs);

      otherwise

         error ( ['Usage: [state(6), lt] = ',               ...
                  'cspice_spkgeo( targ, et, `ref`, obs )'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [state_s] = mice('spkgeo_s', targ, et, ref, obs);
      state  = reshape( [state_s.state ], 6, [] );
      lt     = reshape( [state_s.lt    ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end
