%-Abstract
%
%   CSPICE_SPKSSB returns the state (position and velocity) of a target body
%   relative to the solar system barycenter.
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
%      ref      the name of the reference frame to which the vectors returned
%               by the routine should be rotated.
%
%               [1,c1] = size(ref); char = class(ref)
%
%                  or
%
%               [1,1] = size(ref); cell = class(ref)
%
%               This may be any frame supported by the Mice frame system,
%               including dynamic and other non-inertial frames.
%
%   the call:
%
%      [starg] = cspice_spkssb( targ, et, ref )
%
%   returns:
%
%      starg    contains the position and velocity of the target body,
%               relative to the solar system barycenter, at epoch `et'.
%
%               [6,1] = size(starg); double = class(starg)
%
%               These vectors are rotated into the specified reference
%               frame. Units are always km and km/sec.
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
%   1) In the following example, cspice_spkssb is used to display
%      the distance from Earth (Body 399) to Mars (body 499) at
%      a given epoch.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: spkssb_ex1.tm
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
%            de418.bsp                     Planetary ephemeris
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de418.bsp',
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
%      function spkssb_ex1()
%
%         %
%         % Define parameters for a state lookup:
%         %
%         % Return the state vector of Mars (499)
%         % and Earth (399) with respect to the Solar System
%         % Barycenter in the J2000 frame at epoch
%         % July 4, 2003 11:00 AM PST.
%         %
%         EARTH = 399;
%         EPOCH = 'July 4, 2003 11:00 AM PST';
%         FRAME = 'J2000';
%         MARS  = 499;
%
%         %
%         % Load the required kernels.
%         %
%         cspice_furnsh( 'spkssb_ex1.tm' );
%
%         %
%         % Convert the epoch to ephemeris time.
%         %
%         [et] = cspice_str2et( EPOCH );
%
%         %
%         % Look-up the states for the defined parameters.
%         %
%         [searth] = cspice_spkssb( EARTH, et, FRAME );
%         [smars]  = cspice_spkssb( MARS, et, FRAME );
%
%         %
%         % What measure of distance separates the two bodies
%         % at epoch.
%         %
%         dist = cspice_vdist( searth(1:3), smars(1:3) );
%
%         fprintf( 'The absolute distance (km)     : %23.10f\n', dist )
%         fprintf( 'between Mars and Earth at epoch:  %s\n', EPOCH )
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
%      The absolute distance (km)     :     80854819.7017317861
%      between Mars and Earth at epoch:  July 4, 2003 11:00 AM PST
%
%
%      Note that the two cspice_spkssb calls could be replaced by
%
%         [state, lt] = cspice_spkgeo( EARTH, et, frame, MARS );
%
%      or
%
%         [state, lt] = cspice_spkezr( 'EARTH', et,     ...
%                                      frame,   'NONE', ...
%                                      'MARS'           );
%
%      using the norm of the position components of the `state'
%      vector to compute the distance between the bodies.
%
%-Particulars
%
%   In order to compute the state of one body relative to another,
%   the states of the two bodies must be known relative to a third
%   body. One simple solution is to use the solar system barycenter
%   as the third body.
%
%   Ephemeris data from more than one segment may be required
%   to determine the state of a body relative to the barycenter.
%   cspice_spkssb reads as many segments as necessary, from as many
%   files as necessary, using files that have been loaded by
%   previous calls to cspice_furnsh or cspice_spklef (load ephemeris file).
%
%-Exceptions
%
%   1)  If sufficient information has not been "loaded" via the routine
%       cspice_furnsh, cspice_spklef or the PCK kernel loaders, an error is
%       signaled by a routine in the call tree of this routine.
%
%   2)  If any of the input arguments, `targ', `et' or `ref', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `targ', `et' or `ref', is not
%       of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See -Restrictions.
%
%-Restrictions
%
%   1)  The ephemeris files to be used by cspice_spkssb must be loaded
%       by cspice_furnsh or cspice_spklef before cspice_spkssb is called.
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
%   state relative to solar system barycenter
%
%-&
function [starg] = cspice_spkssb( targ, et, ref )

   switch nargin
      case 3

         targ = zzmice_int(targ);
         et = zzmice_dp(et);
         ref = zzmice_str(ref);

      otherwise

         error ( 'Usage: [starg(6)] = cspice_spkssb( targ, et, `ref` )' )

   end

   %
   % Call the MEX library.
   %
   try
      [starg] = mice('spkssb_c', targ, et, ref);
   catch spiceerr
      rethrow(spiceerr)
   end
