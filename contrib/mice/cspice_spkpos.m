%-Abstract
%
%   CSPICE_SPKPOS returns the position of a target body relative to an
%   observing body, optionally corrected for light time (planetary
%   aberration) and stellar aberration.
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
%      targ     the name of a target body.
%
%               [1,c1] = size(targ); char = class(targ)
%
%                  or
%
%               [1,1] = size(targ); cell = class(targ)
%
%               Optionally, you may supply the integer ID code for the
%               object as an integer string. For example both 'MOON' and
%               '301' are legitimate strings that indicate the moon is the
%               target body.
%
%               The target and observer define a position vector
%               which points from the observer to the target.
%
%      et       the ephemeris time(s), expressed as seconds past J2000 TDB, at
%               which the position of the target body relative to the
%               observer is to be computed.
%
%               [1,n] = size(et); double = class(et)
%
%               `et' refers to time at the observer's location.
%
%      ref      the name of the reference frame relative to which the output
%               position vector should be expressed.
%
%               [1,c2] = size(ref); char = class(ref)
%
%                  or
%
%               [1,1] = size(ref); cell = class(ref)
%
%               This may be any frame supported by the SPICE system,
%               including built-in frames (documented in the Frames Required
%               Reading) and frames defined by a loaded frame kernel (FK).
%
%               When `ref' designates a non-inertial frame, the
%               orientation of the frame is evaluated at an epoch
%               dependent on the selected aberration correction. See
%               the description of the output position vector `ptarg'
%               for details.
%
%      abcorr   indicates the aberration corrections to be applied to the
%               position of the target body to account for one-way light time
%               and stellar aberration.
%
%               [1,c3] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               See the discussion in the -Particulars section for
%               recommendations on how to choose aberration corrections.
%
%               `abcorr' may be any of the following:
%
%                  'NONE'     Apply no correction. Return the
%                             geometric position of the target body
%                             relative to the observer.
%
%               The following values of `abcorr' apply to the
%               "reception" case in which photons depart from the
%               target's location at the light-time corrected epoch
%               et-lt and *arrive* at the observer's location at `et':
%
%                  'LT'       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             yields the position of the target at
%                             the moment it emitted photons arriving
%                             at the observer at `et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation (see -Particulars for details).
%                             The solution invoked by the 'LT' option
%                             uses one iteration.
%
%                  'LT+S'     Correct for one-way light time and
%                             stellar aberration using a Newtonian
%                             formulation. This option modifies the
%                             position obtained with the 'LT' option
%                             to account for the observer's velocity
%                             relative to the solar system
%                             barycenter. The result is the apparent
%                             position of the target---the position
%                             as seen by the observer.
%
%                  'CN'       Converged Newtonian light time
%                             correction. In solving the light time
%                             equation, the 'CN' correction iterates
%                             until the solution converges (three
%                             iterations on all supported platforms).
%                             Whether the 'CN+S' solution is
%                             substantially more accurate than the
%                             'LT' solution depends on the geometry
%                             of the participating objects and on the
%                             accuracy of the input data. In all
%                             cases this routine will execute more
%                             slowly when a converged solution is
%                             computed. See the -Particulars section
%                             below for a discussion of precision of
%                             light time corrections.
%
%                  'CN+S'     Converged Newtonian light time
%                             correction and stellar aberration
%                             correction.
%
%
%               The following values of `abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at `et' and arrive at the
%               target's location at the light-time corrected epoch
%               et+lt:
%
%                  'XLT'      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation. This correction yields the
%                             position of the target at the moment it
%                             receives photons emitted from the
%                             observer's location at `et'.
%
%                  'XLT+S'    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation. This option modifies the
%                             position obtained with the 'XLT' option
%                             to account for the observer's velocity
%                             relative to the solar system
%                             barycenter. The computed target
%                             position indicates the direction that
%                             photons emitted from the observer's
%                             location must be "aimed" to hit the
%                             target.
%
%                  'XCN'      "Transmission" case: converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    "Transmission" case: converged
%                             Newtonian light time correction and
%                             stellar aberration correction.
%
%
%               Neither special nor general relativistic effects are
%               accounted for in the aberration corrections applied
%               by this routine.
%
%               Case and blanks are not significant in the string
%               `abcorr'.
%
%      obs      the name of an observing body.
%
%               [1,c4] = size(obs); char = class(obs)
%
%                  or
%
%               [1,1] = size(obs); cell = class(obs)
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the observer is
%               Earth.
%
%   the call:
%
%      [ptarg, lt] = cspice_spkpos( targ, et, ref, abcorr, obs )
%
%   returns:
%
%      ptarg    a Cartesian 3-vector(s) representing the position of the
%               target body relative to the specified observer.
%
%               [3,n] = size(ptarg); double = class(ptarg)
%
%               `ptarg' is corrected for the specified aberrations, and is
%               expressed with respect to the reference frame specified by
%               `ref'. The three components of `ptarg' represent the x-, y-
%               and z-components of the target's position.
%
%               `ptarg' points from the observer's location at `et' to
%               the aberration-corrected location of the target.
%               Note that the sense of this position vector is
%               independent of the direction of radiation travel
%               implied by the aberration correction.
%
%               Units are always km.
%
%               Non-inertial frames are treated as follows: letting
%               `ltcent' be the one-way light time between the observer
%               and the central body associated with the frame, the
%               orientation of the frame is evaluated at et-ltcent,
%               et+ltcent, or `et' depending on whether the requested
%               aberration correction is, respectively, for received
%               radiation, transmitted radiation, or is omitted.
%               `ltcent' is computed using the method indicated by
%               `abcorr'.
%
%      lt       the one-way light time(s) between the observer and target in
%               seconds.
%
%               [1,n] = size(lt); double = class(lt)
%
%               If the target position is corrected for aberrations, then
%               `lt' is the one-way light time between the observer and the
%               light time corrected target location.
%
%               `ptarg' and `lt' return with the same vectorization
%               measure, N, as `et'.
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
%   1) Load a planetary SPK, and look up the position of Mars
%      as seen from the Earth in the J2000 frame with aberration
%      corrections 'LT+S' (ligth time plus stellar aberration) at
%      different epochs.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: spkpos_ex1.tm
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
%      function spkpos_ex1()
%
%         %
%         %  Load a set of kernels: an SPK file, a PCK
%         %  file and a leapseconds file. Use a meta
%         %  kernel for convenience.
%         %
%         cspice_furnsh( 'spkpos_ex1.tm' )
%
%         %
%         %  Define parameters for a position lookup:
%         %
%         %  Return the position vector of Mars (499) as seen from
%         %  Earth (399) in the J2000 frame
%         %  using aberration correction LT+S (light time plus
%         %  stellar aberration) at the epoch
%         %  July 4, 2003 11:00 AM PST.
%         %
%         target   = 'Mars';
%         epoch    = 'July 4, 2003 11:00 AM PST';
%         frame    = 'J2000';
%         abcorr   = 'LT+S';
%         observer = 'Earth';
%
%         %
%         %  Convert the epoch to ephemeris time.
%         %
%         et = cspice_str2et( epoch );
%
%         %
%         %  Look-up the position for the defined parameters.
%         %
%         [ptarg, lt] = cspice_spkpos( target, et, frame, ...
%                                      abcorr, observer );
%
%         %
%         %  Output...
%         %
%         txt = sprintf( 'The position of    : %s', target);
%         disp( txt )
%
%         txt = sprintf( 'As observed from   : %s', observer );
%         disp( txt )
%
%         txt = sprintf( 'In reference frame : %s', frame );
%         disp( txt )
%         disp( ' ' )
%
%         txt = sprintf( 'Scalar' );
%         disp( txt )
%
%         utc_epoch = cspice_et2utc( et, 'C', 3 );
%
%         txt = sprintf(  'At epoch           : %s', epoch );
%         disp( txt )
%
%         txt = sprintf(  '                   : i.e. %s', utc_epoch );
%         disp( txt )
%
%         txt = sprintf( ['R (kilometers)     : ' ...
%                         '%12.4f %12.4f %12.4f'], ptarg );
%         disp( txt )
%
%         txt = sprintf( 'Light time (secs)  : %12.7f', lt );
%         disp( txt )
%
%         disp(' between observer' )
%         disp(' and target' )
%         disp( ' ' )
%
%         %
%         % Create a vector of et's, starting at `epoch'
%         % in steps of 100000 ephemeris seconds.
%         %
%         vec_et = [0:4]*100000. + et;
%
%         disp( 'Vector' )
%         vec_epoch = cspice_et2utc( vec_et, 'C', 3 );
%
%         %
%         % Look up the position vectors and light time values
%         % `lt'  corresponding to the vector of input
%         % ephemeris time `vec_et'.
%         %
%         [ptarg , lt] = cspice_spkpos( target, vec_et, ...
%                                       frame, abcorr, observer );
%
%         for i=1:5
%
%            txt = sprintf(  'At epoch (UTC)     : %s', vec_epoch(i,:) );
%            disp( txt )
%
%            txt = sprintf( ['R (kilometers)     : ' ...
%                            '%12.4f %12.4f %12.4f'], ptarg(i) );
%            disp( txt )
%
%            txt = sprintf( ['Light time (secs)  : ' ...
%                           '%12.7f'], lt(i) );
%            disp( txt )
%
%            disp(' between observer' )
%            disp(' and target' )
%            disp( ' ' )
%
%         end
%
%         %
%         %  It's always good form to unload kernels after use,
%         %  particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      The position of    : Mars
%      As observed from   : Earth
%      In reference frame : J2000
%
%      Scalar
%      At epoch           : July 4, 2003 11:00 AM PST
%                         : i.e. 2003 JUL 04 19:00:00.000
%      R (kilometers)     : 73822235.3312 -27127919.1784 -18741306.2848
%      Light time (secs)  :  269.6898816
%       between observer
%       and target
%
%      Vector
%      At epoch (UTC)     : 2003 JUL 04 19:00:00.000
%      R (kilometers)     : 73822235.3312
%      Light time (secs)  :  269.6898816
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 05 22:46:40.000
%      R (kilometers)     : -27127919.1784
%      Light time (secs)  :  266.5640396
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 07 02:33:20.000
%      R (kilometers)     : -18741306.2848
%      Light time (secs)  :  263.4803536
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 08 06:20:00.000
%      R (kilometers)     : 73140185.4372
%      Light time (secs)  :  260.4395237
%       between observer
%       and target
%
%      At epoch (UTC)     : 2003 JUL 09 10:06:40.000
%      R (kilometers)     : -26390524.9551
%      Light time (secs)  :  257.4422004
%       between observer
%       and target
%
%
%-Particulars
%
%   A sister version of this routine exists named mice_spkpos that returns
%   the output arguments as fields in a single structure.
%
%   This routine is part of the user interface to the SPICE ephemeris
%   system. It allows you to retrieve position information for any
%   ephemeris object relative to any other in a reference frame that
%   is convenient for further computations.
%
%   Please refer to the Aberration Corrections Required Reading
%   abcorr.req for detailed information describing the nature and
%   calculation of the applied corrections.
%
%-Exceptions
%
%   1)  If name of target or observer cannot be translated to its NAIF
%       ID code, the error SPICE(IDCODENOTFOUND) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the reference frame `ref' is not a recognized reference
%       frame, the error SPICE(UNKNOWNFRAME) is signaled by a routine
%       in the call tree of this routine.
%
%   3)  If the loaded kernels provide insufficient data to compute the
%       requested position vector, an error is signaled by a routine
%       in the call tree of this routine.
%
%   4)  If an error occurs while reading an SPK or other kernel file,
%       the error  is signaled by a routine in the call tree
%       of this routine.
%
%   5)  If any of the input arguments, `targ', `et', `ref', `abcorr'
%       or `obs', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   6)  If any of the input arguments, `targ', `et', `ref', `abcorr'
%       or `obs', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   This routine computes positions using SPK files that have been loaded
%   into the SPICE system, normally via the kernel loading interface routine
%   cspice_furnsh. See the routine cspice_furnsh and the SPK and KERNEL
%   Required Reading for further information on loading (and unloading)
%   kernels.
%
%   If the output position `ptarg' is to be expressed relative to a
%   non-inertial frame, or if any of the ephemeris data used to
%   compute `ptarg' are expressed relative to a non-inertial frame in
%   the SPK files providing those data, additional kernels may be
%   needed to enable the reference frame transformations required to
%   compute the position. These additional kernels may be C-kernels,
%   PCK files or frame kernels. Any such kernels must already be
%   loaded at the time this routine is called.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   ABCORR.REQ
%   FRAMES.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   B.V. Semenov        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 02-NOV-2021 (EDW) (JDR)
%
%       Changed output argument name "pos" to "ptarg".
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and meta-kernel.
%
%       Extended -I/O section.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%       Updated -Particulars to refer to Aberration Corrections
%       Required Reading document, which was added to
%       -Required_Reading list.
%
%   -Mice Version 1.0.4, 03-DEC-2014 (EDW)
%
%       Corrections made to -Version section numbering. 07-NOV-2013
%       notation now numbered as 1.0.2, and Version 1.0.3, 03-JUL-2014.
%
%       Corrections made to author identifiers for Version 1.0.3,
%       03-JUL-2014, and Version 1.0.2, 07-NOV-2013 to indicate institution.
%
%   -Mice Version 1.0.3, 03-JUL-2014 (NJB) (BVS) (EDW)
%
%       Discussion of light time corrections was updated. Assertions
%       that converged light time corrections are unlikely to be
%       useful were removed.
%
%   -Mice Version 1.0.2, 07-NOV-2013 (EDW)
%
%       Added aberration algorithm explanation to -Particulars section.
%
%   -Mice Version 1.0.1, 22-DEC-2008 (EDW)
%
%       Header edits performed to improve argument descriptions.
%       These descriptions should now closely match the descriptions
%       in the corresponding CSPICE routine.
%
%       Corrected typo in -I/O section. Replaced the "ptarg"
%       return argument name with "pos."
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   using body names get position relative to an observer
%   get position relative observer corrected for aberrations
%   read ephemeris data
%   read trajectory data
%
%-&

function [ptarg, lt] = cspice_spkpos(targ, et, ref, abcorr, obs)

   switch nargin
      case 5

         targ   = zzmice_str(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_str(obs);

      otherwise

         error ( ['Usage: [_ptarg(3)_, _lt_] = ' ...
                  'cspice_spkpos( `targ`, _et_, `ref`, `abcorr`, `obs`)'] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [ptarg_s] = mice('spkpos_s', targ, et,ref, abcorr, obs);
      ptarg     = reshape( [ptarg_s.pos], 3, [] );
      lt        = reshape( [ptarg_s.lt],  1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end
