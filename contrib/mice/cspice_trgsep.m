%-Abstract
%
%   CSPICE_TRGSEP computes the angular separation in radians between two
%   spherical or point objects.
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
%      et       the time in ephemeris seconds past J2000 TDB at which the
%               separation is to be measured.
%
%               [1,1] = size(et); double = class(et)
%
%      targ1    the string naming the first body of interest.
%
%               [1,c1] = size(targ1); char = class(targ1)
%
%                  or
%
%               [1,1] = size(targ1); cell = class(targ1)
%
%               You can also supply the integer ID code for the object as
%               an integer string. For example both 'MOON' and '301' are
%               legitimate strings that indicate the moon is the target body.
%
%      shape1   the string naming the geometric model used to represent the
%               shape of the `targ1' body.
%
%               [1,c2] = size(shape1); char = class(shape1)
%
%                  or
%
%               [1,1] = size(shape1); cell = class(shape1)
%
%               Models supported by this routine:
%
%                 'SPHERE'        Treat the body as a sphere with
%                                 radius equal to the maximum value of
%                                 BODYnnn_RADII.
%
%                 'POINT'         Treat the body as a point;
%                                 radius has value zero.
%
%               The `shape1' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      frame1   the string naming the body-fixed reference frame
%               corresponding to `targ1'.
%
%               [1,c3] = size(frame1); char = class(frame1)
%
%                  or
%
%               [1,1] = size(frame1); cell = class(frame1)
%
%               cspice_trgsep does not currently use this argument's value,
%               its use is reserved for future shape models. The value 'NULL'
%               will suffice for 'POINT' and 'SPHERE' shaped bodies.
%
%      targ2    the string naming the second body of interest.
%
%               [1,c4] = size(targ2); char = class(targ2)
%
%                  or
%
%               [1,1] = size(targ2); cell = class(targ2)
%
%               You can also supply the integer ID code for the object as
%               an integer string. For example both 'MOON' and '301' are
%               legitimate strings that indicate the moon is the target body.
%
%      shape2   the string naming the geometric model used to represent the
%               shape of the `targ2'.
%
%               [1,c5] = size(shape2); char = class(shape2)
%
%                  or
%
%               [1,1] = size(shape2); cell = class(shape2)
%
%               Models supported by this routine:
%
%                 'SPHERE'        Treat the body as a sphere with
%                                 radius equal to the maximum value of
%                                 BODYnnn_RADII.
%
%                 'POINT'         Treat the body as a single point;
%                                 radius has value zero.
%
%               The `shape2' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      frame2   the string naming the body-fixed reference frame
%               corresponding to `targ2'.
%
%               [1,c6] = size(frame2); char = class(frame2)
%
%                  or
%
%               [1,1] = size(frame2); cell = class(frame2)
%
%               cspice_trgsep does not currently use this argument's value,
%               its use is reserved for future shape models. The value 'NULL'
%               will suffice for 'POINT' and 'SPHERE' shaped bodies.
%
%      obsrvr   the string naming the observing body.
%
%               [1,c7] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the observer is
%               Earth.
%
%      abcorr   the string description of the aberration corrections to apply
%               to the state evaluations to account for one-way light time
%               and stellar aberration.
%
%               [1,c8] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               This routine accepts the same aberration corrections
%               as does the SPICE routine cspice_spkezr. See the header of
%               cspice_spkezr for a detailed description of the aberration
%               correction options. For convenience, the options are
%               listed below:
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       "Reception" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     "Reception" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'CN'       "Reception" case: converged
%                             Newtonian light time correction.
%
%                  'CN+S'     "Reception" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%                  'XLT'      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'XCN'      "Transmission" case: converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    "Transmission" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               The `abcorr' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%   the call:
%
%      [trgsep] = cspice_trgsep( et,     targ1,  shape1, frame1, targ2, ...
%                                shape2, frame2, obsrvr, abcorr  )
%
%   returns:
%
%      trgsep   the angular separation between two targets, `targ1' and
%               `targ2', as seen from an observer `obsrvr' expressed in
%               radians.
%
%               [1,1] = size(trgsep); double = class(trgsep)
%
%               The observer is the angle's vertex. The angular separation
%               between the targets may be measured between the centers or
%               figures (limbs) of the targets, depending on whether the
%               target shapes are modeled as spheres or points.
%
%               If the target shape is either a spheroid or an ellipsoid, the
%               radius used to compute the limb will be the largest of the
%               radii of the target's tri-axial ellipsoid model.
%
%               If the targets are modeled as points the result ranges from 0
%               to Pi radians or 180 degrees.
%
%               If the target shapes are modeled as spheres or ellipsoids,
%               the function returns a negative value when the bodies overlap
%               (occult). Note that in this situation the function returns 0
%               when the limbs of the bodies start or finish the overlap.
%
%               The positions of the targets may optionally be corrected for
%               light time and stellar aberration.
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
%   1) Calculate the apparent angular separation of the Earth and
%      Moon as observed from the Sun at a TDB time known as a time
%      of maximum separation. Calculate and output the separation
%      modeling the Earth and Moon as point bodies and as spheres.
%      Provide the result in degrees.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: trgsep_ex1.tm
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
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00009.tpc',
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
%      function trgsep_ex1()
%
%         %
%         % Local variables.
%         %
%         frame = {'IAU_MOON','IAU_EARTH'};
%
%         targ  = {'MOON','EARTH'};
%
%         shape = {'POINT','SPHERE'};
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'trgsep_ex1.tm' );
%
%         tdbstr = '2007-JAN-11 11:21:20.213872 (TDB)';
%         obsrvr = 'SUN';
%         abcorr = 'LT+S';
%
%         [et]   = cspice_str2et( tdbstr );
%
%         value  = cspice_trgsep( et,       targ(1), shape(1),             ...
%                                 frame(1), targ(2), shape(1),             ...
%                                 frame(2), obsrvr,  abcorr    );
%
%         fprintf( 'Bodies:          %-6s%-6s\n',                          ...
%                    char(targ(1)), char(targ(2)) )
%         fprintf( 'as seen from:    %-6s\n', obsrvr )
%         fprintf( 'at TDB time:     %-36s\n', tdbstr )
%         fprintf( 'with correction: %s\n', abcorr )
%         fprintf( '\n' )
%
%         fprintf( 'Apparent angular separation:\n' )
%         fprintf( '   point body models  (deg.):  %11.8f\n',              ...
%                                          value * cspice_dpr )
%
%         value  = cspice_trgsep( et,       targ(1), shape(2),             ...
%                                 frame(1), targ(2), shape(2),             ...
%                                 frame(2), obsrvr,  abcorr    );
%
%         fprintf( '   sphere body models (deg.):  %11.8f\n',              ...
%                                          value * cspice_dpr )
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
%      Bodies:          MOON  EARTH
%      as seen from:    SUN
%      at TDB time:     2007-JAN-11 11:21:20.213872 (TDB)
%      with correction: LT+S
%
%      Apparent angular separation:
%         point body models  (deg.):   0.15729276
%         sphere body models (deg.):   0.15413221
%
%
%-Particulars
%
%   This routine determines the apparent separation between the
%   two objects as observed from a third. The value reported is
%   corrected for light time. Moreover, if at the time this routine
%   is called, stellar aberration corrections are enabled, this
%   correction will also be applied to the apparent positions of the
%   centers of the two objects.
%
%   Please refer to the Aberration Corrections Required Reading
%   (abcorr.req) for detailed information describing the nature and
%   calculation of the applied corrections.
%
%-Exceptions
%
%   1)  If the three objects `targ1', `targ2' and `obsrvr' are not
%       distinct, an error is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If the object names for `targ1', `targ2' or `obsrvr' cannot resolve
%       to a NAIF body ID, an error is signaled by a routine in the
%       call tree of this routine.
%
%   3)  If the reference frame associated with `targ1', `frame1', is not
%       centered on `targ1', or if the reference frame associated with
%       `targ2', `frame2', is not centered on `targ2', an error is signaled
%       by a routine in the call tree of this routine. This
%       restriction does not apply to shapes 'SPHERE' and 'POINT', for
%       which the frame input is ignored.
%
%   4)  If the frame name for `frame1' or `frame2' cannot resolve to a
%       NAIF frame ID, an error is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the body shape for `targ1', `shape1', or the body shape for
%       `targ2', `shape2', is not recognized, an error is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If the requested aberration correction `abcorr' is not
%       recognized, an error is signaled by a routine in the call tree
%       of this routine.
%
%   7)  If either one or both targets' shape is modeled as sphere, and
%       the required PCK data has not been loaded, an error is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the ephemeris data required to perform the needed state
%       look-ups are not loaded, an error is signaled by a routine in
%       the call tree of this routine.
%
%   9)  If the observer `obsrvr' is located within either one of the
%       targets, an error is signaled by a routine in the call tree of
%       this routine.
%
%   10) If an error is signaled, the function returns a meaningless
%       result.
%
%   11) If any of the input arguments, `et', `targ1', `shape1',
%       `frame1', `targ2', `shape2', `frame2', `obsrvr' or `abcorr',
%       is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   12) If any of the input arguments, `et', `targ1', `shape1',
%       `frame1', `targ2', `shape2', `frame2', `obsrvr' or `abcorr',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   The following data are required:
%
%   -  An SPK file (or files) containing ephemeris data sufficient to
%      compute the position of each of the targets with respect to the
%      observer. If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%
%   -  A PCK file containing the targets' tri-axial ellipsoid model,
%      if the targets are modeled as spheres.
%
%   -  If non-inertial reference frames are used, then PCK files,
%      frame kernels, C-kernels, and SCLK kernels may be needed.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   ABCORR.REQ
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   M. Costa Sitja      (JPL)
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 07-AUG-2021 (MCS) (JDR)
%
%-Index_Entries
%
%   compute the angular separation between two target bodies
%
%-&
function [trgsep] = cspice_trgsep( et, targ1, shape1, frame1, targ2, ...
                                   shape2, frame2, obsrvr, abcorr )

   switch nargin
      case 9

         et     = zzmice_dp(et);
         targ1  = zzmice_str(targ1);
         shape1 = zzmice_str(shape1);
         frame1 = zzmice_str(frame1);
         targ2  = zzmice_str(targ2);
         shape2 = zzmice_str(shape2);
         frame2 = zzmice_str(frame2);
         obsrvr = zzmice_str(obsrvr);
         abcorr = zzmice_str(abcorr);

      otherwise

         error ( [ 'Usage: [trgsep] = '                                     ...
                   'cspice_trgsep( et, `targ1`, `shape1`, `frame1`, '       ...
                   '`targ2`, `shape2`, `frame2`, `obsrvr`, `abcorr` )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [trgsep] = mice('trgsep_c', et, targ1, shape1, frame1, targ2, shape2, ...
                      frame2, obsrvr, abcorr);
   catch spiceerr
      rethrow(spiceerr)
   end
