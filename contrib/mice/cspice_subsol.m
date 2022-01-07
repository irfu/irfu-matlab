%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   cspice_subslr. This routine is supported for purposes of
%   backward compatibility only.
%
%   CSPICE_SUBSOL determines the coordinates of the sub-solar
%   point on a target  body as seen by a specified observer at a
%   specified epoch, optionally corrected for planetary (light time)
%   and stellar aberration.
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
%      method   a short string specifying the computation method to be used.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               The choices are:
%
%                  'Near point'       The sub-solar point is defined
%                                     as the nearest point on the
%                                     target to the sun.
%
%                  'Intercept'        The sub-observer point is
%                                     defined as the target surface
%                                     intercept of the line
%                                     containing the target's center
%                                     and the sun's center.
%
%               In both cases, the intercept computation treats the
%               surface of the target body as a triaxial ellipsoid.
%               The ellipsoid's radii must be available in the kernel
%               pool.
%
%               Neither case nor white space are significant in
%               `method'. For example, the string ' NEARPOINT' is
%               valid.
%
%      target   the name of the target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               `target' is case-insensitive, and leading and trailing
%               blanks in `target' are not significant. Optionally, you may
%               supply a string containing the integer ID code for the
%               object. For example both 'MOON' and '301' are legitimate
%               strings that indicate the moon is the target body.
%
%               This routine assumes that the target body is modeled
%               by a tri-axial ellipsoid, and that a PCK file
%               containing its radii has been loaded into the kernel
%               pool via cspice_furnsh.
%
%      et       the epoch(s) in ephemeris seconds past J2000 at which the
%               sub-solar point on the target body is to be computed.
%
%               [1,1] = size(et); double = class(et)
%
%      abcorr   indicates the aberration corrections to be applied when
%               computing the observer-target state.
%
%               [1,c3] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               abcorr may be any of the following.
%
%                  'NONE'     Apply no correction. Return the
%                             geometric sub-solar point on the target
%                             body.
%
%                  'LT'       Correct for planetary (light time)
%                             aberration. Both the state and rotation
%                             of the target body are corrected for one
%                             way light time from target to observer.
%
%                             The state of the sun relative to the
%                             target is corrected for one way light
%                             from the sun to the target; this state
%                             is evaluated at the epoch obtained by
%                             retarding `et' by the one way light time
%                             from target to observer.
%
%                  'LT+S'     Correct for planetary (light time) and
%                             stellar aberrations. Light time
%                             corrections are the same as in the 'LT'
%                             case above. The target state is
%                             additionally corrected for stellar
%                             aberration as seen by the observer, and
%                             the sun state is corrected for stellar
%                             aberration as seen from the target.
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
%                             of cspice_spkezr for a discussion of precision
%                             of light time corrections. Light time
%                             corrections are applied as in the 'LT'
%                             case.
%
%                  'CN+S'     Converged Newtonian light time
%                             corrections and stellar aberration
%                             correction. Light time and stellar
%                             aberration corrections are applied as
%                             in the 'LT+S' case.
%
%      obsrvr   the name of the observing body, typically a spacecraft, the
%               earth, or a surface point on the earth.
%
%               [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               `obsrvr' is case-insensitive, and leading and trailing
%               blanks in `obsrvr' are not significant. Optionally, you may
%               supply a string containing the integer ID code for the
%               object. For example both 'EARTH' and '399' are legitimate
%               strings that indicate the earth is the observer.
%
%   the call:
%
%      [spoint] = cspice_subsol( method, target, et, abcorr, obsrvr )
%
%   returns:
%
%      spoint   the sub-solar point(s) on the target body at `et' expressed
%               relative to the body-fixed frame of the target body.
%
%               [3,n] = size(spoint); double = class(spoint)
%
%               The sub-solar point is defined either as the point on
%               the target body that is closest to the sun, or the
%               target surface intercept of the line containing the
%               target's center and the sun's center; the input
%               argument `method' selects the definition to be used.
%
%               The body-fixed frame, which is time-dependent, is
%               evaluated at `et' if `abcorr' is 'NONE'; otherwise the
%               frame is evaluated at et-lt, where `lt' is the one way
%               light time from target to observer.
%
%               The state of the target body is corrected for
%               aberration as specified by `abcorr'; the corrected
%               state is used in the geometric computation. As
%               indicated above, the rotation of the target is
%               retarded by one way light time if `abcorr' specifies
%               that light time correction is to be done.
%
%               The state of the sun as seen from the observing
%               body is also corrected for aberration as specified
%               by `abcorr'. The corrections, when selected, are
%               applied at the epoch et-lt, where `lt' is the one way
%               light time from target to observer.
%
%               `spoint' returns with the same vectorization measure, N,
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
%   1) Find the sub-solar position on the Earth as seen from the Moon at
%      at epoch JAN 1, 2006 using the 'near point' then the 'intercept'
%      options. Apply light time correction to return apparent position.
%
%      Compute the distance between the location of the sub-solar
%      points computed using the two different options, and the
%      angular separation between their position vectors.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: subsol_ex1.tm
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
%      function subsol_ex1()
%
%         %
%         % Load the meta kernel listing the needed SPK, PCK, LSK
%         % kernels.
%         %
%         cspice_furnsh( 'subsol_ex1.tm' )
%
%         et = cspice_str2et( 'JAN 1, 2006' );
%
%         %
%         % First use option 'Near Point'
%         %
%         point1 = cspice_subsol( 'near point', 'earth', et, ...
%                                 'lt+s', 'moon');
%
%         disp( 'Sub solar location coordinates - near point:' )
%         fprintf( '    %15.8f\n', point1 )
%
%         disp(' ')
%
%         %
%         % Now use option 'Intercept'
%         %
%         point2 = cspice_subsol( 'intercept', 'earth', et, ...
%                                 'lt+s', 'moon');
%
%         disp( 'Sub solar location coordinates - intercept' )
%         fprintf( '    %15.8f\n', point2 )
%
%         %
%         % Calculate the Euclidean distance between the two locations
%         % and the angular separation between the position vectors.
%         %
%         dist = norm( point1 - point2);
%         sep  = cspice_vsep(point1, point2 )*cspice_dpr;
%
%         disp(' ')
%
%         fprintf(['Distance between locations', ...
%                  '            (km): %8.5f\n'], dist);
%         fprintf(['Angular separation between', ...
%                  ' locations (deg): %8.5f\n'], sep );
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
%      Sub solar location coordinates - near point:
%           -5872.12721997
%             -91.14113098
%           -2479.72444945
%
%      Sub solar location coordinates - intercept
%           -5866.09273280
%             -91.04746986
%           -2493.87452061
%
%      Distance between locations            (km): 15.38338
%      Angular separation between locations (deg):  0.13826
%
%
%      Note that the difference between the location of the sub-solar
%      points computed using the two different options, results from the
%      non-spherical shape of the Earth.
%
%-Particulars
%
%   cspice_subsol computes the sub-solar point on a target body, as seen by
%   a specified observer.
%
%   There are two different popular ways to define the sub-solar point:
%   "nearest point on target to the sun" or "target surface intercept of
%   line containing target and sun." These coincide when the target is
%   spherical and generally are distinct otherwise.
%
%   When comparing sub-point computations with results from sources
%   other than SPICE, it's essential to make sure the same geometric
%   definitions are used.
%
%-Exceptions
%
%   If any of the listed errors occur, the output arguments are
%   left unchanged.
%
%   1)  If the input argument `method' is not recognized, the error
%       SPICE(DUBIOUSMETHOD) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If either of the input body names `target' or `obsrvr' cannot be
%       mapped to NAIF integer codes, the error SPICE(IDCODENOTFOUND)
%       is signaled by a routine in the call tree of this routine.
%
%   3)  If `obsrvr' and `target' map to the same NAIF integer ID codes,
%       the error SPICE(BODIESNOTDISTINCT) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If frame definition data enabling the evaluation of the state
%       of the target relative to the observer in target body-fixed
%       coordinates have not been loaded prior to calling cspice_subsol, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   5)  If the specified aberration correction is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If insufficient ephemeris data have been loaded prior to
%       calling cspice_subsol, an error is signaled by a
%       routine in the call tree of this routine.
%
%   7)  If the triaxial radii of the target body have not been loaded
%       into the kernel pool prior to calling cspice_subsol, an error is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If the size of the `target' body radii kernel variable is not
%       three, an error is signaled by a routine in the call tree of
%       this routine.
%
%   9)  If any of the three `target' body radii is less-than or equal to
%       zero, an error is signaled by a routine in the call tree of
%       this routine.
%
%   10) If PCK data supplying a rotation model for the target body
%       have not been loaded prior to calling cspice_subsol, an error is
%       signaled by a routine in the call tree of this routine.
%
%   11) If any of the input arguments, `method', `target', `et',
%       `abcorr' or `obsrvr', is undefined, an error is signaled by
%       the Matlab error handling system.
%
%   12) If any of the input arguments, `method', `target', `et',
%       `abcorr' or `obsrvr', is not of the expected type, or it does
%       not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPK, PCK, and frame data must be available to
%   the calling program before this routine is called. Typically
%   the data are made available by loading kernels; however the
%   data may be supplied via routine interfaces if applicable.
%
%   The following data are required:
%
%   -  SPK data: ephemeris data for sun, target, and observer must
%      be loaded. If aberration corrections are used, the states of
%      sun, target, and observer relative to the solar system
%      barycenter must be calculable from the available ephemeris
%      data. Ephemeris data are made available by loading
%      one or more SPK files via cspice_furnsh.
%
%   -  PCK data: triaxial radii for the target body must be loaded
%      into the kernel pool. Typically this is done by loading a
%      text PCK file via cspice_furnsh.
%
%   -  Further PCK data:  a rotation model for the target body must
%      be loaded. This may be provided in a text or binary PCK
%      file which is loaded via cspice_furnsh.
%
%   -  Frame data: if a frame definition is required to convert
%      the sun, observer, and target states to the body-fixed frame
%      of the target, that definition must be available in the
%      kernel pool. Typically the definition is supplied by loading
%      a frame kernel via cspice_furnsh.
%
%   In all cases, kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   1)  The appropriate kernel data must have been loaded before this
%       routine is called. See the -Files section above.
%
%-Required_Reading
%
%   MICE.REQ
%   FRAMES.REQ
%   PCK.REQ
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
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.4, 30-OCT-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.3, 23-JUN-2014 (NJB)
%
%       Updated description of converged Newtonian light time
%       correction.
%
%   -Mice Version 1.0.2, 18-MAY-2010 (BVS)
%
%       Index line now states that this routine is deprecated.
%
%   -Mice Version 1.0.1, 11-NOV-2008 (EDW)
%
%       Edits to header; -Abstract now states that this routine is
%       deprecated.
%
%   -Mice Version 1.0.0, 07-MAR-2007 (EDW)
%
%-Index_Entries
%
%   DEPRECATED sub-solar point
%
%-&

function [spoint] = cspice_subsol( method, target, et, abcorr, obsrvr )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise

         error ( ['Usage: [_spoint(3)_ ] = '            ...
                  'cspice_subsol( `method`, `target`, ' ...
                                  '_et_, `abcorr`, `obsrvr`)']  )

   end

   %
   % Call the MEX library.
   %
   try
      [spoint] = mice('subsol_c', method, target, et, abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end


