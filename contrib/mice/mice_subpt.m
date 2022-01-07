%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   mice_subpnt. This routine is supported for purposes of
%   backward compatibility only.
%
%   MICE_SUBPT determines the coordinates of the sub-observer point
%   on a target body at a particular epoch, optionally corrected
%   for planetary (light time) and stellar aberration. The call also
%   returns the observer's altitude above the target body.
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
%      method   a string providing parameters defining the
%               computation method to use.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               The choices are:
%
%                  'Near point'       The sub-observer point is
%                                     defined as the nearest point on
%                                     the target relative to the
%                                     observer.
%
%                  'Intercept'        The sub-observer point is
%                                     defined as the target surface
%                                     intercept of the line
%                                     containing the observer and the
%                                     target's center.
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
%      target   the name of the observed target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               `target' is case-insensitive, and leading and trailing blanks
%               in `target' are not significant. Optionally, you may supply
%               a string containing the integer ID code for the object.
%               For example both 'MOON' and '301' are legitimate strings
%               that indicate the moon is the target body.
%
%               This routine assumes that the target body is modeled by
%               a tri-axial ellipsoid, and that a PCK file containing
%               its radii has been loaded into the kernel pool via
%               cspice_furnsh.
%
%      et       the epoch(s), expressed as seconds past J2000 TDB, of the
%               observer: `et' is the epoch at which the observer's state
%               is computed.
%
%               [1,n] = size(et); double = class(et)
%
%               When aberration corrections are not used, `et' is also
%               the epoch at which the position and orientation of
%               the target body are computed.
%
%               When aberration corrections are used, `et' is the epoch
%               at which the observer's state relative to the solar
%               system barycenter is computed; in this case the
%               position and orientation of the target body are
%               computed at et-lt or et+lt, where `lt' is the one-way
%               light time between the sub-observer point and the
%               observer, and the sign applied to `lt' depends on the
%               selected correction. See the description of `abcorr'
%               below for details.
%
%      abcorr   the aberration correction to apply
%               when computing the observer-target state and the
%               orientation of the target body.
%
%               [1,c3] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               For remote sensing applications, where the apparent
%               sub-observer point seen by the observer is desired,
%               normally either of the corrections
%
%                     'LT+S'
%                     'CN+S'
%
%               should be used. These and the other supported options
%               are described below. `abcorr' may be any of the
%               following:
%
%                     'NONE'     Apply no correction. Return the
%                                geometric sub-observer point on the
%                                target body.
%
%               Let `lt' represent the one-way light time between the
%               observer and the sub-observer point (note: NOT
%               between the observer and the target body's center).
%               The following values of `abcorr' apply to the
%               "reception" case in which photons depart from the
%               sub-observer point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at `et':
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                yields the location of sub-observer
%                                point at the moment it emitted photons
%                                arriving at the observer at `et'.
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                Both the target position as seen by the
%                                observer, and rotation of the target
%                                body, are corrected for light time.
%
%                     'LT+S'     Correct for one-way light time and
%                                stellar aberration using a Newtonian
%                                formulation. This option modifies the
%                                state obtained with the 'LT' option to
%                                account for the observer's velocity
%                                relative to the solar system
%                                barycenter. The result is the apparent
%                                sub-observer point as seen by the
%                                observer.
%
%                     'CN'       Converged Newtonian light time
%                                correction. In solving the light time
%                                equation, the 'CN' correction iterates
%                                until the solution converges. Both the
%                                position and rotation of the target
%                                body are corrected for light time.
%
%                     'CN+S'     Converged Newtonian light time and
%                                stellar aberration corrections. This
%                                option produces a solution that is at
%                                least as accurate at that obtainable
%                                with the 'LT+S' option. Whether the 'CN+S'
%                                solution is substantially more accurate
%                                depends on the geometry of the
%                                participating objects and on the
%                                accuracy of the input data. In all
%                                cases this routine will execute more
%                                slowly when a converged solution is
%                                computed.
%
%               The following values of `abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at `et' and arrive at the
%               sub-observer point at the light-time corrected epoch
%               et+lt:
%
%                     'XLT'      "Transmission" case: correct for
%                                one-way light time using a Newtonian
%                                formulation. This correction yields the
%                                sub-observer location at the moment it
%                                receives photons emitted from the
%                                observer's location at `et'.
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                Both the target position as seen by the
%                                observer, and rotation of the target
%                                body, are corrected for light time.
%
%                     'XLT+S'    "Transmission" case: correct for
%                                one-way light time and stellar
%                                aberration using a Newtonian
%                                formulation  This option modifies the
%                                sub-observer point obtained with the
%                                'XLT' option to account for the
%                                observer's velocity relative to the
%                                solar system barycenter.
%
%                     'XCN'      Converged Newtonian light time
%                                correction. This is the same as XLT
%                                correction but with further iterations
%                                to a converged Newtonian light time
%                                solution.
%
%                     'XCN+S'    "Transmission" case: converged
%                                Newtonian light time and stellar
%                                aberration corrections.
%
%      obsrvr   the scalar string name of the observing body.
%
%               [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               The observing body is an ephemeris object: it typically
%               is a spacecraft, the earth, or a surface point on the
%               earth. `obsrvr' is case-insensitive, and leading and
%               `obsrvr' are not significant. Optionally, you may
%               trailing blanks in supply a string containing the integer
%               ID code for the object. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               observer.
%
%   the call:
%
%      [spoint] = mice_subpt( method, target, et, abcorr, obsrvr )
%
%   returns:
%
%      spoint   the structure(s) containing the results of the calculation.
%
%               [1,n] = size(spoint); struct = class(spoint)
%
%               Each structure consists of the fields:
%
%                  pos      the array(s) defining the sub-observer point
%                           on the target body.
%
%                           [3,1]  = size(spoint(i).pos);
%                           double = class(spoint(i).pos(i))
%
%                           The sub-observer point is defined either as the
%                           point on the target body that is closest to the
%                           observer, or the target surface intercept of the
%                           line from the observer to the target's center;
%                           the input argument `method' selects the
%                           definition to be used.
%
%                           The body-fixed frame, which is time-dependent, is
%                           evaluated at `et' if `abcorr' is 'NONE'; otherwise
%                           the frame is evaluated at et-lt, where `lt' is the
%                           one-way light time from target to observer.
%
%                           The state of the target body is corrected for
%                           aberration as specified by `abcorr'; the corrected
%                           state is used in the geometric computation. As
%                           indicated above, the rotation of the target is
%                           retarded by one-way light time if `abcorr'
%                           specifies that light time correction is to be
%                           done.
%
%                  alt      the values(s) of the altitude(s) of `obsrvr' above
%                           `target'.
%
%                           [1,1]  = size(spoint(i).alt);
%                           double = class(spoint(i).alt)
%
%                           When `method' specifies a "near point"
%                           computation, `alt' is truly altitude in the
%                           standard geometric sense: the length of a segment
%                           dropped from the observer to the target's surface,
%                           such that the segment is perpendicular to the
%                           surface at the contact point `spoint'.
%
%                           When `method' specifies an "intercept"
%                           computation, `alt' is still the length of the
%                           segment from the observer to the surface point
%                           `spoint', but this segment in general is not
%                           perpendicular to the surface.
%
%               `spoint' returns with the same vectorization measure, N, as
%               `et'.
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
%   1) Find the sub-point position of the Moon on the Earth at
%      epoch JAN 1, 2006 using the "near point" then the "intercept"
%      options. Apply light time correction to return apparent position.
%
%      Compute the distance between the location of the sub-points
%      computed using the two different options, and the angular
%      separation between their position vectors.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: subpt_ex1.tm
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
%      function subpt_ex1()
%
%         %
%         % Find the sub point position of the moon on the earth at
%         % a given time using the 'near point' then the 'intercept'
%         % options.
%         %
%         % Load the meta kernel listing the needed SPK, PCK, LSK
%         % kernels.
%         %
%         cspice_furnsh( 'subpt_ex1.tm' )
%
%         %
%         % Calculate the location of the sub point of the moon as
%         % seen from Earth at epoch JAN 1, 2006. Apply light time
%         % correction to return apparent position.
%         %
%         et = cspice_str2et( 'JAN 1, 2006' );
%
%         %
%         % First use option 'Near Point'
%         %
%         [point1] = mice_subpt( 'near point', 'earth', et, 'lt+s', 'moon');
%
%         disp( 'Sub-point location coordinates - near point:' )
%         fprintf( '    %15.8f  %15.8f  %15.8f\n', point1.pos )
%
%         disp( 'Sub-point observer altitude:' )
%         fprintf( '    %15.8f\n', point1.alt )
%
%         disp(' ')
%
%         %
%         % Now use option 'Intercept'
%         %
%         [point2] = mice_subpt( 'intercept', 'earth', et, 'lt+s', 'moon');
%
%
%         disp( 'Sub-point location coordinates - intercept:' )
%         fprintf( '    %15.8f  %15.8f  %15.8f\n', point2.pos )
%
%         disp( 'Sub-point observer altitude:' )
%         fprintf( '    %15.8f\n', point2.alt )
%
%         %
%         % Calculate the Euclidean distance between the two locations
%         % and the angular separation between the position vectors.
%         %
%         dist = norm( point1.pos - point2.pos);
%         sep  = cspice_vsep(point1.pos, point2.pos )*cspice_dpr;
%
%         disp(' ')
%         fprintf( 'Distance between locations            (km): %8.5f\n',  ...
%                                                                   dist);
%         fprintf( 'Angular separation between locations (deg): %8.5f\n',  ...
%                                                                   sep );
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
%      Sub-point location coordinates - near point:
%           -5532.84459341   -1443.48674452   -2816.23526877
%      Sub-point observer altitude:
%          356091.70835039
%
%      Sub-point location coordinates - intercept:
%           -5525.64307897   -1441.60791159   -2831.19586108
%      Sub-point observer altitude:
%          356091.73073431
%
%      Distance between locations            (km): 16.70961
%      Angular separation between locations (deg):  0.15020
%
%
%      Note that the difference between the location of the sub-points
%      computed using the two different options, results from the
%      non-spherical shape of the Earth.
%
%   2) Find the sub body position of the moon on the earth at
%      at epoch JAN 1, 2006 and for the next 3 months. Use the
%      'near point' option to calculate the physically
%      closest point between the two bodies.
%
%      Use the meta-kernel from the first example.
%
%
%      Example code begins here.
%
%
%      function subpt_ex2()
%
%         %
%         % Find the sub body position of the moon on the earth at
%         % at epoch JAN 1, 2006 and for the next 3 months. Use the
%         % 'near point' option to calculate the physically
%         % closest point between the two bodies.
%         %
%         % Load the meta kernel listing the needed SPK, PCK, LSK
%         % kernels.
%         %
%         cspice_furnsh( 'subpt_ex1.tm' )
%
%         %
%         % Convert the calendar string to ephemeris time.
%         %
%         et0 = cspice_str2et( 'JAN 1, 2006' );
%
%         %
%         % Fill an array with epochs, start with the epoch above then
%         % epochs in steps on one month ( thirty days in seconds)
%         %
%         et  = [0:3]*cspice_spd*30. + et0;
%
%         %
%         % Calculate the nearpoint of the moon with respect to earth at
%         % the epochs defined in 'et'.
%         %
%         [point] = mice_subpt( 'near point', 'earth', et, 'lt+s', 'moon');
%
%         %
%         % Convert the subpoint coordinates to lat/lon expressed in degrees
%         % with the radius.
%         %
%         % Extract from the `point' structure the 3XN array of position data.
%         %
%         position = reshape( [point(:).pos], 3, [] );
%
%         [radius, longitude, latitude] = cspice_reclat(position);
%         longitude                     = longitude * cspice_dpr;
%         latitude                      = latitude  * cspice_dpr;
%
%         %
%         % Convert the 'et' epochs to calendar format.
%         %
%         utc = cspice_et2utc( et, 'C', 3 );
%
%         for i=1:4
%            txt = sprintf( 'Moon subpoint epoch: %s', utc(i,:) );
%            disp( txt )
%
%            txt = sprintf( '   Longitude  (deg): %9.4f', longitude(i) );
%            disp( txt )
%
%            txt = sprintf( '   Latitude   (deg): %9.4f', latitude(i) );
%            disp( txt )
%            disp( ' ' )
%
%         end
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
%      Moon subpoint epoch: 2006 JAN 01 00:00:00.000
%         Longitude  (deg): -165.3778
%         Latitude   (deg):  -26.2210
%
%      Moon subpoint epoch: 2006 JAN 30 23:59:59.999
%         Longitude  (deg): -155.9765
%         Latitude   (deg):  -13.7804
%
%      Moon subpoint epoch: 2006 MAR 01 23:59:59.999
%         Longitude  (deg): -151.3550
%         Latitude   (deg):    3.9954
%
%      Moon subpoint epoch: 2006 MAR 31 23:59:59.998
%         Longitude  (deg): -146.5347
%         Latitude   (deg):   19.9326
%
%
%-Particulars
%
%   A sister version of this routine exists named cspice_subpt that returns
%   the structure field data as separate arguments.
%
%   Alternatively, if needed, the user can extract the field data from
%   vectorized `spoint' structures into separate arrays:
%
%      Extract the N `pos' field data to a 3XN array `position':
%
%         position = reshape( [spoint(:).pos], 3, [] )
%
%      Extract the N `alt' field data to a 1XN array `altitude':
%
%         altitude = reshape( [point(:).alt], 1, [] )
%
%   mice_subpt computes the sub-observer point on a target body.
%   (The sub-observer point is commonly called the sub-spacecraft
%   point when the observer is a spacecraft.) mice_subpt also
%   determines the altitude of the observer above the target body.
%
%   There are two different popular ways to define the sub-observer
%   point: "nearest point on target to observer" or "target surface
%   intercept of line containing observer and target." These
%   coincide when the target is spherical and generally are distinct
%   otherwise.
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
%       coordinates have not been loaded prior to calling cspice_subpt, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   5)  If the specified aberration correction is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If insufficient ephemeris data have been loaded prior to
%       calling cspice_subpt, an error is signaled by a
%       routine in the call tree of this routine.
%
%   7)  If the triaxial radii of the target body have not been loaded
%       into the kernel pool prior to calling cspice_subpt, an error is
%       signaled by a routine in the call tree of this routine.
%
%   8)  If any of the radii of the target body are non-positive, an
%       error is signaled by a routine in the call tree of this
%       routine. The target must be an extended body.
%
%   9)  If PCK data supplying a rotation model for the target body
%       have not been loaded prior to calling cspice_subpt, an error is
%       signaled by a routine in the call tree of this routine.
%
%   10) If any of the input arguments, `method', `target', `et',
%       `abcorr' or `obsrvr', is undefined, an error is signaled by
%       the Matlab error handling system.
%
%   11) If any of the input arguments, `method', `target', `et',
%       `abcorr' or `obsrvr', is not of the expected type, or it does
%       not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPK, PCK, and frame kernels must be loaded
%   prior by the calling program before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: ephemeris data for target and observer must be
%      loaded. If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%      Typically ephemeris data are made available by loading one
%      or more SPK files via cspice_furnsh.
%
%   -  PCK data: triaxial radii for the target body must be loaded
%      into the kernel pool. Typically this is done by loading a
%      text PCK file via cspice_furnsh.
%
%   -  Further PCK data: rotation data for the target body must
%      be loaded. These may be provided in a text or binary PCK
%      file. Either type of file may be loaded via cspice_furnsh.
%
%   -  Frame data: if a frame definition is required to convert
%      the observer and target states to the body-fixed frame of
%      the target, that definition must be available in the kernel
%      pool. Typically the definition is supplied by loading a
%      frame kernel via cspice_furnsh.
%
%   In all cases, kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   None.
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
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-OCT-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel. Reduced the size of the array of times
%       used to generate the output in example #2.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       extended the -Particulars section.
%
%       -Abstract and -Index_Entries now state that this routine is
%       deprecated.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 12-JAN-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.0, 16-DEC-2005 (EDW)
%
%-Index_Entries
%
%   DEPRECATED sub-observer point
%
%-&

function [spoint] = mice_subpt( method, target, et, abcorr, obsrvr )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise
         error ( ['Usage: [_spoint_] = mice_subpt( `method`, `target`,'    ...
                  ' _et_, `abcorr`, `obsrvr`)'] )
   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [spoint] = mice('subpt_s', method, target, et, abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end


