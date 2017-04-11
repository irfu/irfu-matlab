%-Abstract
%
%   CSPICE_SUBPNT computes the rectangular coordinates of the
%   sub-observer point on a target body at a specified epoch,
%   optionally corrected for light time and stellar aberration.
%
%   The surface of the target body may be represented by a triaxial
%   ellipsoid or by topographic data provided by DSK files.
%
%   This routine supersedes cspice_subpt, which does not have an input
%   argument for the target body-fixed frame name.
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
%      method   a short string providing parameters defining
%               the computation method to be used. In the syntax
%               descriptions below, items delimited by brackets
%               are optional.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               `method' may be assigned the following values:
%
%                  'NEAR POINT/ELLIPSOID'
%
%                     The sub-observer point computation uses a
%                     triaxial ellipsoid to model the surface of the
%                     target body. The sub-observer point is defined
%                     as the nearest point on the target relative to
%                     the observer.
%
%                     The word 'NADIR' may be substituted for the phrase
%                     'NEAR POINT' in the string above.
%
%                     For backwards compatibility, the older syntax
%
%                        'Near point: ellipsoid'
%
%                     is accepted as well.
%
%
%                  'INTERCEPT/ELLIPSOID'
%
%                     The sub-observer point computation uses a
%                     triaxial ellipsoid to model the surface of the
%                     target body. The sub-observer point is defined
%                     as the target surface intercept of the line
%                     containing the observer and the target's
%                     center.
%
%                     For backwards compatibility, the older syntax
%
%                        'Intercept: ellipsoid'
%
%                     is accepted as well.
%
%
%                  'NADIR/DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                     The sub-observer point computation uses DSK data
%                     to model the surface of the target body. The
%                     sub-observer point is defined as the intercept, on
%                     the surface represented by the DSK data, of the
%                     line containing the observer and the nearest point
%                     on the target's reference ellipsoid. If multiple
%                     such intercepts exist, the one closest to the
%                     observer is selected.
%
%                     Note that this definition of the sub-observer
%                     point is not equivalent to the "nearest point on
%                     the surface to the observer." The phrase 'NEAR
%                     POINT' may NOT be substituted for 'NADIR' in the
%                     string above.
%
%                     The surface list specification is optional. The
%                     syntax of the list is
%
%                        <surface 1> [, <surface 2>...]
%
%                     If present, it indicates that data only for the
%                     listed surfaces are to be used; however, data
%                     need not be available for all surfaces in the
%                     list. If absent, loaded DSK data for any surface
%                     associated with the target body are used.
%
%                     The surface list may contain surface names or
%                     surface ID codes. Names containing blanks must
%                     be delimited by double quotes, for example
%
%                        'SURFACES = "Mars MEGDR 128 PIXEL/DEG"'
%
%                     If multiple surfaces are specified, their names
%                     or IDs must be separated by commas.
%
%                     See the Particulars section below for details
%                     concerning use of DSK data.
%
%
%                  'INTERCEPT/DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                     The sub-observer point computation uses DSK data
%                     to model the surface of the target body. The
%                     sub-observer point is defined as the target
%                     surface intercept of the line containing the
%                     observer and the target's center.
%
%                     If multiple such intercepts exist, the one closest
%                     to the observer is selected.
%
%                     The surface list specification is optional. The
%                     syntax of the list is identical to that for the
%                     NADIR option described above.
%
%
%                  Neither case nor white space are significant in
%                  `method', except within double-quoted strings. For
%                  example, the string ' eLLipsoid/nearpoint ' is valid.
%
%                  Within double-quoted strings, blank characters are
%                  significant, but multiple consecutive blanks are
%                  considered equivalent to a single blank. Case is
%                  not significant. So
%
%                     "Mars MEGDR 128 PIXEL/DEG"
%
%                  is equivalent to
%
%                     " mars megdr  128  pixel/deg "
%
%                  but not to
%
%                     "MARS MEGDR128PIXEL/DEG"
%
%      target   the name of the target body. The target
%               body is an ephemeris object (its trajectory is given by
%               SPK data), and is an extended object.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               The string 'target' is case-insensitive, and leading
%               and trailing blanks in 'target' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object. For example both
%               'MOON' and '301' are legitimate strings that indicate
%               the moon is the target body.
%
%               When the target body's surface is represented by a
%               tri-axial ellipsoid, this routine assumes that a
%               kernel variable representing the ellipsoid's radii is
%               present in the kernel pool. Normally the kernel
%               variable would be defined by loading a PCK file.
%
%      et       the epoch(s), expressed as seconds past J2000 TDB, of the
%               observer: 'et' is the epoch at which the observer's state
%               is computed.
%
%               [1,n] = size(et); double = class(et)
%
%               When aberration corrections are not used, 'et' is also
%               the epoch at which the position and orientation of
%               the target body are computed.
%
%               When aberration corrections are used, 'et' is the epoch
%               at which the observer's state relative to the solar
%               system barycenter is computed; in this case the
%               position and orientation of the target body are
%               computed at et-lt or et+lt, where 'lt' is the one-way
%               light time between the sub-observer point and the
%               observer, and the sign applied to 'lt' depends on the
%               selected correction. See the description of 'abcorr'
%               below for details.
%
%      fixref   the name of a body-fixed reference frame centered
%               on the target body. `fixref' may be any such frame
%               supported by the SPICE system, including built-in
%               frames (documented in the Frames Required Reading)
%               and frames defined by a loaded frame kernel (FK). The
%               string `fixref' is case-insensitive, and leading and
%               trailing blanks in `fixref' are not significant.
%
%               [1,c3] = size(fixref); char = class(fixref)
%
%                  or
%
%               [1,1] = size(fixref); cell = class(fixref)
%
%               The output sub-observer point `spoint' and the
%               observer-to-sub-observer point vector `srfvec' will be
%               expressed relative to this reference frame.
%
%      abcorr   the aberration correction to apply when computing the
%               observer-target state and the orientation of the target body.
%
%               [1,c4] = size(abcorr); char = class(abcorr)
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
%               are described below. 'abcorr' may be any of the
%               following:
%
%                     'NONE'     Apply no correction. Return the
%                                geometric sub-observer point on the
%                                target body.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the sub-observer point (note: NOT
%               between the observer and the target body's center).
%               The following values of 'abcorr' apply to the
%               "reception" case in which photons depart from the
%               sub-observer point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at 'et':
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                yields the location of sub-observer
%                                point at the moment it emitted photons
%                                arriving at the observer at 'et'.
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
%               The following values of 'abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at 'et' and arrive at the
%               sub-observer point at the light-time corrected epoch
%               et+lt:
%
%                     'XLT'      "Transmission" case: correct for
%                                one-way light time using a Newtonian
%                                formulation. This correction yields the
%                                sub-observer location at the moment it
%                                receives photons emitted from the
%                                observer's location at 'et'.
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
%      obsrvr   the name of the observing body. The
%               observing body is an ephemeris object: it typically
%               is a spacecraft, the earth, or a surface point on the
%               earth. 'obsrvr' is case-insensitive, and leading and
%               'obsrvr' are not significant. Optionally, you may
%               trailing blanks in supply a string containing the integer
%               ID code for the object. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               observer.
%
%               [1,c5] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%   the call:
%
%      [spoint, trgepc, srfvec] = cspice_subpnt( method, target, ...
%                                                et,     fixref, ...
%                                                abcorr, obsrvr)
%
%   returns:
%
%      spoint   the array(s) defining the sub-observer point
%               on the target body.
%
%               [3,n] = size(spoint); double = class(spoint)
%
%               For target shapes modeled by ellipsoids, the
%               sub-observer point is defined either as the point on
%               the target body that is closest to the observer, or
%               the target surface intercept of the line from the
%               observer to the target's center.
%
%               For target shapes modeled by topographic data
%               provided by DSK files, the sub-observer point is
%               defined as the target surface intercept of the line
%               from the observer to either the nearest point on the
%               reference ellipsoid, or to the target's center. If
%               multiple such intercepts exist, the one closest to
%               the observer is selected.
%
%               The input argument `method' selects the target shape
%               model and sub-observer point definition to be used.
%
%               `spoint' is expressed in Cartesian coordinates,
%               relative to the body-fixed target frame designated by
%               `fixref'. The body-fixed target frame is evaluated at
%               the sub-observer epoch `trgepc' (see description below).
%
%               When light time correction is used, the duration of
%               light travel between `spoint' to the observer is
%               considered to be the one way light time.
%
%               When aberration corrections are used, `spoint' is
%               computed using target body position and orientation
%               that have been adjusted for the corrections
%               applicable to `spoint' itself rather than to the target
%               body's center. In particular, if the stellar
%               aberration correction applicable to `spoint' is
%               represented by a shift vector S, then the light-time
%               corrected position of the target is shifted by S
%               before the sub-observer point is computed.
%
%               The components of `spoint' have units of km.
%
%      trgepc   the "sub-observer point epoch(s)." 'trgepc' is
%               defined as follows: letting 'lt' be the one-way
%               light time between the observer and the sub-observer point,
%               'trgepc' is the epoch et-lt, et+lt, or 'et' depending on
%               whether the requested aberration correction is,
%               respectively, for received radiation, transmitted
%               radiation, or omitted. 'lt' is computed using the
%               method indicated by 'abcorr'.
%
%               [1,n] = size(trgepc); double = class(trgepc)
%
%              'trgepc' is expressed as seconds past J2000 TDB.
%
%      srfvec   the array(s) defining the position vector from
%               the observer at 'et' to 'spoint'. 'srfvec'
%               is expressed in the target body-fixed  reference frame
%               designated by 'fixref', evaluated at  'trgepc'.
%
%               [3,n] = size(spoint); double = class(spoint)
%
%               The components of 'srfvec' are given in units of km.
%
%               One can use the CSPICE function vnorm_c to obtain the
%               distance between the observer and 'spoint':
%
%                  dist = norm( srfvec )
%
%               The observer's position 'obspos', relative to the
%               target body's center, where the center's position is
%               corrected for aberration effects as indicated by
%               'abcorr', can be computed with:
%
%                  obspos = spoint - srfvec
%
%               To transform the vector 'srfvec' from a reference frame
%               'fixref' at time 'trgepc' to a time-dependent reference
%               frame 'ref' at time 'et', the routine 'cspice_pxfrm2' should be
%               called. Let 'xform' be the 3x3 matrix representing the
%               rotation from the reference frame 'fixref' at time
%               'trgepc' to the reference frame 'ref' at time 'et'. Then
%               'srfvec' can be transformed to the result 'refvec' as
%               follows:
%
%                  xform  = cspice_pxfrm2 ( fixref, ref, trgepc, et )
%                  refvec = xform * srfvec
%
%               'spoint', 'trgepc', and 'srfvec' return with the same
%               vectorization measure, N, as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Find the sub-Earth point on Mars for a specified time. Perform
%      the computation twice, using both the "intercept" and "near
%      point" options. Display the location of both the Earth and the
%      sub-Earth point using both planetocentric and planetographic
%      coordinates.
%
%      %
%      % Load kernel files via the meta-kernel.
%      %
%      cspice_furnsh( '/kernels/standard.tm' );
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et0 = cspice_str2et( '2008 aug 11 00:00:00' );
%
%      %
%      % Create a vector of times. The code will also run for 'et'
%      % a scalar.
%      %
%      et = [0:10]*cspice_spd + et0;
%
%      %
%      % Look up the target body's radii. We'll use these to
%      % convert Cartesian to planetodetic coordinates. Use
%      % the radii to compute the flattening coefficient of
%      % the reference ellipsoid.
%      %
%      radii  = cspice_bodvrd( 'MARS', 'RADII', 3 );
%
%      %
%      % Let RE and RP be, respectively, the equatorial and
%      % polar radii of the target.
%      %
%      re = radii(1);
%      rp = radii(3);
%      f = ( re-rp)/re;
%
%      %
%      % Compute sub-observer point using light time and stellar
%      % aberration corrections. Use the "target surface intercept"
%      % definition of the sub-observer point on the first loop
%      % iteration, and use the "near point" definition on the
%      % second.
%      %
%
%      method = { 'Intercept:  ellipsoid', 'Near point: ellipsoid' };
%
%      for i=1:2
%
%         [spoint, trgepc, srfvec] = cspice_subpnt( method(i), ...
%                         'MARS', et, 'IAU_MARS', 'LT+S', 'EARTH' );
%
%         N = size(spoint, 2);
%
%         %
%         % Compute the observer's distance from SPOINT.
%         %
%         odist = norm(srfvec);
%
%         %
%         % Convert the sub-observer point's rectangular coordinates
%         % to planetographic longitude, latitude and altitude.
%         % Convert radians to degrees.
%         %
%         [ spglon, spglat, spgalt] = cspice_recpgr( 'mars', spoint, re, f );
%
%         spglon = spglon * cspice_dpr;
%         spglat = spglat * cspice_dpr;
%
%         %
%         % Convert sub-observer point's rectangular coordinates to
%         % planetocentric radius, longitude, and latitude. Convert
%         % radians to degrees.
%         %
%         [ spcrad, spclon, spclat ] =cspice_reclat( spoint ) ;
%
%         spclon = spclon * cspice_dpr;
%         spclat = spclat * cspice_dpr;
%
%         %
%         % Compute the observer's position relative to the center of the
%         % target, where the center's location has been adjusted using
%         % the aberration corrections applicable to the sub-point.
%         % Express the observer's location in geodetic coordinates.
%         %
%         obspos = spoint - srfvec;
%
%         [ opglon, opglat, opgalt] = cspice_recpgr( 'mars', obspos, re, f );
%
%         opglon = opglon * cspice_dpr;
%         opglat = opglat * cspice_dpr;
%
%         %
%         % Convert the observer's rectangular coordinates to planetodetic
%         % longitude, latitude and altitude. Convert radians to degrees.
%         %
%         [opcrad, opclon, opclat] = cspice_reclat( obspos ) ;
%
%         opclon = opclon * cspice_dpr;
%         opclat = opclat * cspice_dpr;
%
%         utcstr = cspice_et2utc( et, 'C', 6);
%
%         for j=1:N
%
%           fprintf( 'Computational Method %s\n\n', char(method(i)) )
%
%           fprintf( 'Time (UTC):                          %s\n',  ...
%                                                        utcstr(j,:) )
%
%           fprintf(                                                  ...
%           'Observer altitude                      (km) = %21.9f\n', ...
%                                                        opgalt(j) )
%
%           fprintf(                                                  ...
%           'Length of SRFVEC                       (km) = %21.9f\n', ...
%                                               norm(srfvec(:,j))  )
%
%           fprintf(                                                  ...
%           'Sub-observer point altitude            (km) = %21.9f\n', ...
%                                                        spgalt(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetographic longitude (deg) = %21.9f\n', ...
%                                                        spglon(j) )
%
%           fprintf(                                                  ...
%           'Observer planetographic longitude     (deg) = %21.9f\n', ...
%                                                        opglon(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetographic latitude  (deg) = %21.9f\n', ...
%                                                        spglat(j) )
%
%           fprintf(                                                  ...
%           'Observer planetographic latitude      (deg) = %21.9f\n', ...
%                                                        opglat(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetocentric longitude (deg) = %21.9f\n', ...
%                                                        spclon(j) )
%
%           fprintf(                                                  ...
%           'Observer planetocentric longitude     (deg) = %21.9f\n', ...
%                                                        opclon(j) )
%
%           fprintf(                                                  ...
%           'Sub-observer planetocentric latitude  (deg) = %21.9f\n', ...
%                                                        spclat(j) )
%
%           fprintf(                                                  ...
%           'Observer planetocentric latitude      (deg) = %21.9f\n', ...
%                                                        opclat(j) )
%
%           fprintf( '\n')
%
%         end
%
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Computational Method Intercept:  ellipsoid
%
%      Observer altitude                      (km) =   349199089.604656994
%      Length of SRFVEC                       (km) =   349199089.641352594
%      Sub-observer point altitude            (km) =          -0.000000000
%      Sub-observer planetographic longitude (deg) =         199.302304818
%      Observer planetographic longitude     (deg) =         199.302304818
%      Sub-observer planetographic latitude  (deg) =          26.262401078
%      Observer planetographic latitude      (deg) =          25.994936593
%      Sub-observer planetocentric longitude (deg) =         160.697695182
%      Observer planetocentric longitude     (deg) =         160.697695182
%      Sub-observer planetocentric latitude  (deg) =          25.994934013
%      Observer planetocentric latitude      (deg) =          25.994934013
%
%      Computational Method Near point: ellipsoid
%
%      Observer altitude                      (km) =   349199089.604648590
%      Length of SRFVEC                       (km) =   349199089.604648590
%      Sub-observer point altitude            (km) =          -0.000000000
%      Sub-observer planetographic longitude (deg) =         199.302304819
%      Observer planetographic longitude     (deg) =         199.302304819
%      Sub-observer planetographic latitude  (deg) =          25.994936593
%      Observer planetographic latitude      (deg) =          25.994936593
%      Sub-observer planetocentric longitude (deg) =         160.697695181
%      Observer planetocentric longitude     (deg) =         160.697695181
%      Sub-observer planetocentric latitude  (deg) =          25.729407071
%      Observer planetocentric latitude      (deg) =          25.994934013
%
%-Particulars
%
%   A sister version of this routine exists named mice_subpnt that returns
%   the output arguments as fields in a single structure.
%
%   For ellipsoidal target bodies, there are two different popular
%   ways to define the sub-observer point: "nearest point on the
%   target to the observer" or "target surface intercept of the line
%   containing observer and target." These coincide when the target
%   is spherical and generally are distinct otherwise.
%
%   For target body shapes modeled using topographic data provided by
%   DSK files, the "surface intercept" notion is valid, but the
%   "nearest point on the surface" computation is both inefficient to
%   execute and may fail to yield a result that is "under" the
%   observer in an intuitively clear way. The NADIR option for DSK
%   shapes instead finds the surface intercept of a ray that passes
%   through the nearest point on the target reference ellipsoid. For
%   shapes modeled using topography, there may be multiple
%   ray-surface intercepts; the closest one to the observer is
%   selected.
%
%   The NADIR definition makes sense only if the target shape is
%   reasonably close to the target's reference ellipsoid. If the
%   target is very different---the nucleus of comet
%   Churyumov-Gerasimenko is an example---the intercept definition
%   should be used.
%
%   This routine computes light time corrections using light time
%   between the observer and the sub-observer point, as opposed to
%   the center of the target. Similarly, stellar aberration
%   corrections done by this routine are based on the direction of
%   the vector from the observer to the light-time corrected
%   sub-observer point, not to the target center. This technique
%   avoids errors due to the differential between aberration
%   corrections across the target body. Therefore it's valid to use
%   aberration corrections with this routine even when the observer
%   is very close to the sub-observer point, in particular when the
%   observer to sub-observer point distance is much less than the
%   observer to target center distance.
%
%   When comparing sub-observer point computations with results from
%   sources other than SPICE, it's essential to make sure the same
%   geometric definitions are used.
%
%
%   Using DSK data
%   ==============
%
%      DSK loading and unloading
%      -------------------------
%
%      DSK files providing data used by this routine are loaded by
%      calling cspice_furnsh and can be unloaded by calling cspice_unload or
%      cspice_kclear. See the documentation of cspice_furnsh for limits on
%      numbers of loaded DSK files.
%
%      For run-time efficiency, it's desirable to avoid frequent
%      loading and unloading of DSK files. When there is a reason to
%      use multiple versions of data for a given target body---for
%      example, if topographic data at varying resolutions are to be
%      used---the surface list can be used to select DSK data to be
%      used for a given computation. It is not necessary to unload
%      the data that are not to be used. This recommendation presumes
%      that DSKs containing different versions of surface data for a
%      given body have different surface ID codes.
%
%
%      DSK data priority
%      -----------------
%
%      A DSK coverage overlap occurs when two segments in loaded DSK
%      files cover part or all of the same domain---for example, a
%      given longitude-latitude rectangle---and when the time
%      intervals of the segments overlap as well.
%
%      When DSK data selection is prioritized, in case of a coverage
%      overlap, if the two competing segments are in different DSK
%      files, the segment in the DSK file loaded last takes
%      precedence. If the two segments are in the same file, the
%      segment located closer to the end of the file takes
%      precedence.
%
%      When DSK data selection is unprioritized, data from competing
%      segments are combined. For example, if two competing segments
%      both represent a surface as sets of triangular plates, the
%      union of those sets of plates is considered to represent the
%      surface.
%
%      Currently only unprioritized data selection is supported.
%      Because prioritized data selection may be the default behavior
%      in a later version of the routine, the UNPRIORITIZED keyword is
%      required in the `method' argument.
%
%
%      Syntax of the `method' input argument
%      -----------------------------------
%
%      The keywords and surface list in the `method' argument
%      are called "clauses." The clauses may appear in any
%      order, for example
%
%         'NADIR/DSK/UNPRIORITIZED/<surface list>'
%         'DSK/NADIR/<surface list>/UNPRIORITIZED'
%         'UNPRIORITIZED/<surface list>/DSK/NADIR'
%
%      The simplest form of the `method' argument specifying use of
%      DSK data is one that lacks a surface list, for example:
%
%         'NADIR/DSK/UNPRIORITIZED'
%         'INTERCEPT/DSK/UNPRIORITIZED'
%
%      For applications in which all loaded DSK data for the target
%      body are for a single surface, and there are no competing
%      segments, the above strings suffice. This is expected to be
%      the usual case.
%
%      When, for the specified target body, there are loaded DSK
%      files providing data for multiple surfaces for that body, the
%      surfaces to be used by this routine for a given call must be
%      specified in a surface list, unless data from all of the
%      surfaces are to be used together.
%
%      The surface list consists of the string
%
%         'SURFACES = '
%
%      followed by a comma-separated list of one or more surface
%      identifiers. The identifiers may be names or integer codes in
%      string format. For example, suppose we have the surface
%      names and corresponding ID codes shown below:
%
%         Surface Name                              ID code
%         ------------                              -------
%         "Mars MEGDR 128 PIXEL/DEG"                1
%         "Mars MEGDR 64 PIXEL/DEG"                 2
%         "Mars_MRO_HIRISE"                         3
%
%      If data for all of the above surfaces are loaded, then
%      data for surface 1 can be specified by either
%
%         'SURFACES = 1'
%
%      or
%
%         'SURFACES = "Mars MEGDR 128 PIXEL/DEG"'
%
%      Double quotes are used to delimit the surface name
%      because it contains blank characters.
%
%      To use data for surfaces 2 and 3 together, any
%      of the following surface lists could be used:
%
%         'SURFACES = 2, 3'
%
%         'SURFACES = "Mars MEGDR  64 PIXEL/DEG", 3'
%
%         'SURFACES = 2, Mars_MRO_HIRISE'
%
%         'SURFACES = "Mars MEGDR 64 PIXEL/DEG", Mars_MRO_HIRISE'
%
%      An example of a `method' argument that could be constructed
%      using one of the surface lists above is
%
%      'NADIR/DSK/UNPRIORITIZED/SURFACES= "Mars MEGDR 64 PIXEL/DEG",3'
%
%
%      Aberration corrections
%      ----------------------
%
%      For irregularly shaped target bodies, the distance between the
%      observer and the nearest surface intercept need not be a
%      continuous function of time; hence the one-way light time
%      between the intercept and the observer may be discontinuous as
%      well. In such cases, the computed light time, which is found
%      using iterative algorithm, may converge slowly or not at all.
%      In all cases, the light time computation will terminate, but
%      the result may be less accurate than expected.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine subpnt_c.
%
%   MICE.REQ
%   DSK.REQ
%   FRAMES.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 2.0.0, 04-APR-2017, EDW (JPL), NJB (JPL)
%
%       Header update to reflect support for use of DSKs.
%
%       Vectorized interface on input 'et'.
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.2, 25-OCT-2011, SCK (JPL)
%
%       References to the new 'cspice_pxfrm2' routine were
%       added to the 'I/O returns' section.
%
%   -Mice Version 1.0.1, 12-MAY-2009, EDW (JPL)
%
%       Corrected type in I/O call description. The call description
%       lacked the 'fixref' argument.
%
%   -Mice Version 1.0.0, 30-JAN-2008, EDW (JPL)
%
%-Index_Entries
%
%   find sub-observer point on target body 
%   find sub-spacecraft point on target body 
%   find nearest point to observer on target body 
%
%-&

function [spoint, trgepc, srfvec] = cspice_subpnt( method, target, et, ...
                                                   fixref, abcorr, obsrvr )

   switch nargin
      case 6

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);

      otherwise

         error ( ['Usage: [_spoint_, _trgepc_, _srfvec_] = ' ...
                  'cspice_subpnt( `method`, `target`,'       ...
                  ' _et_, `fixref`, `abcorr`, `obsrvr`)']  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [subpnt] = mice('subpnt_s', method, target, et, fixref, abcorr, obsrvr);
      spoint   = reshape( [subpnt.spoint], 3, [] );
      trgepc   = reshape( [subpnt.trgepc], 1, [] );
      srfvec   = reshape( [subpnt.srfvec], 3, [] );
   catch
      rethrow(lasterror)
   end



