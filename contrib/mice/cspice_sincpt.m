%-Abstract
%
%   CSPICE_SINCPT computes the surface intercept of the ray on a target
%   body at a specified epoch, optionally corrected for light time and stellar
%   aberration, given an observer and a direction vector defining a ray,
%
%   This routine supersedes cspice_srfxpt, which does not have an input
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
%                  'ELLIPSOID'
%
%                     The intercept computation uses a triaxial
%                     ellipsoid to model the surface of the target
%                     body. The ellipsoid's radii must be available
%                     in the kernel pool.
%
%
%                  'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                     The intercept computation uses topographic data
%                     to model the surface of the target body. These
%                     data must be provided by loaded DSK files.
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
%               Neither case nor white space are significant in
%               `method', except within double-quoted strings. For
%               example, the string " eLLipsoid " is valid.
%
%               Within double-quoted strings, blank characters are
%               significant, but multiple consecutive blanks are
%               considered equivalent to a single blank. Case is
%               not significant. So
%
%                  "Mars MEGDR 128 PIXEL/DEG"
%
%               is equivalent to
%
%                  " mars megdr  128  pixel/deg "
%
%               but not to
%
%                  "MARS MEGDR128PIXEL/DEG"
%
%      target   the name of the target body. `target' is
%               case-insensitive, and leading and trailing blanks in
%               `target' are not significant. Optionally, you may
%               supply a string containing the integer ID code
%               for the object. For example both 'MOON' and '301'
%               are legitimate strings that indicate the moon is the
%               target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               When the target body's surface is represented by a
%               tri-axial ellipsoid, this routine assumes that a
%               kernel variable representing the ellipsoid's radii is
%               present in the kernel pool. Normally the kernel
%               variable would be defined by loading a PCK file.
%
%      et       the epoch, expressed as seconds past J2000 TDB, of the
%               observer: 'et' is the epoch at which the observer's state
%               is computed.
%
%               [1,1] = size(et); double = class(et)
%
%               When aberration corrections are not used, 'et' is also
%               the epoch at which the position and orientation of the
%               target body are computed.
%
%               When aberration corrections are used, 'et' is the epoch
%               at which the observer's state relative to the solar
%               system barycenter is computed; in this case the
%               position and orientation of the target body are
%               computed at et-lt or et+lt, where 'lt' is the one-way
%               light time between the intercept point and the
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
%               The output intercept point `spoint' and the observer-to-
%               intercept vector `srfvec' will be expressed relative to
%               this reference frame.
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
%               'abcorr' may be any of the following.
%
%                  "NONE"     Apply no correction. Return the
%                             geometric surface intercept point on the
%                             target body.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the surface intercept point (note: NOT
%               between the observer and the target body's center).
%               The following values of 'abcorr' apply to the
%               "reception" case in which photons depart from the
%               intercept point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at 'et':
%
%                  "LT"       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             yields the location of the surface
%                             intercept point at the moment it
%                             emitted photons arriving at the
%                             observer at 'et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             "LT" option uses one iteration.
%
%                             Both the target state as seen by the
%                             observer, and rotation of the target
%                             body, are corrected for light time.
%
%                  "LT+S"     Correct for one-way light time and
%                             stellar aberration using a Newtonian
%                             formulation. This option modifies the
%                             state obtained with the "LT" option to
%                             account for the observer's velocity
%                             relative to the solar system
%                             barycenter. The result is the apparent
%                             surface intercept point as seen by the
%                             observer.
%
%                  "CN"       Converged Newtonian light time
%                             correction. In solving the light time
%                             equation, the "CN" correction iterates
%                             until the solution converges. Both the
%                             state and rotation of the target body
%                             are corrected for light time.
%
%                  "CN+S"     Converged Newtonian light time
%                             and stellar aberration corrections.
%
%               The following values of 'abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at 'et' and arrive at the
%               intercept point at the light-time corrected epoch
%               ET+LT:
%
%                  "XLT"      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation. This correction yields the
%                             intercept location at the moment it
%                             receives photons emitted from the
%                             observer's location at 'et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             "LT" option uses one iteration.
%
%                             Both the target state as seen by the
%                             observer, and rotation of the target
%                             body, are corrected for light time.
%
%                  "XLT+S"    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation  This option modifies the
%                             intercept obtained with the "XLT"
%                             option to account for the observer's
%                             velocity relative to the solar system
%                             barycenter.
%
%                  "XCN"      Converged Newtonian light time
%                             correction. This is the same as XLT
%                             correction but with further iterations
%                             to a converged Newtonian light time
%                             solution.
%
%                  "XCN+S"    "Transmission" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               Case and embedded blanks are not significant in 'abcorr'.
%
%      obsrvr   the name of the observing body. This is typically
%               a spacecraft, the earth, or a surface point on the earth.
%
%               [1,c5] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               'obsrvr' is case-insensitive, and leading and
%               trailing blanks in 'obsrvr' are not significant.
%               Optionally, you may supply a string containing the
%               integer ID code for the object. For example both
%               "MOON" and "301" are legitimate strings that indicate
%               the moon is the observer.
%
%      dref     the name of the reference frame relative to which
%               the ray's direction vector is expressed. This may be
%               any frame supported by the SPICE system, including
%               built-in frames (documented in the Frames Required
%               Reading) and frames defined by a loaded frame kernel
%               (FK). The string `dref' is case-insensitive, and
%               leading and trailing blanks in `dref' are not
%               significant.
%
%               [1,c6] = size(dref); char = class(dref)
%
%                  or
%
%               [1,1] = size(dref); cell = class(dref)
%
%               When `dref' designates a non-inertial frame, the
%               orientation of the frame is evaluated at an epoch
%               dependent on the frame's center and, if the center is
%               not the observer, on the selected aberration
%               correction. See the description of the direction
%               vector `dvec' for details.
%
%      dvec     the pointing vector emanating from the observer. The
%               intercept with the target body's surface of the ray
%               defined by the observer and 'dvec' is sought.
%
%               [3,1] = size(dvec); double = class(dvec)
%
%               'dvec' is specified relative to the reference frame
%               designated by 'dref'.
%
%               Non-inertial reference frames are treated as follows:
%               if the center of the frame is at the observer's
%               location, the frame is evaluated at 'et'. If the
%               frame's center is located elsewhere, then letting
%               'ltcent' be the one-way light time between the observer
%               and the central body associated with the frame, the
%               orientation of the frame is evaluated at et-ltcent,
%               et+ltcent, or 'et' depending on whether the requested
%               aberration correction is, respectively, for received
%               radiation, transmitted radiation, or is omitted.
%               'ltcent' is computed using the method indicated by
%               'abcorr'.
%
%   the call:
%
%      [ spoint, trgepc, srfvec, found] = cspice_sincpt( method, target, ...
%                                                        et,     fixref, ...
%                                                        abcorr, obsrvr, ...
%                                                        dref,   dvec)
%
%   returns:
%
%      spoint   the surface intercept point on the target body of the ray
%               defined by the observer and the direction vector. If the ray
%               intersects the target body in multiple points, the selected
%               intersection point is the one closest to the observer. The
%               output argument 'found' (see below) indicates whether an
%               intercept was found.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               'spoint' is expressed in Cartesian coordinates,
%               relative to the target body-fixed frame designated by
%               FIXFRM. The body-fixed target frame is evaluated at
%               the intercept epoch 'trgepc' (see description below).
%
%               When light time correction is used, the duration of
%               light travel between 'spoint' to the observer is
%               considered to be the one way light time. When both
%               light time and stellar aberration corrections are
%               used, 'spoint' is selected such that, when 'spoint' is
%               corrected for light time and the vector from the
%               observer to the light-time corrected location of
%               'spoint' is corrected for stellar aberration, the
%               resulting vector is parallel to the ray defined by
%               the observer's location and 'dvec'.
%
%               The components of 'spoint' are given in units of km.
%
%      trgepc   the "intercept epoch."
%
%               [1,1] = size(trgepc); double = class(trgepc)
%
%               This is the epoch at which the ray defined by 'obsrvr'
%               and 'dvec' intercepts the target surface at 'spoint'.
%               'trgepc' is defined as follows: letting 'lt' be the one-way
%               light time between the observer and the intercept point,
%               'trgepc' is the epoch et-lt, et+lt, or 'et' depending on
%               whether the requested aberration correction is, respectively,
%               for  received radiation, transmitted radiation, or
%               omitted. 'lt' is computed using the method indicated by
%               'abcorr'.
%
%               'trgepc' is expressed as seconds past J2000 TDB.
%
%      srfvec   the vector from the observer's position at 'et' to
%               'spoint'. 'srfvec' is expressed in the target body-fixed
%               reference frame designated by 'fixref', evaluated at
%               'trgepc'.
%
%               [3,1] = size(srfvec); double = class(srfvec)
%
%               The components of 'srfvec' are given in units of km.
%
%               One can use the CSPICE function cspice_vnorm to obtain the
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
%               frame 'ref' at time 'et', call the routine 'cspice_pxfrm2'.
%               Let 'xform' be the 3x3 matrix representing the
%               rotation from the reference frame 'fixref' at time
%               'trgepc' to the reference frame 'ref' at time 'et'. Then
%               'srfvec' can be transformed to the result 'refvec' as
%               follows:
%
%                  xform  = cspice_pxfrm2 ( fixref, ref, trgepc, et )
%                  refvec = xform * srfvec
%
%      found    a flag indicating whether or not the ray
%               intersects the target. If an intersection exists
%               'found' will return true If the ray misses
%               the target, 'found' will return false.
%
%               [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      The following program computes surface intercept points on
%      Mars for the boresight and FOV boundary vectors of the MGS MOC
%      narrow angle camera. The intercepts are computed for a single
%      observation epoch. Light time and stellar aberration corrections
%      are used. For simplicity, camera distortion is ignored.
%
%      %
%      %  Local variables
%      %
%      abcorr  = 'CN+S';
%      camera  = 'MGS_MOC_NA';
%      fixref  = 'IAU_MARS';
%      method  = 'Ellipsoid';
%      obsrvr  = 'MGS';
%      target  = 'Mars';
%      utc     = '2003 OCT 13 06:00:00 UTC';
%      NCORNR  = 4;
%
%      %
%      % Load kernel files.
%      %
%      cspice_furnsh( 'standard.tm' );
%      cspice_furnsh( { '/kernels/MGS/ik/moc20.ti',                 ...
%                       '/kernels/MGS/sclk/MGS_SCLKSCET.00061.tsc'  ,...
%                       '/kernels/MGS/spk/mgs_ext12_ipng_mgs95j.bsp',...
%                       '/kernels/MGS/ck/mgs_sc_ext12.bc' } )
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et = cspice_str2et( utc );
%
%      %
%      % Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
%      % ID code. Then look up the field of view (FOV)
%      % parameters.
%      %
%      [ camid, found ] = cspice_bodn2c( camera );
%
%      if ( ~found )
%         txt = sprintf( [ 'SPICE(NOTRANSLATION) ' ...
%                         'Could not find ID code for instrument %s.' ], ...
%                          camera);
%         error( txt )
%      end
%
%      %
%      % cspice_getfov will return the name of the camera-fixed frame
%      % in the string 'dref', the camera boresight vector in
%      % the array 'bsight', and the FOV corner vectors in the
%      % array 'bounds'.
%      %
%      [shape, dref, bsight, bounds] = cspice_getfov( camid, NCORNR);
%
%      fprintf (  ['\n' ...
%                  'Surface Intercept Locations for Camera\n'  ...
%                  'FOV Boundary and Boresight Vectors\n'      ...
%                  '\n'                                        ...
%                  '   Instrument:             %s\n'           ...
%                  '   Epoch:                  %s\n'           ...
%                  '   Aberration correction:  %s\n'           ...
%                  '\n'],                                      ...
%                  camera, utc, abcorr )
%
%      for i=1:NCORNR+1
%
%         if( i <= NCORNR )
%            fprintf( 'Corner vector %d\n\n', i)
%            dvec = bounds(:,i);
%         end
%
%         if ( i == (NCORNR + 1) )
%            fprintf( 'Boresight vector\n\n' )
%            dvec = bsight;
%         end
%
%         %
%         % Compute the surface intercept point using
%         % the specified aberration corrections.
%         %
%         [ spoint, trgepc, srfvec, found ] =                   ...
%                        cspice_sincpt( method, target,         ...
%                                       et,     fixref, abcorr, ...
%                                       obsrvr, dref,   dvec );
%         if( found )
%
%            %
%            % Compute range from observer to apparent intercept.
%            %
%            dist = vnorm( srfvec );
%
%            %
%            % Convert rectangular coordinates to planetocentric
%            % latitude and longitude. Convert radians to degrees.
%            %
%            [ radius, lon, lat ] = cspice_reclat( spoint );
%
%            lon = lon * cspice_dpr;
%            lat = lat * cspice_dpr;
%
%            %
%            % Display the results.
%            %
%            fprintf( '  Vector in %s frame = \n', dref )
%            fprintf( '   %18.10e %18.10e %18.10e\n', dvec );
%
%            fprintf( [ '\n'                                              ...
%                       '  Intercept:\n'                                  ...
%                       '\n'                                              ...
%                       '     Radius                   (km)  = %18.10e\n' ...
%                       '     Planetocentric Latitude  (deg) = %18.10e\n' ...
%                       '     Planetocentric Longitude (deg) = %18.10e\n' ...
%                       '     Range                    (km)  = %18.10e\n' ...
%                       '\n' ],                                           ...
%                        radius,  lat,  lon,  dist                          )
%         else
%            disp( 'Intercept not found.' )
%         end
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
%      Surface Intercept Locations for Camera
%      FOV Boundary and Boresight Vectors
%
%         Instrument:             MGS_MOC_NA
%         Epoch:                  2003 OCT 13 06:00:00 UTC
%         Aberration correction:  CN+S
%
%      Corner vector 1
%
%        Vector in MGS_MOC_NA frame =
%           1.8571383810e-06  -3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849411359e+03
%           Planetocentric Latitude  (deg) =  -4.8477481924e+01
%           Planetocentric Longitude (deg) =  -1.2347407905e+02
%           Range                    (km)  =   3.8898310366e+02
%
%      Corner vector 2
%
%        Vector in MGS_MOC_NA frame =
%           1.8571383810e-06   3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849396987e+03
%           Planetocentric Latitude  (deg) =  -4.8481636340e+01
%           Planetocentric Longitude (deg) =  -1.2339882297e+02
%           Range                    (km)  =   3.8897512129e+02
%
%      Corner vector 3
%
%        Vector in MGS_MOC_NA frame =
%          -1.8571383810e-06   3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849396899e+03
%           Planetocentric Latitude  (deg) =  -4.8481661910e+01
%           Planetocentric Longitude (deg) =  -1.2339882618e+02
%           Range                    (km)  =   3.8897466238e+02
%
%      Corner vector 4
%
%        Vector in MGS_MOC_NA frame =
%          -1.8571383810e-06  -3.8015622659e-03   9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849411271e+03
%           Planetocentric Latitude  (deg) =  -4.8477507498e+01
%           Planetocentric Longitude (deg) =  -1.2347408220e+02
%           Range                    (km)  =   3.8898264472e+02
%
%      Boresight vector
%
%        Vector in MGS_MOC_NA frame =
%           0.0000000000e+00   0.0000000000e+00   1.0000000000e+00
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849404102e+03
%           Planetocentric Latitude  (deg) =  -4.8479579822e+01
%           Planetocentric Longitude (deg) =  -1.2343645396e+02
%           Range                    (km)  =   3.8897573572e+02
%
%-Particulars
%
%   Given a ray defined by a direction vector and the location of an
%   observer, cspice_sincpt computes the surface intercept point of the ray
%   on a specified target body. cspice_sincpt also determines the vector
%   from the observer to the surface intercept point. If the ray
%   intersects the target in multiple locations, the intercept
%   closest to the observer is selected.
%
%   When aberration corrections are used, this routine finds the
%   value of `spoint' such that, if `spoint' is regarded as an ephemeris
%   object, after the selected aberration corrections are applied to
%   the vector from the observer to `spoint', the resulting vector is
%   parallel to the direction vector `dvec'.
%
%   This routine computes light time corrections using light time
%   between the observer and the surface intercept point, as opposed
%   to the center of the target. Similarly, stellar aberration
%   corrections done by this routine are based on the direction of
%   the vector from the observer to the light-time corrected
%   intercept point, not to the target center. This technique avoids
%   errors due to the differential between aberration corrections
%   across the target body. Therefore it's valid to use aberration
%   corrections with this routine even when the observer is very
%   close to the intercept point, in particular when the
%   observer-intercept point distance is much less than the
%   observer-target center distance. It's also valid to use stellar
%   aberration corrections even when the intercept point is near or
%   on the limb (as may occur in occultation computations using a
%   point target).
%
%   When comparing surface intercept point computations with results
%   from sources other than SPICE, it's essential to make sure the
%   same geometric definitions are used.
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
%      both represent a surface as a set of triangular plates, the
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
%         'DSK/<surface list>/UNPRIORITIZED'
%         'DSK/UNPRIORITIZED/<surface list>'
%         'UNPRIORITIZED/<surface list>/DSK'
%
%      The simplest form of the `method' argument specifying use of
%      DSK data is one that lacks a surface list, for example:
%
%         'DSK/UNPRIORITIZED'
%
%      For applications in which all loaded DSK data for the target
%      body are for a single surface, and there are no competing
%      segments, the above string suffices. This is expected to be
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
%      Double quotes are used to delimit the surface name because
%      it contains blank characters.
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
%         'DSK/UNPRIORITIZED/SURFACES = "Mars MEGDR 64 PIXEL/DEG", 3'
%
%
%      Round-off errors and mitigating algorithms
%      ------------------------------------------
%
%      When topographic data are used to represent the surface of a
%      target body, round-off errors can produce some results that
%      may seem surprising.
%
%      Note that, since the surface in question might have mountains,
%      valleys, and cliffs, the points of intersection found for
%      nearly identical sets of inputs may be quite far apart from
%      each other: for example, a ray that hits a mountain side in a
%      nearly tangent fashion may, on a different host computer, be
%      found to miss the mountain and hit a valley floor much farther
%      from the observer, or even miss the target altogether.
%
%      Round-off errors can affect segment selection: for example, a
%      ray that is expected to intersect the target body's surface
%      near the boundary between two segments might hit either
%      segment, or neither of them; the result may be
%      platform-dependent.
%
%      A similar situation exists when a surface is modeled by a set
%      of triangular plates, and the ray is expected to intersect the
%      surface near a plate boundary.
%
%      To avoid having the routine fail to find an intersection when
%      one clearly should exist, this routine uses two "greedy"
%      algorithms:
%
%         1) If the ray passes sufficiently close to any of the
%            boundary surfaces of a segment (for example, surfaces of
%            maximum and minimum longitude or latitude), that segment
%            is tested for an intersection of the ray with the
%            surface represented by the segment's data.
%
%            This choice prevents all of the segments from being
%            missed when at least one should be hit, but it could, on
%            rare occasions, cause an intersection to be found in a
%            segment other than the one that would be found if higher
%            precision arithmetic were used.
%
%         2) For type 2 segments, which represent surfaces as
%            sets of triangular plates, each plate is expanded very
%            slightly before a ray-plate intersection test is
%            performed. The default plate expansion factor is
%
%               1 + 1.e-10
%
%            In other words, the sides of the plate are lengthened by
%            1/10 of a micron per km. The expansion keeps the centroid
%            of the plate fixed.
%
%            Plate expansion prevents all plates from being missed
%            in cases where clearly at least one should be hit.
%
%            As with the greedy segment selection algorithm, plate
%            expansion can occasionally cause an intercept to be
%            found on a different plate than would be found if higher
%            precision arithmetic were used. It also can occasionally
%            cause an intersection to be found when the ray misses
%            the target by a very small distance.
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
%   the CSPICE routine sincpt_c.
%
%   MICE.REQ
%   DSK.REQ
%   FRAMES.REQ
%   NAIF_IDS.REQ
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
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.3, 12-MAR-2012, SCK (JPL)
%
%       References to the new 'cspice_pxfrm2' routine were
%       added to the 'I/O returns' section. A problem description was
%       added to the 'Examples' section, and the references to
%       'cspice_srfxpt' and the second example were removed.
%
%   -Mice Version 1.0.2, 14-JUL-2010, EDW (JPL)
%
%       Corrected minor typo in header.
%
%   -Mice Version 1.0.1, 23-FEB-2009, EDW (JPL)
%
%       Added proper markers for usage string variable types.
%
%   -Mice Version 1.0.0, 11-FEB-2008, EDW (JPL)
%
%-Index_Entries
%
%   find surface intercept point
%   find intersection of ray and target body surface
%   find intercept of ray on target body surface
%
%-&

function [ spoint, trgepc, srfvec, found] = ...
          cspice_sincpt( method, target, et, fixref, abcorr, obsrvr, dref, dvec)

   switch nargin
      case 8

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);

      otherwise

         error( [ 'Usage: [ spoint, trgepc, srfvec, found] =  ' ...
                  'cspice_sincpt( `method`, `target`, et, `fixref`, ' ...
                                 '`abcorr`, `obsrvr`, `dref`, dvec)' ]  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [sincpt] = mice('sincpt_s', method, target, ...
                                  et, fixref, abcorr, obsrvr, dref, dvec);
      spoint = reshape( [sincpt.spoint], 3, [] );
      trgepc = reshape( [sincpt.trgepc], 1, [] );
      srfvec = reshape( [sincpt.srfvec], 3, [] );
      found  = reshape( [sincpt.found] , 1, [] );
   catch
      rethrow(lasterror)
   end

