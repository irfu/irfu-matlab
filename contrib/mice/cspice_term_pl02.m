%-Abstract
%
%   Deprecated: This routine has been superseded by the CSPICE routine
%   cspice_termpt. This routine is supported for purposes of backward
%   compatibility only.
%
%   CSPICE_TERM_PL02 computes a set of points on the umbral or penumbral
%   terminator of a specified target body, where the target body's surface
%   is represented by a triangular plate model contained in a type 2 DSK
%   segment.
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
%      handle      the DAS file handle of a DSK file open for read
%                  access. This kernel must contain a type 2 segment
%                  that provides a plate model representing the entire
%                  surface of the target body.
%
%                  [1,1] = size(handle); int32 = class(handle)
%
%      dladsc      the DLA descriptor of a DSK segment representing
%                  the surface of a target body.
%
%                 [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                 int32 = class(dladsc)
%
%      trmtyp      a string indicating the type of terminator to
%                  compute:  umbral or penumbral. The umbral terminator is
%                  the boundary of the portion of the target surface in
%                  total shadow. The penumbral terminator is the boundary
%                  of the portion of the surface that is completely
%                  illuminated. Note that in astronomy references, the
%                  unqualified word "terminator" refers to the umbral
%                  terminator. Here, the unqualified word refers to either
%                  type of terminator.
%
%                  [1,c1] = size(trmtyp); char = class(trmtyp)
%
%                     or
%
%                  [1,1] = size(trmtyp); cell = class(trmtyp)
%
%                  To compute the terminator points, this routine first
%                  computes a set of points on the terminator of the
%                  indicated type on the surface of a reference ellipsoid
%                  for the target body. Each such point defines the
%                  direction of a ray emanating from the target center and
%                  associated with a terminator point on the actual surface
%                  defined by the plate model. The outermost surface
%                  intercept of each such ray is a considered to be a
%                  terminator point of the surface defined by the plate
%                  model.
%
%                  Possible values of `trmtyp' are
%
%                     'UMBRAL'
%                     'PENUMBRAL'
%
%                  Case and leading or trailing blanks in `trmtyp' are
%                  not significant.
%
%      source      the name of the body acting as a light source.
%                  `source' is case-insensitive, and leading and trailing
%                  blanks in `target' are not significant. Optionally, you
%                  may supply a string containing the integer ID code
%                  for the object. For example both 'SUN' and '10' are
%                  legitimate strings that indicate the Sun is the light
%                  source.
%
%                  [1,c2] = size(source); char = class(source)
%
%                     or
%
%                  [1,1] = size(source); cell = class(source)
%
%                  This routine assumes that a kernel variable
%                  representing the light source's radii is present in
%                  the kernel pool. Normally the kernel variable would
%                  be defined by loading a PCK file.
%
%                  The shape of the light source is always modeled as a
%                  sphere, regardless of whether radii defining a
%                  triaxial ellipsoidal shape model are available in the
%                  kernel pool. The maximum radius of the body is used
%                  as the radius of the sphere.
%
%      target      the name of the target body. `target' is
%                  case-insensitive, and leading and trailing blanks in
%                  `target' are not significant. Optionally, you may supply
%                  a string containing the integer ID code for the object.
%                  For example both 'MOON' and '301' are legitimate strings
%                  that indicate the moon is the target body.
%
%                  [1,c3] = size(target); char = class(target)
%
%                     or
%
%                  [1,1] = size(target); cell = class(target)
%
%                  This routine assumes that a kernel variable representing
%                  the target's radii is present in the kernel pool.
%                  Normally the kernel variable would be defined by loading
%                  a PCK file.
%
%      et          the epoch of participation of the observer,
%                  expressed as ephemeris seconds past J2000 TDB: `et' is
%                  the epoch at which the observer's position is
%                  computed.
%
%                  [1,1] = size(et); double = class(et)
%
%                  When aberration corrections are not used, `et' is also
%                  the epoch at which the position and orientation of the
%                  target body and position of the light source are
%                  computed.
%
%                  When aberration corrections are used, `et' is the epoch
%                  at which the observer's position relative to the solar
%                  system barycenter is computed; in this case the
%                  position and orientation of the target body are
%                  computed at et-lt or et+lt, where `lt' is the one-way
%                  light time between the target body's center and the
%                  observer, and the sign applied to `lt' depends on the
%                  selected correction. See the description of `abcorr'
%                  below for details.
%
%      fixfrm      the name of the reference frame relative to which
%                  the output terminator points are expressed. This must
%                  a body-centered, body-fixed frame associated with the
%                  target. The frame's axes must be compatible with the
%                  triaxial ellipsoidal shape model associated with the
%                  target body (normally provide via a PCK): this
%                  routine assumes that the first, second, and third
%                  axis lengths correspond, respectively, to the x, y,
%                  and z-axes of the frame designated by `fixfrm'.
%
%                  [1,c4] = size(fixfrm); char = class(fixfrm)
%
%                     or
%
%                  [1,1] = size(fixfrm); cell = class(fixfrm)
%
%                  `fixfrm' may refer to a built-in frame (documented in
%                  the Frames Required Reading) or a frame defined by a
%                  loaded frame kernel (FK).
%
%                  The orientation of the frame designated by `fixfrm' is
%                  evaluated at epoch of participation of the target
%                  body. See the descriptions of `et' and `abcorr' for
%                  details.
%
%
%      abcorr      indicates the aberration correction to be applied
%                  when computing the observer-target position, the
%                  orientation of the target body, and the target-
%                  source position vector.
%
%                  [1,c5] = size(abcorr); char = class(abcorr)
%
%                     or
%
%                  [1,1] = size(abcorr); cell = class(abcorr)
%
%
%                  `abcorr' may be any of the following.
%
%                     'NONE'     Apply no correction. Compute the
%                                terminator points using the position
%                                of the light source and target, and
%                                the orientation of the target, at `et'.
%
%                  Let `lt' represent the one-way light time between the
%                  observer and the target body's center. The following
%                  values of `abcorr' apply to the "reception" case in
%                  which photons depart from the target body's center at
%                  the light-time corrected epoch et-lt and *arrive* at
%                  the observer's location at `et':
%
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                yields the location of the terminator
%                                points at the approximate time they
%                                emitted photons arriving at the
%                                observer at `et' (the difference between
%                                light time to the target center and
%                                light time to the terminator points
%                                is ignored).
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                The target position as seen by the
%                                observer, the position of the light
%                                source as seen from the target at
%                                et-lt, and the rotation of the target
%                                body, are corrected for light time.
%
%                     'LT+S'     Correct for one-way light time and
%                                stellar aberration using a Newtonian
%                                formulation. This option modifies the
%                                positions obtained with the 'LT' option
%                                to account for the observer's velocity
%                                relative to the solar system
%                                barycenter. This correction also
%                                applies to the position of the light
%                                source relative to the target. The
%                                result is the apparent terminator as
%                                seen by the observer.
%
%                     'CN'       Converged Newtonian light time
%                                correction. In solving the light time
%                                equation, the 'CN' correction iterates
%                                until the solution converges. The
%                                position and rotation of the target
%                                body and the position of the light
%                                source relative to the target are
%                                corrected for light time.
%
%                     'CN+S'     Converged Newtonian light time
%                                and stellar aberration corrections.
%
%
%      obsrvr      the name of the observing body. This is typically
%                  a spacecraft, the Earth, or a surface point on the
%                  Earth. `obsrvr' is case-insensitive, and leading and
%                  trailing blanks in `obsrvr' are not significant.
%                  Optionally, you may supply a string containing the
%                  integer ID code for the object. For example both
%                  'EARTH' and '399' are legitimate strings that indicate
%                  the Earth is the observer.
%
%                  [1,c6] = size(obsrvr); char = class(obsrvr)
%
%                     or
%
%                  [1,1] = size(obsrvr); cell = class(obsrvr)
%
%      npoints     the number of terminator points to compute.
%
%   the call:
%
%      [trgepc, obspos, termpts, plateids] =                 ...
%                  cspice_term_pl02( handle, dladsc,         ...
%                                    trmtyp, source, target, ...
%                                    et,     fixref, abcorr, ...
%                                    obsrvr, npoints )
%
%   returns:
%
%      trgepc      the "target epoch."  `trgepc' is defined as follows:
%                  letting `lt' be the one-way light time between the
%                  target center and observer, `trgepc' is either the
%                  epoch et-lt or `et' depending on whether the requested
%                  aberration correction is, respectively, for received
%                  radiation or omitted. `lt' is computed using the
%                  method indicated by `abcorr'.
%
%                  [1,1] = size(trgepc); double = class(trgepc)
%
%                  `trgepc' is expressed as seconds past J2000 TDB.
%
%      obspos      the vector from the center of the target body at
%                  epoch `trgepc' to the observer at epoch `et'.  `obspos' is
%                  expressed in the target body-fixed reference frame
%                  `fixfrm', which is evaluated at `trgepc'.
%
%                  [3,1] = size(obspos); double = class(obspos)
%
%                  `obspos' is returned to simplify various related
%                  computations that would otherwise be cumbersome.  For
%                  example, the vector `xvec' from the observer to the
%                  Ith terminator point can be calculated via the
%                  expression
%
%                     xvec = termpts(:,i) - obspos
%
%                  The components of `obspos' are given in units of km.
%
%      termpts     an array of points on the umbral or penumbral
%                  terminator of the target, as specified by the input
%                  argument `trmtyp'.
%
%                  [3,npoints] = size(termpts); double = class(termpts)
%
%                  The ith point is contained in the array elements
%
%                      termpts(:,i)
%
%                  As described above, each terminator point lies on a ray
%                  emanating from the center of the target and passing
%                  through a terminator point on the target's reference
%                  ellipsoid. Each terminator point *on the reference
%                  ellipsoid* is the point of tangency of a plane that is
%                  also tangent to the light source. These associated
%                  points of tangency on the light source have uniform
%                  distribution in longitude when expressed in a
%                  cylindrical coordinate system whose Z-axis is `obspos'.
%                  The magnitude of the separation in longitude between the
%                  tangency points on the light source is
%
%                     2*Pi / npoints
%
%                  If the reference ellipsoid for the target is spherical,
%                  the terminator points also are uniformly distributed in
%                  longitude in the cylindrical system described above.  If
%                  the reference ellipsoid of the target is non-spherical,
%                  the longitude distribution of the points generally is
%                  not uniform.
%
%                  The terminator points are expressed in the body-fixed
%                  reference frame designated by `fixfrm'. Units are km.
%
%      plateids    an array of integer ID codes of the plates on which
%                  the terminator points are located.  The ith plate ID
%                  corresponds to the ith terminator point. These ID codes can
%                  be use to look up data associated with the plate, such
%                  as the plate's vertices or outward normal vector.
%
%                  [1,npoints] = size(plateIDs); int32 = class(plateIDs)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Compute sets of umbral and penumbral terminator points on Phobos
%   as seen from Mars. Perform a consistency check using the solar
%   incidence angle at each point, where the solar incidence angle
%   is computed using both a reference ellipsoid and the actual
%   plate model surface and surface normal. We expect to see a
%   solar incidence angle of approximately 90 degrees. Since the
%   solar incidence angle is measured between the local outward
%   normal and the direction to the Sun, the solar incidence angle
%   at an umbral or penumbral terminator point should be,
%   respectively, greater than or less than 90 degrees by
%   approximately the angular radius of the Sun as seen from each
%   terminator point.
%
%   In the following example program, the file
%
%      phobos_3_3.dsk
%
%   is a DSK file containing a type 2 segment that provides a plate model
%   representation of the surface of Phobos. The file
%
%      mar097.bsp
%
%   is a binary SPK file containing data for Phobos, the Earth, and the
%   Sun for a time interval bracketing the date
%
%      2006 NOV 3 00:00:00 UTC.
%
%   pck00010.tpc is a planetary constants kernel file containing radii
%   and rotation model constants. naif0010.tls is a leapseconds kernel.
%
%   All of the kernels other than the DSK file should be loaded via
%   a meta-kernel. An example of the contents of such a kernel is:
%
%         KPL/MK
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0010.tls'
%                                'pck00010.tpc'
%                                'mar097.bsp' )
%         \begintext
%
%
%      function term_pl02_t
%
%         %
%         % Constants
%         %
%         NPOINTS     = 3;
%         NTYPES      = 2;
%         TIMLEN      = 40;
%         ILUM_METHOD = 'ELLIPSOID';
%         TOL         = 1.d-12;
%         UTCSTR      = '2007 FEB 9 00:00:00 UTC';
%
%         %
%         % Initial values
%         %
%         target      = 'Phobos';
%         fixfrm      = 'IAU_PHOBOS';
%         abcorr      = 'CN+S';
%         fixfrm      = 'IAU_PHOBOS';
%         obsrvr      = 'Mars';
%         trmtypes    = { 'Umbral', 'Penumbral' };
%         utcstr      = '2007 FEB 9 00:00:00 UTC';
%
%         %
%         % Prompt for the name of a meta-kernel specifying
%         % all of the other kernels we need. Load the
%         % metakernel.
%         %
%         meta = input( 'Enter meta-kernel name > ','s');
%         cspice_furnsh( meta )
%
%         %
%         % Prompt for the name of a DSK file.
%         %
%         dsk = input( 'Enter DSK name         > ','s');
%
%         %
%         % Open the DSK file for read access.
%         % We use the DAS-level interface for
%         % this function.
%         %
%         handle = cspice_dasopr( dsk );
%
%         %
%         % Begin a forward search through the
%         % kernel, treating the file as a DLA.
%         % In this example, it's a very short
%         % search.
%         %
%         [dladsc, found] = cspice_dlabfs( handle );
%
%         if ~found
%
%            %
%            % We arrive here only if the kernel
%            % contains no segments. This is
%            % unexpected, but we're prepared for it.
%            %
%            fprintf( 'No segments found in DSK file %s\n', dsk )
%            return
%
%         end
%
%         %
%         % If we made it this far, `dladsc' is the
%         % DLA descriptor of the first segment.
%         %
%         % Now compute sub-points using both computation
%         % methods. We'll vary the aberration corrections
%         % and the epochs.
%         %
%
%         et = cspice_str2et( UTCSTR );
%
%         timstr = cspice_timout( et,                         ...
%                        'YYYY-MON-DD HR:MN:SC.### ::TDB(TDB)' );
%
%
%         fprintf( '\n\n' )
%         fprintf( '   Observer:                %s\n', obsrvr )
%         fprintf( '   Target:                  %s\n', target )
%         fprintf( '   Observation epoch:       %s\n', timstr )
%         fprintf( '   Aberration correction:   %s\n', abcorr )
%         fprintf( '   Body-fixed frame:        %s\n', fixfrm )
%
%
%         %
%         % Look up the radii of the Sun. We'll use these as
%         % part of a computation to check the solar incidence
%         % angles at the terminator points.
%         %
%         sunRadii = cspice_bodvrd( 'SUN', 'RADII', 3 );
%
%         %
%         % Now compute grids of terminator points using both
%         % terminator types.
%         %
%
%         for  typidx = 1:NTYPES
%
%            %
%            % Select the terminator type.
%            %
%            trmtyp = trmtypes( typidx );
%
%            fprintf( '\n   Terminator type: %s\n', char(trmtyp) )
%
%            %
%            % Compute the terminator point set.
%            %
%            [trgepc, obpos, termpts, plateids] = ...
%                     cspice_term_pl02( handle,  dladsc,         ...
%                                       trmtyp,  'Sun',  target, ...
%                                       et,      fixfrm, abcorr, ...
%                                       obsrvr,  NPOINTS           );
%
%            %
%            % Display the terminator points.
%            %
%            for  i = 1:NPOINTS
%
%               [radius, lon, lat] = cspice_reclat( termpts(:,i) );
%
%               fprintf( '\n      Terminator point: %d\n', i )
%               fprintf( '         Radius                     (km): %f\n', ...
%                                                                 radius )
%
%               fprintf( '         Planetocentric longitude   (deg): %f\n', ...
%                                                     lon * cspice_dpr() )
%
%               fprintf( '         Planetocentric latitude    (deg): %f\n', ...
%                                                     lat * cspice_dpr() )
%
%               fprintf( '         Plate ID:                         %d\n', ...
%                                                            plateids(i) )
%
%               %
%               % Compute the angular radius of the Sun as seen from
%               % the current terminator point. Subtracting (adding)
%               % this value from (to) the solar incidence angle for
%               % umbral (penumbral) terminator points should yield a
%               % value close to 90 degrees. This provides a sanity
%               % check on the locations of the terminator points.
%               %
%               % First find the position of the Sun relative to the
%               % target's center at the light time corrected epoch
%               % trgepc.
%               %
%               [sunPos, ltime] = cspice_spkpos( 'Sun',  trgepc, fixfrm, ...
%                                                abcorr, target );
%
%               sunVec    = sunPos - termpts(:,i);
%
%               sunAngRad = asin( sunRadii(1) / cspice_vnorm(sunVec) );
%
%               %
%               % Compute the delta by which we adjust the solar
%               % incidence angles.
%               %
%               if  ( typidx == 1 )
%
%                  %
%                  % Umbral
%                  %
%                  delta = -sunAngRad;
%
%               else
%
%
%                  %
%                  % Penumbral
%                  %
%                  delta =  sunAngRad;
%
%               end
%
%               %
%               % Compute the illumination angles using an ellipsoidal
%               % representation of the target's surface. The role of
%               % this representation is to provide an outward surface
%               % normal.
%               %
%               [trgepc, srfvec, phase,  solar,  emissn] =        ...
%                         cspice_ilumin( ILUM_METHOD, target, et, ...
%                                        fixfrm, abcorr, obsrvr,  ...
%                                        termpts(:,i) );
%
%               fprintf( '            Solar incidence angle derived using\n' )
%               fprintf( '            '                           )
%               fprintf( '   - an ellipsoidal reference surface'  )
%               fprintf( '        (deg): %f\n', solar * cspice_dpr() )
%
%               fprintf( '            '                           )
%               fprintf( '        > adjusted for Solar angular '  )
%               fprintf( 'radius  (deg): %f\n', (solar+delta) * cspice_dpr() )
%
%               %
%               % Compute the illumination angles at the terminator point
%               % using the actual plate model surface normal.
%               %
%               [phase, solar, emissn] = cspice_illum_pl02(      ...
%                                        handle, dladsc, target, ...
%                                        et,     abcorr, obsrvr, ...
%                                        termpts(:,i)  );
%
%               fprintf( '            '                            )
%               fprintf( '   - plate model''s surface and normal ' )
%               fprintf( 'vector (deg): %f\n', solar * cspice_dpr() )
%
%            end
%
%         end
%
%         %
%         % Close the DSK file. Unload all other kernels as well.
%         %
%         cspice_dascls( handle )
%
%         cspice_kclear
%
%   MATLAB outputs:
%
%         Observer:                Mars
%         Target:                  Phobos
%         Observation epoch:       2007-FEB-09 00:01:05.184 (TDB)
%         Aberration correction:   CN+S
%         Body-fixed frame:        IAU_PHOBOS
%
%         Terminator type: Umbral
%
%            Terminator point: 1
%               Radius                     (km): 12.111257
%               Planetocentric longitude   (deg): 34.584501
%               Planetocentric latitude    (deg): -0.001298
%               Plate ID:                         200400
%                  Solar incidence angle derived using
%                    - an ellipsoidal reference surface        (deg): 90.182028
%                         > adjusted for Solar angular radius  (deg): 89.999999
%                    - plate model's surface and normal vector (deg): 90.240660
%
%            Terminator point: 2
%               Radius                     (km): 9.774665
%               Planetocentric longitude   (deg): -143.659941
%               Planetocentric latitude    (deg): 43.397190
%               Plate ID:                         156958
%                  Solar incidence angle derived using
%                    - an ellipsoidal reference surface        (deg): 90.182028
%                         > adjusted for Solar angular radius  (deg): 90.000000
%                    - plate model's surface and normal vector (deg): 87.138686
%
%            Terminator point: 3
%               Radius                     (km): 11.500619
%               Planetocentric longitude   (deg): -146.128151
%               Planetocentric latitude    (deg): -43.082379
%               Plate ID:                         25552
%                  Solar incidence angle derived using
%                    - an ellipsoidal reference surface        (deg): 90.182028
%                         > adjusted for Solar angular radius  (deg): 90.000000
%                    - plate model's surface and normal vector (deg): 91.404206
%
%         Terminator type: Penumbral
%
%            Terminator point: 1
%               Radius                     (km): 12.859785
%               Planetocentric longitude   (deg): -145.415505
%               Planetocentric latitude    (deg): 0.001299
%               Plate ID:                         86763
%                  Solar incidence angle derived using
%                    - an ellipsoidal reference surface        (deg): 89.817971
%                         > adjusted for Solar angular radius  (deg): 90.000000
%                    - plate model's surface and normal vector (deg): 89.055489
%
%            Terminator point: 2
%               Radius                     (km): 10.327413
%               Planetocentric longitude   (deg): 36.340069
%               Planetocentric latitude    (deg): -43.397192
%               Plate ID:                         76977
%                  Solar incidence angle derived using
%                    - an ellipsoidal reference surface        (deg): 89.817971
%                         > adjusted for Solar angular radius  (deg): 90.000000
%                    - plate model's surface and normal vector (deg): 77.351956
%
%            Terminator point: 3
%               Radius                     (km): 10.086025
%               Planetocentric longitude   (deg): 33.871859
%               Planetocentric latitude    (deg): 43.082380
%               Plate ID:                         282136
%                  Solar incidence angle derived using
%                    - an ellipsoidal reference surface        (deg): 89.817971
%                         > adjusted for Solar angular radius  (deg): 90.000000
%                    - plate model's surface and normal vector (deg): 88.997322
%
%-Particulars
%
%   In this routine, we use the term "umbral terminator" to denote
%   the curve usually called the "terminator":  this curve is the
%   boundary of the portion of the target body's surface that lies in
%   total shadow. We use the term "penumbral terminator" to denote
%   the boundary of the completely illuminated portion of the
%   surface.
%
%   Boundaries of illuminated regions on an arbitrary surface are often
%   complicated point sets:  boundaries of shadows of mountains and
%   craters, if present, all contribute to the overall set. To make the
%   terminator computation tractable, we simplify the problem by using a
%   reference ellipsoid for guidance. We compute a set of terminator
%   points on the reference ellipsoid for the target body, then use
%   those points to define the latitudes and longitudes of terminator
%   points on the surface defined by the specified triangular shape
%   model. As such, the set of terminator points found by this routine
%   is just an approximation.
%
%   Below we discuss the computation of terminator points on the target
%   body's reference ellipsoid.
%
%   This routine assumes a spherical light source. Light rays are
%   assumed to travel along straight lines; refraction is not modeled.
%
%   Points on the reference ellipsoid at which the entire cap of
%   the light source is visible are considered to be completely
%   illuminated. Points on the ellipsoid at which some portion
%   (or all) of the cap of the light source are blocked are
%   considered to be in partial (or total) shadow.
%
%   In general, the terminator on an ellipsoid is a more complicated
%   curve than the limb (which is always an ellipse). Aside from
%   various special cases, the terminator does not lie in a plane.
%
%   However, the condition for a point X on the ellipsoid to lie on
%   the terminator is simple:  a plane tangent to the ellipsoid at X
%   must also be tangent to the light source. If this tangent plane
%   does not intersect the vector from the center of the ellipsoid to
%   the center of the light source, then X lies on the umbral
%   terminator; otherwise X lies on the penumbral terminator.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine term_pl02.
%
%   MICE.REQ
%   ABCORR.REQ
%   DSK.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 25-JUL-2016, NJB, EDW (JPL)
%
%-Index_Entries
%
%   find terminator on plate model
%   find terminator on triangular shape model
%   find terminator on dsk type_2 shape model
%   find umbral terminator on plate model
%   find umbral terminator on triangular shape model
%   find umbral terminator on dsk type_2 shape model
%   find penumbral terminator on plate model
%   find penumbral terminator on triangular shape model
%   find penumbral terminator on dsk type_2 shape model
%
%-&

function [trgepc, obspos, termpts, plateIDs] = ...
                     cspice_term_pl02( handle, dladsc, ...
                                       trmtyp, source, target, ...
                                       et,     fixref, abcorr, ...
                                       obsrvr, npoints           )

   switch nargin
      case 10

         handle  = zzmice_int( handle );
         dladsc  = zzmice_int( dladsc );
         trmtyp  = zzmice_str( trmtyp );
         source  = zzmice_str( source );
         target  = zzmice_str( target );
         et      = zzmice_dp( et );
         fixref  = zzmice_str( fixref );
         abcorr  = zzmice_str( abcorr );
         obsrvr  = zzmice_str( obsrvr );
         npoints = zzmice_int( npoints );

      otherwise

         error ( [ 'Usage: [trgepc, obspos(3), '                         ...
                  'termpts(3,npoints), plateIDs(npoints)] = '            ...
                  'cspice_term_pl02( handle, dladsc(SPICE_DLA_DSCSIZ), ' ...
                  '`trmtyp`, `source`, `target`, '                       ...
                  'et, `fixref`, `abcorr`, `obsrvr`, npoints )' ] )

   end

   %
   % Call the MEX library.
   %
   try

      [trgepc, obspos, termpts, plateIDs] = mice( 'term_pl02', ...
                                   handle, dladsc, trmtyp, source, target, ...
                                   et,     fixref, abcorr, obsrvr, npoints );
   catch
      rethrow(lasterror)
   end


