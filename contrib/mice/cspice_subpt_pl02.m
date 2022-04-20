%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   cspice_subpnt. This routine is supported for purposes of backward
%   compatibility only.
%
%   CSPICE_SUBPT_PL02 returns the rectangular coordinates of the sub-observer
%   point on a target body in the body-fixed frame associated with the body, 
%   at a particular epoch, optionally corrected for light time and stellar 
%   aberration. In addition, it returns the observer's altitude above the 
%   target body.
%
%   The target body's surface is represented by a triangular plate model 
%   contained in a type 2 DSK segment.
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
%      handle   the DAS file handle of a DSK file open for read
%               access.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               This kernel must contain a type 2 segment that
%               provides a plate model representing the entire surface
%               of the target body.
%
%      dladsc   the DLA descriptor of a DSK segment representing the
%               surface of the target body.
%
%               [SPICE_DLA_DSCSIZ,1] = size(dladsc); int32 = class(dladsc)
%
%      method   a short string specifying the computation method
%               to be used.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               The choices are:
%
%                  'Intercept'       The sub-observer point is defined as
%                                    the plate model surface intercept of
%                                    the ray starting at the observer and
%                                    passing through the target's center.
%
%                  'Ellipsoid
%                   near point'      The sub-observer point is defined
%                                    as the plate model surface
%                                    intercept of the ray starting at
%                                    the observer and passing through
%                                    the nearest point to the observer
%                                    on a reference ellipsoid
%                                    associated with the target body.
%
%                                    This option requires that the
%                                    reference ellipsoid's radii be
%                                    available in the kernel pool.
%
%               For both computation methods, this routine finds a
%               sub-point on the same side of the target body as the
%               observer. If the observer is inside the target body,
%               the 'sub-point' will actually be above the observer.
%               In the case of multiple intercepts, the outermost one
%               (that is, the one farthest from the target center) is
%               selected.
%
%               Neither case nor white space are significant in the
%               string `method'. For example, the string
%
%                  '  ellipsoidNEARPOINT'
%
%               is valid.
%
%      target   the name of the target body.
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
%               This routine assumes that the target body's surface is
%               represented by a plate model, and that a DSK file
%               containing the plate model has been loaded via cspice_dasopr.
%
%      et       the epoch, represented  as seconds past J2000 TDB, at
%               which the sub-observer point on the target body is to be
%               computed.
%
%               [1,1] = size(et); double = class(et)
%
%               When aberration corrections are used, `et' refers to the epoch
%               at which radiation is received at the observer.
%
%      abcorr   indicates the aberration corrections to be applied to
%               the position and orientation of the target body and the
%               position of the Sun to account for one-way light time
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
%                  'NONE'     Apply no correction. Use the geometric
%                             positions of the Sun and target body
%                             relative to the observer; evaluate the
%                             target body's orientation at `et'.
%
%               The following values of `abcorr' apply to the
%               "reception" case in which photons depart from the
%               target's location at the light-time corrected epoch
%               et-lt and *arrive* at the observer's location at
%               `et':
%
%                  'LT'       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             uses the position and orientation of the
%                             target at the moment it emitted photons
%                             arriving at the observer at `et'. The
%                             position of the Sun relative to the
%                             target is corrected for the one-way light
%                             time from the Sun to the target.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation (see -Particulars for details).
%                             The solution invoked by the 'LT' option
%                             uses one iteration.
%
%                  'LT+S'     Correct for one-way light time and stellar
%                             aberration using a Newtonian formulation.
%                             This option modifies the positions
%                             obtained with the 'LT' option to account
%                             for the observer's velocity relative to
%                             the solar system barycenter (note the
%                             target plays the role of "observer" in the
%                             computation of the aberration-corrected
%                             target-Sun vector). The result is the
%                             sub-solar point computed using apparent
%                             position and orientation of the target as
%                             seen by the observer and the apparent
%                             position of the Sun as seen by the target.
%
%                  'CN'       Converged Newtonian light time correction.
%                             In solving the light time equation, the
%                             'CN' correction iterates until the
%                             solution converges (three iterations on
%                             all supported platforms).
%
%                  'CN+S'     Converged Newtonian light time
%                             and stellar aberration corrections.
%
%      obsrvr   the name of the observing body.
%
%               [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               This is typically a spacecraft, the earth, or a surface point
%               on the earth. `obsrvr' is case-insensitive, and leading and
%               trailing blanks in `obsrvr' are not significant. Optionally,
%               you may supply a string containing the integer ID code for
%               the object. For example both 'EARTH' and '399' are
%               legitimate strings that indicate the earth is the
%               observer.
%
%   the call:
%
%      [spoint, alt, pltid] = cspice_subpt_pl02( handle, dladsc, method,   ...
%                                                target, et,     abcorr,   ...
%                                                obsrvr                )
%
%   returns:
%
%      spoint   the sub-observer point on the target body expressed
%               relative to the body-fixed reference frame of the target
%               body.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               The definition of sub-observer point depends on the
%               selected computation method. See the description of the
%               input argument `method' for details.
%
%               The target body-fixed frame, which is time-dependent, is
%               evaluated at `et' if `abcorr' is 'NONE'; otherwise the
%               frame is evaluated at et-lt, where `lt' is the one-way
%               light time from target to observer.
%
%               The position and orientation of the target body are
%               corrected for aberration as specified by `abcorr'; the
%               corrected position and orientation are used in the
%               computation of `spoint'.
%
%      alt      the signed distance between the observer and the
%               sub-point.
%
%               [1,1] = size(alt); double = class(alt)
%
%               When the observer is outside the body `alt' is positive; when
%               the observer is inside, `alt' is negative.
%
%               Note that `alt' is not truly an "altitude" unless the
%               observer-to-sub-point vector happens to be perpendicular
%               to the target body's surface at the sub-point. In
%               general this condition should not be expected to hold,
%               unless the plate model representation of the target
%               body's surface very nearly matches the target body's
%               reference ellipsoid and the 'ellipsoid near point'
%               computation method is selected.
%
%      pltid    the integer ID code of the plate on which the
%               sub-observer point is located.
%
%               [1,1] = size(pltid); int32 = class(pltid)
%
%               This ID code can be use to look up data associated with the
%               plate, such as the plate's vertices or outward normal vector.
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
%   1) Find the sub-observer point on Phobos as seen from Earth for a
%      specified sequence of times. Perform the computation twice,
%      using both the "intercept" and "ellipsoid near point"
%      options. Compute the corresponding sub-subobserver point values
%      using an ellipsoidal surface for comparison.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: subpt_pl02_ex1.tm
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
%            mar097.bsp                       Mars satellite ephemeris
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0010.tls                     Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'mar097.bsp',
%                                'pck00010.tpc',
%                                'naif0010.tls' )
%         \begintext
%
%         End of meta-kernel
%
%
%      Use the DSK kernel below to provide the plate model representation
%      of the surface of Phobos.
%
%         phobos_3_3.bds
%
%
%
%      Example code begins here.
%
%
%      function subpt_pl02_ex1( meta, dsknam )
%
%         %
%         % Constants
%         %
%         NCORR       = 2;
%         NSAMP       = 3;
%         NMETHOD     = 2;
%         TIMLEN      = 40;
%         FIXREF      = 'IAU_PHOBOS';
%         TOL         = 1.d-12;
%
%         %
%         % Initial values
%         %
%         abcorrs     = { 'NONE', 'CN+S' };
%         emethods     = { 'Intercept: ellipsoid', 'Near point: ellipsoid' };
%         plmethods   = { 'Intercept', 'Ellipsoid near point' };
%
%         obsrvr      = 'Earth';
%         target      = 'Phobos';
%
%         %
%         % Load the metakernel.
%         %
%         cspice_furnsh( meta )
%
%         %
%         % Open the DSK file for read access.
%         % We use the DAS-level interface for
%         % this function.
%         %
%         handle = cspice_dasopr( dsknam );
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
%            fprintf( 'No segments found in DSK file %s\n', dsknam )
%            return
%
%         end
%
%         %
%         % If we made it this far, `dladsc' is the
%         % DLA descriptor of the first segment.
%         %
%         % Now compute sub-points using both computation
%         % methods.  We'll vary the aberration corrections
%         % and the epochs.
%         %
%
%         et0      = 0.0;
%         stepsize = 1.d6;
%
%
%         for  i = 0:(NSAMP-1)
%
%            %
%            % Set the computation time for the ith sample.
%            %
%
%            et = et0 + i * stepsize;
%
%            timstr = cspice_timout( et,                                   ...
%                                    'YYYY-MON-DD HR:MN:SC.### ::TDB(TDB)' );
%
%
%            fprintf( '\n\nObservation epoch:  %s\n', timstr )
%
%
%            for  coridx = 1:NCORR
%
%               abcorr = abcorrs( coridx );
%
%               fprintf( '\n  abcorr = %s\n', char(abcorr) );
%
%               for  midx = 1:NMETHOD
%
%                  %
%                  % Select the computation method.
%                  %
%                  method = plmethods( midx );
%
%                  fprintf( '\n    Method = %s\n ', char(method) )
%
%                  %
%                  % Compute the sub-observer point using a plate
%                  % model representation of the target's surface.
%                  %
%                  [xpt, alt, plid] = cspice_subpt_pl02( handle, dladsc,   ...
%                                                        method, target,   ...
%                                                        et,     abcorr,   ...
%                                                        obsrvr        );
%
%                  %
%                  % Represent the surface point in latitudinal
%                  % coordinates.
%                  %
%                  [ xr, xlon, xlat] = cspice_reclat( xpt );
%
%                  fprintf(                                                ...
%                  '\n    Sub-observer point on plate model surface:\n' )
%                  fprintf( '      Planetocentric Longitude (deg):  %f\n', ...
%                                                   xlon * cspice_dpr() )
%                  fprintf( '      Planetocentric Latitude  (deg):  %f\n', ...
%                                                   xlat * cspice_dpr() )
%                  fprintf( '      Radius                    (km):  %f\n', ...
%                                                                    xr )
%                  fprintf( '      Observer altitude         (km):  %f\n', ...
%                                                                  alt )
%                  fprintf( '      ID of surface point plate:       %d\n', ...
%                                                                  plid )
%
%                  %
%                  % Compute the sub-observer point using an ellipsoidal
%                  % representation of the target's surface.
%                  %
%                  method = emethods( midx );
%
%                  [xpt, trgepc, srfvec] = cspice_subpnt( method, target,  ...
%                                                         et,     FIXREF,  ...
%                                                         abcorr, obsrvr );
%
%                  %
%                  % Represent the surface point in latitudinal
%                  % coordinates.
%                  %
%                  [xr, xlon, xlat] = cspice_reclat( xpt );
%
%                  fprintf(                                                ...
%                  '    Sub-observer point on ellipsoidal surface:\n' )
%                  fprintf( '      Planetocentric Longitude (deg):  %f\n', ...
%                                                    xlon * cspice_dpr() )
%                  fprintf( '      Planetocentric Latitude  (deg):  %f\n', ...
%                                                    xlat * cspice_dpr() )
%                  fprintf( '      Radius                    (km):  %f\n', ...
%                                                                     xr )
%                  fprintf( '      Observer altitude         (km):  %f\n', ...
%                                                   cspice_vnorm(srfvec) )
%
%               end
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
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, with the following variables as inputs
%
%         meta   = 'subpt_pl02_ex1.tm';
%         dsknam = 'phobos_3_3.bds';
%
%      the output was:
%
%
%      Observation epoch:  2000-JAN-01 12:00:00.000 (TDB)
%
%        abcorr = NONE
%
%          Method = Intercept
%
%          Sub-observer point on plate model surface:
%            Planetocentric Longitude (deg):  68.184433
%            Planetocentric Latitude  (deg):  -22.013938
%            Radius                    (km):  10.904558
%            Observer altitude         (km):  276700025.580291
%            ID of surface point plate:       154969
%          Sub-observer point on ellipsoidal surface:
%            Planetocentric Longitude (deg):  68.184433
%            Planetocentric Latitude  (deg):  -22.013938
%            Radius                    (km):  11.111635
%            Observer altitude         (km):  276700025.373214
%
%          Method = Ellipsoid near point
%
%          Sub-observer point on plate model surface:
%            Planetocentric Longitude (deg):  62.358996
%            Planetocentric Latitude  (deg):  -13.611694
%            Radius                    (km):  11.189797
%            Observer altitude         (km):  276700025.467226
%            ID of surface point plate:       182987
%          Sub-observer point on ellipsoidal surface:
%            Planetocentric Longitude (deg):  62.501544
%            Planetocentric Latitude  (deg):  -13.828255
%            Radius                    (km):  11.480112
%            Observer altitude         (km):  276700025.172491
%
%        abcorr = CN+S
%
%          Method = Intercept
%
%          Sub-observer point on plate model surface:
%            Planetocentric Longitude (deg):  80.397649
%            Planetocentric Latitude  (deg):  -22.012145
%            Radius                    (km):  11.102641
%            Observer altitude         (km):  276710248.420206
%            ID of surface point plate:       161027
%          Sub-observer point on ellipsoidal surface:
%            Planetocentric Longitude (deg):  80.397649
%            Planetocentric Latitude  (deg):  -22.012145
%            Radius                    (km):  10.997901
%            Observer altitude         (km):  276710248.524541
%
%          Method = Ellipsoid near point
%
%          Sub-observer point on plate model surface:
%            Planetocentric Longitude (deg):  77.596717
%            Planetocentric Latitude  (deg):  -14.326159
%            Radius                    (km):  11.278719
%            Observer altitude         (km):  276710248.357560
%            ID of surface point plate:       184540
%          Sub-observer point on ellipsoidal surface:
%            Planetocentric Longitude (deg):  77.592521
%            Planetocentric Latitude  (deg):  -14.314136
%            Radius                    (km):  11.261264
%            Observer altitude         (km):  276710248.374842
%
%
%      Observation epoch:  2000-JAN-13 01:46:40.000 (TDB)
%
%        abcorr = NONE
%
%          Method = Intercept
%
%          Sub-observer point on plate model surface:
%            Planetocentric Longitude (deg):  -28.168808
%            Planetocentric Latitude  (deg):  -23.838824
%            Radius                    (km):  12.939782
%            Observer altitude         (km):  286106844.028053
%            ID of surface point plate:       105656
%          Sub-observer point on ellipsoidal surface:
%            Planetocentric Longitude (deg):  -28.168808
%            Planetocentric Latitude  (deg):  -23.838824
%            Radius                    (km):  11.740687
%            Observer altitude         (km):  286106845.227148
%
%          Method = Ellipsoid near point
%
%          Sub-observer point on plate model surface:
%            Planetocentric Longitude (deg):  -22.941997
%            Planetocentric Latitude  (deg):  -13.930327
%            Radius                    (km):  13.832136
%            Observer altitude         (km):  286106843.393083
%            ID of surface point plate:       135034
%          Sub-observer point on ellipsoidal surface:
%            Planetocentric Longitude (deg):  -22.381480
%            Planetocentric Latitude  (deg):  -12.794856
%            Radius                    (km):  12.437687
%            Observer altitude         (km):  286106844.817034
%
%        abcorr = CN+S
%
%          Method = Intercept
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 189 lines have been
%      provided.
%
%
%-Particulars
%
%   cspice_subpt_pl02 computes the sub-observer point (abbreviated as
%   "sub-point") on a target body. cspice_subpt_pl02 also determines the
%   distance from the observer to the sub-observer point.
%
%   Sub-point Definitions
%   =====================
%
%   This routine offers two ways of defining the sub-point:
%
%      - The 'intercept' method. In general, this definition
%        calls for defining a ray emanating from the observer and
%        passing through the center of the target body. The intercept
%        on the first plate (the one closest to the observer) hit by this
%        ray is the sub-point.
%
%        The details of this definition are a bit more complex, because
%        this routine handles the case where the observer is inside the
%        target. In such cases, the sub-point is actually the point
%        that would be obtained by scaling up the target center-
%        observer vector so as to place the observer outside the target,
%        then computing the sub-point in the usual way. This handling
%        of the interior observer case prevents an observer location
%        that is slightly below the surface from accidentally being
%        associated with a sub-point on the opposite side of the target.
%        However, the possibility that the "sub-point" may be "above"
%        the observer may seem counterintuitive.
%
%      - The 'ellipsoid near point' method. When a target's surface is
%        modeled by a set of triangular plates, the notion of "dropping
%        a perpendicular segment to the surface," which makes sense
%        for convex surfaces, becomes problematic: there need not be
%        any plate whose normal vector is parallel to a segment from
%        the observer to some point on that plate, or there could be
%        more than one such plate. If such a plate exists, it might
%        be located anywhere on the visible surface---not necessarily
%        "below" the observer.
%
%        To work around these problems, the ellipsoid near point method
%        uses a reference ellipsoid to define a preliminary sub-point:
%        for an exterior observer, this is the unique point on the
%        ellipsoid's surface at which the outward surface normal points
%        toward the observer. Then the plate model sub-point is defined
%        as the plate intercept closest to the observer of a ray
%        emanating from the observer and passing through the preliminary
%        sub-point on the ellipsoid.
%
%   For a large target such as Mars, or for any target whose reference
%   ellipsoid deviates significantly from spherical, the results
%   obtained using the two sub-point definitions can be quite different.
%   The example program provided below demonstrates this fact; Phobos is
%   the target body in this case. Some analysis on the user's part will
%   be needed to select the "best" definition for a given application.
%
%   When comparing sub-point computations with results from sources
%   other than SPICE, it's essential to make sure the same geometric
%   definitions are used.
%
%
%   Aberration Corrections
%   ======================
%
%   Below, we indicate the aberration corrections to use for some
%   common applications:
%
%      1) Compute the sub-observer point using the apparent direction
%         and orientation of a target. This is the most common case for
%         a remote-sensing observation. When the observer's altitude
%         is more than one target radius above the surface:
%
%            Use 'CN+S': apply both converged light time and stellar
%            aberration corrections.
%
%         Note that when the observer is close to the target surface,
%         this choice may yield inaccurate results, since light time is
%         measured between the observer and the target center. When the
%         observer has altitude of less than one target radius above the
%         surface, aberration corrections should be omitted, so in this
%         case abcorr should be set to:
%
%            'NONE'
%
%      2) Use a geometric position vector and uncorrected target
%         orientation as low-accuracy estimates for an application where
%         execution speed is critical.
%
%            Use 'NONE'.
%
%   See the header of the CSPICE routine spkezr_c for a detailed
%   discussion of aberration corrections.
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
%   3)  If `obsrvr' and `target' map to the same NAIF integer ID codes, the
%       error SPICE(BODIESNOTDISTINCT) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If frame definition data enabling the evaluation of the state
%       of the target relative to the observer in the target
%       body-fixed frame have not been loaded prior to calling
%       cspice_subpt_pl02, an error is signaled by a routine in the call tree
%       of this routine.
%
%   5)  If the specified aberration correction is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If insufficient ephemeris data have been loaded prior to
%       calling cspice_subpt_pl02, an error is signaled by a
%       routine in the call tree of this routine.
%
%   7)  If a DSK providing a DSK type 2 plate model has not been
%       loaded prior to calling cspice_subpt_pl02, an error is signaled by a
%       routine in the call tree of this routine.
%
%   8)  If the computation method is "ellipsoid near point" and radii
%       of the target body have not been loaded into the kernel pool,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   9)  If PCK data supplying a rotation model for the target body
%       have not been loaded prior to calling cspice_subpt_pl02, an error is
%       signaled by a routine in the call tree of this routine.
%
%   10) If the segment associated with the input DLA descriptor does not
%       contain data for the designated target, the error
%       SPICE(TARGETMISMATCH) is signaled by a routine in the call tree
%       of this routine.
%
%   11) If the segment associated with the input DLA descriptor is not
%       of data type 2, the error SPICE(WRONGDATATYPE) is signaled by
%       a routine in the call tree of this routine.
%
%   12) If the sub-point cannot be computed because the ray from the
%       observer to the aim point designated by `method' fails to
%       intersect the target surface as defined by the plate model,
%       the error SPICE(NOINTERCEPT) is signaled by a routine in the
%       call tree of this routine.
%
%   13) Use of transmission-style aberration corrections is not
%       permitted. If abcorr specified such a correction, the
%       error SPICE(NOTSUPPORTED) is signaled by a routine in the call
%       tree of this routine.
%
%   14) If any of the input arguments, `handle', `dladsc', `method', `target'
%       `et', `abcorr' or `obsrvr', is undefined, an error is signaled by the
%       Matlab error handling system.
%
%   15) If any of the input arguments, `handle', `dladsc', `method', `target'
%       `et', `abcorr' or `obsrvr', is not of the expected type, or it does
%       not have the expected dimensions and size, an error is signaled by
%       the Mice interface.
%
%-Files
%
%   Appropriate DSK, SPK, PCK, and frame data must be available to
%   the calling program before this routine is called. Typically
%   the data are made available by loading kernels; however the
%   data may be supplied via subroutine interfaces if applicable.
%
%   The following data are required:
%
%   -  DSK data:  a DSK file containing a plate model representing the
%      target body's surface must be loaded. This kernel must contain
%      a type 2 segment that contains data for the entire surface of
%      the target body.
%
%   -  SPK data: ephemeris data for target and observer must be
%      loaded. If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%      Typically ephemeris data are made available by loading one
%      or more SPK files via cspice_furnsh.
%
%   -  PCK data: triaxial radii for the target body must be loaded
%      into the kernel pool if the "Ellipsoid Near Point" method is
%      selected. Typically these data are made available by loading a
%      text PCK file via cspice_furnsh.
%
%   -  Further PCK data: rotation data for the target body must
%      be loaded. These may be provided in a text or binary PCK file.
%      Either type of file may be loaded via cspice_furnsh.
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
%   1)  This routine assumes that the origin of the body-fixed reference
%       frame associated with the target body is located in the interior
%       of that body.
%
%-Required_Reading
%
%   MICE.REQ
%   DSK.REQ
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-OCT-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Changed the output
%       argument name "plateID" to "pltid" for consistency with other
%       functions.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%       Index lines now state that this routine is deprecated.
%
%   -Mice Version 1.0.0, 25-JUL-2016 (NJB) (EDW)
%
%-Index_Entries
%
%   DEPRECATED sub-observer point using plate_model
%   DEPRECATED sub-observer point using type_2 DSK
%   DEPRECATED sub-spacecraft point using plate_model
%   DEPRECATED sub-spacecraft point using type_2 DSK
%
%-&

function [spoint, alt, pltid] = cspice_subpt_pl02( handle, dladsc, method, ...
                                                   target, et,     abcorr, ...
                                                   obsrvr                )

   switch nargin
      case 7

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         method = zzmice_str( method );
         target = zzmice_str( target );
         et     = zzmice_dp( et );
         abcorr = zzmice_str( abcorr );
         obsrvr = zzmice_str( obsrvr );

      otherwise

         error ( [ 'Usage: [spoint(3), alt, pltid] = '                     ...
            'cspice_subpt_pl02( handle, dladsc(SPICE_DLA_DSCSIZ), '        ...
            '`method`, `target`, et, `abcorr`, `obsrvr` )' ] )

   end

   %
   % Call the MEX library.
   %
   try

      [spoint, alt, pltid] = mice( 'subpt_pl02',                           ...
                                   handle, dladsc, method,                 ...
                                   target, et,     abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end


