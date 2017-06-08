%-Abstract
%
%   Deprecated: This routine has been superseded by the CSPICE routine
%   cspice_subslr. This routine is supported for purposes of backward
%   compatibility only.
%
%   CSPICE_SUBSOL_PL02 returns the rectangular coordinates of the sub-solar
%   point on a target body at a particular epoch, optionally corrected for
%   light time and stellar aberration.  The target body's surface is
%   represented by a triangular plate model contained in a type 2 DSK
%   segment. Return the sub-solar point's coordinates expressed in the
%   body-fixed frame associated with the target body.  Also, return the
%   observer's distance from the sub-solar point.
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
%                  access. This kernel must contain a type 2 segment that
%                  provides a plate model representing the entire surface
%                  of the target body.
%
%                  [1,1] = size(handle); int32 = class(handle)
%
%      dladsc      the DLA descriptor of a DSK segment representing the
%                  surface of the target body.
%
%                 [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                 int32 = class(dladsc)
%
%      method      a short string specifying the computation method
%                  to be used.
%
%                  [1,c1] = size(method); char = class(method)
%
%                     or
%
%                  [1,1] = size(method); cell = class(method)
%
%                  The choices are:
%
%                     'Intercept'        The sub-solar point is defined as
%                                        the plate model surface intercept
%                                        of the ray starting at the Sun and
%                                        passing through the target's
%                                        center.
%
%                     'Ellipsoid
%                      near point'       The sub-solar point is defined as
%                                        the plate model surface intercept
%                                        of the ray starting at the Sun and
%                                        passing through the nearest point
%                                        to the observer on a reference
%                                        ellipsoid associated with the
%                                        target body.
%
%                                        This option requires that the
%                                        reference ellipsoid's radii be
%                                        available in the kernel pool.
%
%                  Neither case nor white space are significant in the
%                  string "method". For example, the string
%
%                     '  ellipsoidNEARPOINT'
%
%                  is valid.
%
%      target      the name of the target body. `target' is
%                  case-insensitive, and leading and trailing blanks in
%                  `target' are not significant. Optionally, you may supply
%                  a string containing the integer ID code for the object.
%                  For example both 'MOON' and '301' are legitimate strings
%                  that indicate the moon is the target body.
%
%                  This routine assumes that the target body's surface is
%                  represented by a plate model, and that a DSK file
%                  containing the plate model has been loaded via dasopr_c.
%
%                  [1,c2] = size(target); char = class(target)
%
%                     or
%
%                  [1,1] = size(target); cell = class(target)
%
%      et          the epoch, represented  as seconds past J2000 TDB, at
%                  which the sub-solar point on the target body is to be
%                  computed. When aberration corrections are used, `et'
%                  refers to the epoch at which radiation is received at
%                  the observer.
%
%                  [1,1] = size(et); double = class(et)
%
%      abcorr      indicates the aberration corrections to be applied to
%                  the position and orientation of the target body and the
%                  position of the Sun to account for one-way light time
%                  and stellar aberration. See the discussion in the
%                  Particulars section for recommendations on how to choose
%                  aberration corrections.
%
%                  [1,c3] = size(abcorr); char = class(abcorr)
%
%                     or
%
%                  [1,1] = size(abcorr); cell = class(abcorr)
%
%                  `abcorr' may be any of the following:
%
%                     'NONE'     Apply no correction. Use the geometric
%                                positions of the Sun and target body
%                                relative to the observer; evaluate the
%                                target body's orientation at `et'.
%
%                  The following values of `abcorr' apply to the
%                  "reception" case in which photons depart from the
%                  target's location at the light-time corrected epoch
%                  et-lt and *arrive* at the observer's location at
%                  `et':
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                uses the position and orientation of the
%                                target at the moment it emitted photons
%                                arriving at the observer at `et'. The
%                                position of the Sun relative to the
%                                target is corrected for the one-way light
%                                time from the Sun to the target.
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation (see Particulars for details).
%                                The solution invoked by the 'LT' option
%                                uses one iteration.
%
%                     'LT+S'     Correct for one-way light time and stellar
%                                aberration using a Newtonian formulation.
%                                This option modifies the positions
%                                obtained with the 'LT' option to account
%                                for the observer's velocity relative to
%                                the solar system barycenter (note the
%                                target plays the role of "observer" in the
%                                computation of the aberration-corrected
%                                target-Sun vector). The result is the
%                                sub-solar point computed using apparent
%                                position and orientation of the target as
%                                seen by the observer and the apparent
%                                position of the Sun as seen by the target.
%
%                     'CN'       Converged Newtonian light time correction.
%                                In solving the light time equation, the
%                                'CN' correction iterates until the
%                                solution converges (three iterations on
%                                all supported platforms).
%
%                     'CN+S'     Converged Newtonian light time
%                                and stellar aberration corrections.
%
%
%      obsrvr      the name of the observing body. This is typically a
%                  spacecraft, the earth, or a surface point on the earth.
%                  `obsrvr' is case-insensitive, and leading and trailing
%                  blanks in `obsrvr' are not significant. Optionally, you
%                  may supply a string containing the integer ID code for
%                  the object. For example both 'EARTH' and '399' are
%                  legitimate strings that indicate the earth is the
%                  observer.
%
%                  [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                     or
%
%                  [1,1] = size(obsrvr); cell = class(obsrvr)
%
%   the call:
%
%      [spoint, dist, plateid] =                               ...
%                  cspice_subsol_pl02( handle, dladsc, method, ...
%                                      target, et,     abcorr, ...
%                                      obsrvr  )
%
%   returns:
%
%      spoint      the sub-solar point on the target body expressed
%                  relative to the body-fixed reference frame of the target
%                  body.
%
%                  [3,1] = size(spoint); double = class(spoint)
%
%                  The definition of sub-solar point depends on the
%                  selected computation method. See the description of the
%                  input argument `method' for details.
%
%                  The target body-fixed frame, which is time-dependent, is
%                  evaluated at `et' if `abcorr' is 'NONE'; otherwise the
%                  frame is evaluated at et-lt, where `lt' is the one-way
%                  light time from target to observer.
%
%                  The position and orientation of the target body and the
%                  position of the Sun are corrected for aberration as
%                  specified by `abcorr'; the corrected positions and
%                  orientation are used in the computation of `spoint'.
%
%      dist        the distance between the observer and the sub-solar
%                  point. The observer is presumed to be outside the
%                  target body, so `dist' is always non-negative.
%
%                  [1,1] = size(dist); double = class(dist)
%
%      plateID     the integer ID code of the plate on which the
%                  sub-solar point is located. This ID code can be
%                  use to look up data associated with the plate, such
%                  as the plate's vertices or outward normal vector.
%
%                  [1,1] = size(plateID); int32 = class(plateID)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Find the sub-solar point on Phobos as seen from Earth for a
%   specified sequence of times. Perform the computation twice,
%   using both the "intercept" and "ellipsoid near point"
%   options. Compute the corresponding sub-solar point values
%   using an ellipsoidal surface for comparison.
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
%      function subsol_pl02_t( meta, dsk )
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
%         % Now compute sub-solar points using both computation
%         % methods. We'll vary the aberration corrections
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
%            timstr = cspice_timout( et,                                    ...
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
%               fprintf( '\n   abcorr = %s\n', char(abcorr) );
%
%               for  midx = 1:NMETHOD
%
%                  %
%                  % Select the computation method.
%                  %
%                  method = plmethods( midx );
%
%                  fprintf( '\n     Method =%s\n ', char(method) )
%
%                  %
%                  % Compute the sub-solar point using a plate
%                  % model representation of the target's surface.
%                  %
%                  [xpt, dist, plid] = ...
%                        cspice_subsol_pl02( handle, dladsc, method, ...
%                                            target, et,     abcorr, ...
%                                            obsrvr                    );
%
%                  %
%                  % Represent the surface point in latitudinal
%                  % coordinates.
%                  %
%                  [ xr, xlon, xlat] = cspice_reclat( xpt );
%
%                  fprintf(                                                 ...
%                  '\n     Sub-solar point on plate model surface:\n' )
%                  fprintf( '       Planetocentric Longitude (deg):  %f\n', ...
%                                                    xlon * cspice_dpr() )
%                  fprintf( '       Planetocentric Latitude  (deg):  %f\n', ...
%                                                    xlat * cspice_dpr() )
%                  fprintf( '       Radius                    (km):  %f\n', ...
%                                                                     xr )
%                  fprintf( '       Observer distance         (km):  %f\n', ...
%                                                                   dist )
%                  fprintf( '       ID of surface point plate:       %d\n', ...
%                                                                   plid )
%
%                  %
%                  % Compute the sub-solar point using an ellipsoidal
%                  % representation of the target's surface.
%                  %
%                  method = emethods( midx );
%
%                  [xpt, trgepc, srfvec] = cspice_subslr( method, target, ...
%                                                         et,     FIXREF, ...
%                                                         abcorr, obsrvr );
%
%                  %
%                  % Represent the surface point in latitudinal
%                  % coordinates.
%                  %
%                  [xr, xlon, xlat] = cspice_reclat( xpt );
%
%                  fprintf( '     Sub-solar point on ellipsoidal surface:\n' )
%                  fprintf( '       Planetocentric Longitude (deg):  %f\n', ...
%                                                     xlon * cspice_dpr() )
%                  fprintf( '       Planetocentric Latitude  (deg):  %f\n', ...
%                                                     xlat * cspice_dpr() )
%                  fprintf( '       Radius                    (km):  %f\n', ...
%                                                                      xr )
%                  fprintf( '       Observer distance         (km):  %f\n', ...
%                                                    cspice_vnorm(srfvec) )
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
%   MATLAB outputs:
%
%      >> subsol_pl02_t( 'illum_pl02.tm' ,'phobos_3_3.bds')
%
%      Observation epoch:  2000-JAN-01 12:00:00.000 (TDB)
%
%         abcorr = NONE
%
%           Method =Intercept
%
%           Sub-solar point on plate model surface:
%             Planetocentric Longitude (deg):  102.413905
%             Planetocentric Latitude  (deg):  -24.533127
%             Radius                    (km):  11.612325
%             Observer distance         (km):  276700026.580116
%             ID of surface point plate:       164811
%           Sub-solar point on ellipsoidal surface:
%             Planetocentric Longitude (deg):  102.413905
%             Planetocentric Latitude  (deg):  -24.533127
%             Radius                    (km):  10.922580
%             Observer distance         (km):  276700027.168434
%
%           Method =Ellipsoid near point
%
%           Sub-solar point on plate model surface:
%             Planetocentric Longitude (deg):  105.857346
%             Planetocentric Latitude  (deg):  -16.270558
%             Radius                    (km):  11.645162
%             Observer distance         (km):  276700027.058857
%             ID of surface point plate:       192093
%           Sub-solar point on ellipsoidal surface:
%             Planetocentric Longitude (deg):  105.973365
%             Planetocentric Latitude  (deg):  -15.976232
%             Radius                    (km):  11.249340
%             Observer distance         (km):  276700027.400706
%
%         abcorr = CN+S
%
%           Method =Intercept
%
%           Sub-solar point on plate model surface:
%             Planetocentric Longitude (deg):  114.623420
%             Planetocentric Latitude  (deg):  -24.533628
%             Radius                    (km):  11.411417
%             Observer distance         (km):  276710249.789163
%             ID of surface point plate:       170492
%           Sub-solar point on ellipsoidal surface:
%             Planetocentric Longitude (deg):  114.623420
%             Planetocentric Latitude  (deg):  -24.533628
%             Radius                    (km):  11.046740
%             Observer distance         (km):  276710250.099485
%
%           Method =Ellipsoid near point
%
%           Sub-solar point on plate model surface:
%             Planetocentric Longitude (deg):  120.870428
%             Planetocentric Latitude  (deg):  -15.247903
%             Radius                    (km):  11.350346
%             Observer distance         (km):  276710250.680859
%             ID of surface point plate:       205719
%           Sub-solar point on ellipsoidal surface:
%             Planetocentric Longitude (deg):  120.795481
%             Planetocentric Latitude  (deg):  -15.366726
%             Radius                    (km):  11.494153
%             Observer distance         (km):  276710250.555135
%
%                ...
%
%-Particulars
%
%   cspice_subsol_pl02 computes the sub-solar point on a target body.
%   cspice_subsol_pl02 also determines the distance of the observer from the
%   sub-solar point.
%
%   Sub-point Definitions
%   =====================
%
%   This routine offers two ways of defining the sub-solar point:
%
%      - The 'intercept' method. In general, this definition
%        calls for defining a ray emanating from the Sun and
%        passing through the center of the target body. The intercept
%        on the first plate (the one closest to the observer) hit by this
%        ray is the sub-point.
%
%      - The 'ellipsoid near point' method. When a target's surface is
%        modeled by a set of triangular plates, the notion of "dropping
%        a perpendicular segment to the surface," which makes sense for
%        convex surfaces, becomes problematic:  there need not be any
%        plate whose normal vector is parallel to a segment from the Sun
%        to some point on that plate, or there could be more than one
%        such plate. If such a plate exists, it might be located
%        anywhere on the visible surface---not necessarily "below" the
%        Sun.
%
%        To work around these problems, the ellipsoid near point method
%        uses a reference ellipsoid to define a preliminary sub-solar
%        point: this is the unique point on the ellipsoid's surface at
%        which the outward surface normal points toward the Sun. Then
%        the plate model sub-solar point is defined as the plate
%        intercept closest to the Sun of a ray emanating from the Sun
%        and passing through the preliminary sub-solar point on the
%        ellipsoid.
%
%   For a large target such as Mars, or for any target whose reference
%   ellipsoid deviates significantly from spherical, the results
%   obtained using the two sub-point definitions can be quite different.
%   The example program provided below demonstrates this fact; Phobos is
%   the target body in this case. Some analysis on the user's part will
%   be needed to select the "best" definition for a given application.
%
%   When comparing sub-solar point computations with results from
%   sources other than SPICE, it's essential to make sure the same
%   geometric definitions are used.
%
%
%   Aberration Corrections
%   ======================
%
%   Below, we indicate the aberration corrections to use for some
%   common applications:
%
%      1) Compute the sub-solar point using the apparent direction
%         and orientation of a target. This is the most common case for
%         a remote-sensing observation. When the observer's altitude
%         is more than one target radius above the surface:
%
%            Use 'LT+S':  apply both light time and stellar
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
%         Note that this selection calls for using the geometric position
%         of the Sun.
%
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
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine subsol_pl02.
%
%   MICE.REQ
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
%   sub-solar point using triangular plate_model
%   sub-solar point using type_2 dsk
%
%-&

function [spoint, dist, plateID] = ...
                  cspice_subsol_pl02( handle, dladsc, method, ...
                                      target, et,     abcorr, ...
                                      obsrvr                      )

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

         error ( [ 'Usage: [spoint(3), dist, plateID] = '          ...
                            'cspice_subsol_pl02( handle, '         ...
                            'dladsc(SPICE_DLA_DSCSIZ), `method`, ' ...
                            '`target`, et, `abcorr`, `obsrvr` )' ] )

   end

   %
   % Call the MEX library.
   %
   try

      [spoint, dist, plateID] = mice( 'subsol_pl02', ...
                                      handle, dladsc, method, ...
                                      target, et,     abcorr, obsrvr);
   catch
      rethrow(lasterror)
   end


