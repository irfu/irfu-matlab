%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   cspice_illumf. This routine is supported for purposes of backward
%   compatibility only.
%
%   CSPICE_ILLUM_PLID_PL02 computes the illumination angles---phase, 
%   solar incidence, and emission---at a specified point on a target 
%   body at a particular epoch, optionally corrected for light time 
%   and stellar aberration. In addition, it returns logical flags indicating 
%   whether the surface point is shadowed or occulted by the target body.
%
%   The target body's surface is represented by a triangular plate model
%   contained in a type 2 DSK segment. The ID of the plate on which the
%   point is located must be provided by the caller.
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
%               This kernel must contain a type 2 segment that provides a
%               plate model representing the entire surface of the target
%               body.
%
%      dladsc   the DLA descriptor of a DSK segment representing
%               the surface of a target body.
%
%               [SPICE_DLA_DSCSIZ,1] = size(dladsc); int32 = class(dladsc)
%
%      target   the name of the target body.
%
%               [1,c1] = size(target); char = class(target)
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
%               represented using a plate model, and that a DSK file
%               containing the plate model has been loaded via cspice_dasopr.
%
%      et       the epoch, represented as seconds past J2000 TDB, at
%               which the illumination angles are to be computed.
%
%               [1,1] = size(et); double = class(et)
%
%               When aberration corrections are used, `et' refers to the
%               epoch at which radiation is received at the observer.
%
%      abcorr   indicates the aberration corrections to be applied to
%               the position and orientation of the target body and the
%               position of the Sun to account for one-way light time
%               and stellar aberration.
%
%               [1,c2] = size(abcorr); char = class(abcorr)
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
%                             target-Sun vector). The result is that the
%                             illumination angles are computed using
%                             apparent position and orientation of the
%                             target as seen by the observer and the
%                             apparent position of the Sun as seen by
%                             the target.
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
%               [1,c3] = size(obsrvr); char = class(obsrvr)
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
%      spoint   a surface point on the target body, expressed in
%               rectangular body-fixed (body equator and prime meridian)
%               coordinates.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               `spoint' need not be visible from the observer's location at
%               time `et'.
%
%      plid     the integer ID of the plate on which `spoint' is
%               located.
%
%               [1,1] = size(plid); int32 = class(plid)
%
%               If `spoint' was found by calling any of the routines
%
%                  cspice_dskx02
%                  cspice_subpt_pl02
%                  cspice_subsol_pl02
%
%               `plid' is the plate ID returned by the called routine.
%
%   the call:
%
%      [trgepc, srfvec, phase, solar, emissn, visibl, lit] =               ...
%                cspice_illum_plid_pl02( handle, dladsc, target,           ...
%                                        et,     abcorr, obsrvr,           ...
%                                        spoint, plid )
%
%   returns:
%
%      trgepc   the "surface point epoch."
%
%               [1,1] = size(trgepc); double = class(trgepc)
%
%               `trgepc' is defined as follows: letting `lt' be the one-way
%               light time between the observer and the input surface point
%               `spoint', `trgepc' is either the epoch et-lt or `et' depending
%               on whether the requested aberration correction is,
%               respectively, for received radiation or omitted. `lt' is
%               computed using the method indicated by `abcorr'.
%
%               `trgepc' is expressed as seconds past J2000 TDB.
%
%      srfvec   the vector from the observer's position at `et' to
%               the aberration-corrected (or optionally, geometric)
%               position of `spoint', where the aberration corrections
%               are specified by `abcorr'.
%
%               [3,1] = size(phase); double = class(phase)
%
%               `srfvec' is expressed in the target body-fixed reference frame
%               designated by `fixref', evaluated at `trgepc'.
%
%               The components of `srfvec' are given in units of km.
%
%               One can use the function norm to obtain the
%               distance between the observer and `spoint':
%
%                  dist = norm( srfvec );
%
%               The observer's position `obspos', relative to the
%               target body's center, where the center's position is
%               corrected for aberration effects as indicated by
%               `abcorr', can be computed with:
%
%                  obspos = spoint - srfvec
%
%               To transform the vector `srfvec' from a reference frame
%               `fixref' at time `trgepc' to a time-dependent reference
%               frame `ref' at time `et', the routine cspice_pxfrm2 should be
%               called. Let `xform' be the 3x3 matrix representing the
%               rotation from the reference frame `fixref' at time
%               `trgepc' to the reference frame `ref' at time `et'. Then
%               `srfvec' can be transformed to the result `refvec' as
%               follows:
%
%                  xform  = cspice_pxfrm2 ( fixref, ref, trgepc, et )
%                  refvec = xform * srfvec
%
%      phase    the phase angle at `spoint', as seen from `obsrvr' at
%               time `et'.
%
%               [1,1] = size(phase); double = class(phase)
%
%               This is the angle between the spoint-obsrvr vector and the
%               spoint-sun vector. Units are radians. The range of `phase'
%               is [0, pi].
%
%               See -Particulars below for a detailed discussion of the
%               definitions of this angle.
%
%      solar    the solar incidence angle at `spoint', as seen from
%               `obsrvr' at time `et'.
%
%               [1,1] = size(solar); double = class(solar)
%
%               This is the angle between the surface normal vector at
%               `spoint' and the spoint-sun vector. Units are radians.
%               The range of `solar' is [0, pi].
%
%               Note that if the target surface is non-convex, a solar
%               incidence angle less than pi/2 radians does not imply
%               the surface point is illuminated. See the description of
%               `lit' below.
%
%               See -Particulars below for a detailed discussion of the
%               definitions of this angle.
%
%      emissn   the emission angle at `spoint', as seen from `obsrvr'
%               at time `et'.
%
%               [1,1] = size(emissn); double = class(emissn)
%
%               This is the angle between the surface normal vector at
%               `spoint' and the spoint-obsrvr vector. Units are radians.
%               The range of `emissn' is is [0, pi].
%
%               Note that if the target surface is non-convex, an emission
%               angle less than pi/2 radians does not imply the surface
%               point is visible to the observer. See the description of
%               `visibl' below.
%
%               See -Particulars below for a detailed discussion of the
%               definitions of this angle.
%
%      visibl   a logical flag indicating whether the surface point
%               is visible to the observer.
%
%               [1,1] = size(visibl); logical = class(visibl)
%
%               `visibl' takes into account whether the target surface occults
%               `spoint', regardless of the emission angle at `spoint'.
%               `visibl' is returned with the value true if `spoint' is
%               visible; otherwise it is false.
%
%      lit      a logical flag indicating whether the surface point
%               is illuminated; the point is considered to be
%               illuminated if the vector from the point to the center
%               of the sun doesn't intersect the target surface.
%
%               [1,1] = size(lit); logical = class(lit)
%
%               `lit' takes into account whether the target surface casts a
%               shadow on `spoint', regardless of the solar incidence
%               angle at `spoint'. `lit' is returned with the value
%               true if `spoint' is illuminated; otherwise it is
%               false.
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
%   1) Find the illumination angles at both the sub-observer point and
%      sub-solar point on Phobos as seen from Mars for a specified
%      sequence of times. Perform each computation twice, using both the
%      "intercept" and "ellipsoid near point" options for the sub-observer
%      point and sub-solar point computations. Compute the corresponding
%      illumination angles using an ellipsoidal surface for comparison.
%      (Note that the surface points on the plate model generally will
%      not lie on the ellipsoid's surface, so the emission and solar
%      incidence angles won't generally be zero at the sub-observer
%      and sub-solar points, respectively.)
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: illum_plid_pl02_ex1.tm
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
%      function illum_plid_pl02_ex1
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see MiceUser.m and
%         % MiceDSK.m.
%         %
%         MiceUser
%
%         %
%         % Constants
%         %
%         NCORR       = 2;
%         NSAMP       = 3;
%         NMETHOD     = 2;
%         FIXREF      = 'IAU_PHOBOS';
%         ILUM_METHOD = 'ELLIPSOID';
%         TOL         = 1.d-12;
%
%         %
%         % Initial values
%         %
%         abcorrs     = { 'NONE', 'CN+S' };
%         methods     = { 'Intercept', 'Ellipsoid near point' };
%
%         obsrvr      = 'Mars';
%         target      = 'Phobos';
%
%         %
%         % Prompt for the name of a meta-kernel specifying
%         % all of the other kernels we need. Load the
%         % meta kernel.
%         %
%         meta = input( 'Enter meta-kernel name > ','s');
%         cspice_furnsh( meta )
%
%         %
%         % Prompt for the name of a DSK file.
%         %
%         dsknam = input( 'Enter DSK name         > ','s');
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
%            cspice_kclear
%            fprintf( 'No segments found in DSK file %s\n', dsknam )
%            return
%
%         end
%
%         %
%         % If we made it this far, `dladsc' is the
%         % DLA descriptor of the first segment.
%         %
%
%         %
%         % Get the DSK descriptor of the segment; from this
%         % descriptor we can obtain the ID of body-fixed frame
%         % associated with the segment. We'll need this frame
%         % later to compute illumination angles on the target
%         % body's reference ellipsoid.
%         %
%         dskdsc = cspice_dskgd( handle, dladsc );
%
%         fixref = cspice_frmnam( dskdsc(SPICE_DSK_FRMIDX) );
%
%         if ( strcmp(fixref, ' ') )
%
%            cspice_kclear
%            fprintf( ['Frame ID code # could not be mapped to '           ...
%                      'a frame name %d\n'], dskdsc(SPICE_DSK_FRMIDX) )
%            return
%
%         end
%
%
%         %
%         % Now compute sub-points using both computation
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
%               fprintf( '\n   abcorr = %s\n', char(abcorr) );
%
%               for  midx = 1:NMETHOD
%
%                  %
%                  % Select the computation method.
%                  %
%                  method = methods( midx );
%
%                  fprintf( '\n     Method =%s\n ', char(method) )
%
%                  %
%                  % Compute the sub-observer point using a plate
%                  % model representation of the target's surface.
%                  %
%                  [xpt, alt, plid] = ...
%                        cspice_subpt_pl02( handle, dladsc, method,        ...
%                                           target, et,     abcorr,        ...
%                                           obsrvr                    );
%
%                  %
%                  % Compute the illumination angles at the sub-observer
%                  % point. Also compute the light-of-sight visibility and
%                  % shadowing flags.
%                  %
%                  [trgepc, srfvec, phase,  solar,                         ...
%                           emissn, visibl, lit]   =                       ...
%                       cspice_illum_plid_pl02( handle, dladsc, target,    ...
%                                               et,     abcorr, obsrvr,    ...
%                                               xpt,    plid          );
%
%                  %
%                  %  Represent the surface point in latitudinal
%                  % coordinates.
%                  %
%                  [ xr, xlon, xlat] = cspice_reclat( xpt );
%
%                  fprintf(                                                ...
%                  '\n     Sub-observer point on plate model surface:\n' )
%                  fprintf(                                                ...
%                       '       Planetocentric Longitude (deg):  %f\n',    ...
%                                                 xlon * cspice_dpr() )
%                  fprintf(                                                ...
%                       '       Planetocentric Latitude  (deg):  %f\n',    ...
%                                                   xlat * cspice_dpr() )
%
%                  fprintf(                                                ...
%                   '\n         Illumination angles derived using a\n' )
%                  fprintf( '         plate model surface:\n' )
%                  fprintf(['             Phase angle'                     ...
%                           '              (deg): %f\n'],                  ...
%                                                phase  * cspice_dpr() )
%                  fprintf(['             Solar incidence angle',          ...
%                           '    (deg): %f\n'],  solar  * cspice_dpr() )
%                  fprintf(['             Illumination flag',              ...
%                           '             : %d\n'],  lit               )
%                  fprintf(['             Emission angle',                 ...
%                           '           (deg): %f\n'],                     ...
%                                                emissn * cspice_dpr() )
%                  fprintf(['             Visibility flag',                ...
%                           '               : %d\n'], visibl           )
%                  fprintf(['             Range to surface point',         ...
%                           '    (km): %f\n'],  norm( srfvec )         )
%
%                  %
%                  % Compute the illumination angles using an ellipsoidal
%                  % representation of the target's surface. The role of
%                  % this representation is to provide an outward surface
%                  % normal.
%                  %
%
%                  [trgepc, srfvec, phase,  solar,  emissn] =              ...
%                                           cspice_ilumin( ILUM_METHOD,    ...
%                                               target, et,     FIXREF,    ...
%                                               abcorr, obsrvr, xpt);
%
%                  fprintf(                                                ...
%                   '         Illumination angles derived using an\n' )
%                  fprintf( '         ellipsoidal reference surface:\n' )
%                  fprintf(['             Phase angle',                    ...
%                           '              (deg): %f\n'],                  ...
%                                                phase  * cspice_dpr() )
%                  fprintf(['             Solar incidence angle',          ...
%                           '    (deg): %f\n'],  solar  * cspice_dpr() )
%                  fprintf(['             Emission angle',                 ...
%                           '           (deg): %f\n'],                     ...
%                                                emissn * cspice_dpr() )
%
%
%                  %
%                  % Now repeat our computations using the
%                  % sub-solar point.
%                  %
%                  % Compute the sub-solar point using a plate model
%                  % representation of the target's surface.
%                  %
%
%                  [xpt, dist, plid] =                                     ...
%                           cspice_subsol_pl02( handle, dladsc, method,    ...
%                                               target, et,     abcorr,    ...
%                                               obsrvr                );
%
%                  %
%                  % Compute the illumination angles at the sub-solar point.
%                  % Also compute the light-of-sight visibility and
%                  % shadowing flags.
%                  %
%                  [trgepc, srfvec, phase,  solar,                         ...
%                           emissn, visibl, lit] =                         ...
%                       cspice_illum_plid_pl02( handle, dladsc, target,    ...
%                                               et,     abcorr, obsrvr,    ...
%                                               xpt,    plid          );
%
%                  %
%                  %  Represent the surface point in latitudinal
%                  % coordinates.
%                  %
%                  [ xr, xlon, xlat] = cspice_reclat( xpt );
%
%                  fprintf(                                                ...
%                   '\n     Sub-solar point on plate model surface:\n' )
%                  fprintf(                                                ...
%                       '       Planetocentric Longitude (deg):  %f\n',    ...
%                                                 xlon * cspice_dpr() )
%                  fprintf(                                                ...
%                       '       Planetocentric Latitude  (deg):  %f\n',    ...
%                                                   xlat * cspice_dpr() )
%
%                  fprintf(                                                ...
%                   '\n         Illumination angles derived using a\n' )
%                  fprintf( '         plate model surface:\n' )
%                  fprintf(['             Phase angle'                     ...
%                           '              (deg): %f\n'],                  ...
%                                                phase  * cspice_dpr() )
%                  fprintf(['             Solar incidence angle',          ...
%                           '    (deg): %f\n'],  solar  * cspice_dpr() )
%                  fprintf(['             Illumination flag',              ...
%                           '             : %d\n'],  lit               )
%                  fprintf(['             Emission angle',                 ...
%                           '           (deg): %f\n'],                     ...
%                                                emissn * cspice_dpr() )
%                  fprintf(['             Visibility flag',                ...
%                           '               : %d\n'], visibl           )
%                  fprintf(['             Range to surface point',         ...
%                           '    (km): %f\n'],  norm( srfvec )         )
%
%
%                  %
%                  % Compute the illumination angles using an ellipsoidal
%                  % representation of the target's surface. The role of
%                  % this representation is to provide an outward surface
%                  % normal.
%                  %
%                  [ etrgep, esrfvc, phase, solar, emissn ] =              ...
%                          cspice_ilumin( ILUM_METHOD, target, et, fixref, ...
%                                         abcorr, obsrvr, xpt );
%
%                  fprintf(                                                ...
%                   '         Illumination angles derived using an\n' )
%                  fprintf( '         ellipsoidal reference surface:\n' )
%                  fprintf(['             Phase angle',                    ...
%                           '              (deg): %f\n'],                  ...
%                                                phase  * cspice_dpr() )
%                  fprintf(['             Solar incidence angle',          ...
%                           '    (deg): %f\n'],  solar  * cspice_dpr() )
%                  fprintf(['             Emission angle',                 ...
%                           '           (deg): %f\n'],                     ...
%                                                emissn * cspice_dpr() )
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
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the meta-kernel file named illum_plid_pl02_ex1.tm and
%      the DSK file named phobos_3_3.bds, the output was:
%
%
%      Enter meta-kernel name > illum_plid_pl02_ex1.tm
%      Enter DSK name         > phobos_3_3.bds
%
%
%      Observation epoch:  2000-JAN-01 12:00:00.000 (TDB)
%
%         abcorr = NONE
%
%           Method =Intercept
%
%           Sub-observer point on plate model surface:
%             Planetocentric Longitude (deg):  -0.348118
%             Planetocentric Latitude  (deg):  0.008861
%
%               Illumination angles derived using a
%               plate model surface:
%                   Phase angle              (deg): 101.596824
%                   Solar incidence angle    (deg): 98.376877
%                   Illumination flag             : 0
%                   Emission angle           (deg): 9.812914
%                   Visibility flag               : 1
%                   Range to surface point    (km): 9501.835727
%               Illumination angles derived using an
%               ellipsoidal reference surface:
%                   Phase angle              (deg): 101.596824
%                   Solar incidence angle    (deg): 101.695444
%                   Emission angle           (deg): 0.104977
%
%           Sub-solar point on plate model surface:
%             Planetocentric Longitude (deg):  102.413905
%             Planetocentric Latitude  (deg):  -24.533127
%
%               Illumination angles derived using a
%               plate model surface:
%                   Phase angle              (deg): 101.665306
%                   Solar incidence angle    (deg): 13.068798
%                   Illumination flag             : 1
%                   Emission angle           (deg): 98.408735
%                   Visibility flag               : 0
%                   Range to surface point    (km): 9516.720964
%               Illumination angles derived using an
%               ellipsoidal reference surface:
%                   Phase angle              (deg): 101.665306
%                   Solar incidence angle    (deg): 11.594741
%                   Emission angle           (deg): 98.125499
%
%           Method =Ellipsoid near point
%
%           Sub-observer point on plate model surface:
%             Planetocentric Longitude (deg):  -0.264850
%             Planetocentric Latitude  (deg):  0.004180
%
%               Illumination angles derived using a
%               plate model surface:
%                   Phase angle              (deg): 101.596926
%                   Solar incidence angle    (deg): 98.376877
%                   Illumination flag             : 0
%                   Emission angle           (deg): 9.812985
%                   Visibility flag               : 1
%                   Range to surface point    (km): 9501.837763
%               Illumination angles derived using an
%               ellipsoidal reference surface:
%                   Phase angle              (deg): 101.596926
%                   Solar incidence angle    (deg): 101.593324
%                   Emission angle           (deg): 0.003834
%
%           Sub-solar point on plate model surface:
%             Planetocentric Longitude (deg):  105.857346
%             Planetocentric Latitude  (deg):  -16.270558
%
%               Illumination angles derived using a
%               plate model surface:
%                   Phase angle              (deg): 101.663675
%                   Solar incidence angle    (deg): 16.476730
%                   Illumination flag             : 1
%                   Emission angle           (deg): 118.124981
%                   Visibility flag               : 0
%                   Range to surface point    (km): 9517.506732
%               Illumination angles derived using an
%               ellipsoidal reference surface:
%                   Phase angle              (deg): 101.663675
%                   Solar incidence angle    (deg): 0.422781
%                   Emission angle           (deg): 101.541470
%
%         abcorr = CN+S
%
%           Method =Intercept
%
%           Sub-observer point on plate model surface:
%             Planetocentric Longitude (deg):  -0.348101
%             Planetocentric Latitude  (deg):  0.008861
%
%               Illumination angles derived using a
%               plate model surface:
%                   Phase angle              (deg): 101.592246
%                   Solar incidence angle    (deg): 98.372348
%                   Illumination flag             : 0
%                   Emission angle           (deg): 9.812902
%                   Visibility flag               : 1
%                   Range to surface point    (km): 9502.655917
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 479 lines have been
%      provided.
%
%
%-Particulars
%
%   The term "illumination angles" refers to following set of
%   angles:
%
%
%      solar incidence angle    Angle between the surface normal at the
%                               specified surface point and the vector
%                               from the surface point to the Sun.
%
%      emission angle           Angle between the surface normal at the
%                               specified surface point and the vector
%                               from the surface point to the observer.
%
%      phase angle              Angle between the vectors from the
%                               surface point to the observing body and
%                               from the surface point to the Sun.
%
%
%   The diagram below illustrates the geometric relationships defining
%   these angles. The labels for the solar incidence, emission, and
%   phase angles are "s.i.", "e.", and "phase".
%
%
%                                                    *
%                                                   Sun
%
%                  surface normal vector
%                            ._                 _.
%                            |\                 /|  Sun vector
%                              \    phase      /
%                               \   .    .    /
%                               .            .
%                                 \   ___   /
%                            .     \/     \/
%                                  _\ s.i./
%                           .    /   \   /
%                           .   |  e. \ /
%       *             <--------------- *  surface point on
%    viewing            vector            target body
%    location           to viewing
%    (observer)         location
%
%
%   Note that if the target-observer vector, the target normal vector
%   at the surface point, and the target-sun vector are coplanar, then
%   phase is the sum of incidence and emission. This is rarely true;
%   usually
%
%      phase angle  <  solar incidence angle + emission angle
%
%
%   All of the above angles can be computed using light time
%   corrections, light time and stellar aberration corrections, or
%   no aberration corrections. The way aberration corrections
%   are used is described below.
%
%   Care must be used in computing light time corrections. The
%   guiding principle used here is "describe what appears in
%   an image."
%
%
%      Observer-target body surface point vector
%      -----------------------------------------
%
%      Let `et' be the epoch at which an observation or remote
%      sensing measurement is made, and let et - lt ("lt" stands
%      for "light time") be the epoch at which the photons received
%      at `et' were emitted from the body (we use the term "emitted"
%      loosely here).
%
%      The correct observer-target vector points from the observer's
%      location at `et' to the surface point location at et - lt.
%      The target-observer vector points in the opposite direction.
%
%      Since light time corrections are not anti-symmetric, the correct
%      target-observer vector CANNOT be found by negating the light
%      time corrected position of the observer as seen from the
%      target body.
%
%
%      Target body's orientation
%      -------------------------
%
%      Using the definitions of `et' and `lt' above, the target
%      body's orientation at et-lt is used. The surface
%      normal is dependent on the target body's orientation, so
%      the body's orientation model must be evaluated for the correct
%      epoch.
%
%
%      Target body -- Sun vector
%      -------------------------
%
%      All surface features on the target body will appear in a
%      measurement made at `et' as they were at the target at epoch
%      et-lt. In particular, lighting on the target body is dependent
%      on the apparent location of the Sun as seen from the target body
%      at et-lt. So, a second light time correction is used in finding
%      the apparent location of the Sun.
%
%
%   Stellar aberration corrections, when used, are applied as follows:
%
%
%      Observer-target body vector
%      ---------------------------
%
%      In addition to light time correction, stellar aberration is used
%      in computing the apparent target surface point position as seen
%      from the observer's location at time `et'. This apparent position
%      defines the observer-target surface point vector.
%
%
%      Target body-Sun vector
%      ----------------------
%
%      The target body-Sun vector is the apparent position of the Sun,
%      corrected for light time and stellar aberration, as seen from
%      the target body at time et-lt. Note that the target body's
%      position is not affected by the stellar aberration correction
%      applied in finding its apparent position as seen by the
%      observer.
%
%   Once all of the vectors, as well as the target body's orientation,
%   have been computed with the proper aberration corrections, the
%   element of time is eliminated from the computation. The problem
%   becomes a purely geometric one and is described by the diagram above.
%
%-Exceptions
%
%   If any of the listed errors occur, the output arguments are
%   left unchanged.
%
%   1)  If `plid' is not a valid plate ID, an error is signaled
%       by a routine in the call tree of this routine.
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
%       cspice_illum_plid_pl02, an error is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the specified aberration correction is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If insufficient ephemeris data have been loaded prior to
%       calling cspice_illum_plid_pl02, an error is signaled by a
%       routine in the call tree of this routine.
%
%   7)  If a DSK providing a DSK type 2 plate model has not been
%       loaded prior to calling cspice_illum_plid_pl02, an error is signaled
%       by a routine in the call tree of this routine.
%
%   8)  If PCK data supplying a rotation model for the target body
%       have not been loaded prior to calling cspice_illum_plid_pl02, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   9)  If the segment associated with the input DLA descriptor does not
%       contain data for the designated target, the error
%       SPICE(TARGETMISMATCH) is signaled by a routine in the call tree
%       of this routine. The target body of the DSK segment is determined
%       from the `center' member of the segment's DSK descriptor.
%
%   10) If the segment associated with the input DLA descriptor is not
%       of data type 2, the error SPICE(WRONGDATATYPE) is signaled by a
%       routine in the call tree of this routine.
%
%   11) Use of transmission-style aberration corrections is not
%       permitted. If abcorr specified such a correction, the
%       error SPICE(NOTSUPPORTED) is signaled by a routine in the call
%       tree of this routine.
%
%   12) The observer is presumed to be outside the target body; no
%       checks are made to verify this.
%
%   13) If the DSK segment's coordinate system is not latitudinal
%       (aka planetocentric), the error SPICE(BADCOORDSYSTEM) is signaled
%       by a routine in the call tree of this routine.
%
%   14) If any of the input arguments, `handle', `dladsc', `target', `et',
%       `abcorr', `obsrvr', `spoint' or `plit', is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   15) If any of the input arguments, `handle', `dladsc', `target', `et',
%       `abcorr', `obsrvr', `spoint' or `plid', is not of the expected type,
%       or it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%   -  SPK data: ephemeris data for target, observer, and Sun must be
%      loaded. If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data. Typically
%      ephemeris data are made available by loading one or more SPK
%      files via cspice_furnsh.
%
%   -  PCK data: rotation data for the target body must
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
%   1)  The solar illumination state indicated by the output argument `lit'
%       is computed treating the sun as a point light source. Surface
%       points that are illuminated by part of the sun's disc are
%       classified as "lit" or not depending on whether the center of the
%       sun is visible from those points.
%
%-Required_Reading
%
%   MICE.REQ
%   ABCORR.REQ
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
%   -Mice Version 1.1.0, 05-DEC-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
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
%   -Mice Version 1.0.0, 04-APR-2017 (NJB) (EDW)
%
%-Index_Entries
%
%   DEPRECATED plate model point visibility and shadowing
%   DEPRECATED illumination angles using DSK type 2
%   DEPRECATED lighting angles using DSK type 2 plate model
%   DEPRECATED phase angle using DSK type 2 plate model
%   DEPRECATED emission angle using DSK type 2 plate model
%   DEPRECATED solar incidence angle using DSK type 2
%
%-&

function [trgepc, srfvec, phase, solar, emissn, visibl, lit] =             ...
                cspice_illum_plid_pl02( handle, dladsc, target,            ...
                                        et,     abcorr, obsrvr,            ...
                                        spoint, plid                )

   switch nargin
      case 8

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         target = zzmice_str( target );
         et     = zzmice_dp( et );
         abcorr = zzmice_str( abcorr );
         obsrvr = zzmice_str( obsrvr );
         spoint = zzmice_dp( spoint );
         plid   = zzmice_int( plid );

      otherwise

         error ( [ 'Usage: [trgepc, srfvec(3), phase, solar, emissn, '     ...
                   'visibl, lit] = cspice_illum_plid_pl02( handle, '       ...
                   ' dladsc(SPICE_DLA_DSCSIZ), `target`, '                 ...
                   'et, `abcorr`, `obsrvr`, spoint(3), plid )' ] )

   end

   %
   % Call the MEX library.illum_plid_pl02_s
   %
   try

      [ilumin, visibl, lit]  = mice( 'illum_plid_pl02_s',                  ...
                                     handle, dladsc, target,               ...
                                     et,     abcorr, obsrvr, spoint, plid);
      trgepc   = reshape( [ilumin(:).trgepc], 1, [] );
      srfvec   = reshape( [ilumin(:).srfvec], 3, [] );
      phase    = reshape( [ilumin(:).phase ], 1, [] );
      solar    = reshape( [ilumin(:).incdnc], 1, [] );
      emissn   = reshape( [ilumin(:).emissn], 1, [] );
      visibl   = zzmice_logical(visibl);
      lit      = zzmice_logical(lit);
   catch spiceerr
      rethrow(spiceerr)
   end



