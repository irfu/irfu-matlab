%-Abstract
%
%   CSPICE_ILLUMG computes the illumination angles (phase, incidence,
%   and emission) at a specified surface point of a target body.
%
%   The surface of the target body may be represented by a triaxial
%   ellipsoid or by topographic data provided by DSK files.
%
%   The illumination source is a specified ephemeris object
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA  INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED
%   'AS-IS' TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING
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
%      method   is a short string providing parameters defining
%               the computation method to be used. In the syntax
%               descriptions below, items delimited by brackets
%               are optional.
%
%               [1,c1] = size(method); char = class(method)
%
%                or
%
%               [1,1] = size(method); cell = class(method)
%
%               `method' may be assigned the following values:
%
%                  'ELLIPSOID'
%
%                     The illumination angle computation uses a
%                     triaxial ellipsoid to model the surface of the
%                     target body. The ellipsoid's radii must be
%                     available in the kernel pool.
%
%
%                  'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                     The illumination angle computation uses
%                     topographic data to model the surface of the
%                     target body. These data must be provided by
%                     loaded DSK files.
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
%                        'SURFACES = 'Mars MEGDR 128 PIXEL/DEG''
%
%                     If multiple surfaces are specified, their names
%                     or IDs must be separated by commas.
%
%                     See the Particulars section below for details
%                     concerning use of DSK data.
%
%
%               Neither case nor white space are significant in
%               `method', except within double-quoted strings
%               representing surface names. For example, the string
%               ' eLLipsoid ' is valid.
%
%               Within double-quoted strings representing surface names,
%               blank characters are significant, but multiple
%               consecutive blanks are considered equivalent to a single
%               blank. Case is not significant. So
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
%      target   is the name of the target body. `target' is
%               case-insensitive, and leading and trailing blanks in
%               `target' are not significant. Optionally, you may
%               supply a string containing the integer ID code for
%               the object. For example both 'MOON' and '301' are
%               legitimate strings that indicate the moon is the
%               target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%      ilusrc   the name of the illumination source. This source
%               may be any ephemeris object. Case, blanks, and
%               numeric values are treated in the same way as for the
%               input `target'.
%
%               [1,c3] = size(ilusrc); char = class(ilusrc)
%
%                  or
%
%               [1,1] = size(ilusrc); cell = class(ilusrc)
%
%      et       is the epoch, expressed as seconds past J2000 TDB,
%               for which the apparent illumination angles at the
%               specified surface point on the target body, as seen
%               from the observing body, are to be computed.
%
%               [1,1] = size(et); double = class(et)
%
%      fixref   is the name of the body-fixed, body-centered
%               reference frame associated with the target body. The
%               input surface point `spoint' and the output vector
%               `srfvec' are expressed relative to this reference
%               frame. The string `fixref' is case-insensitive, and
%               leading and trailing blanks in `fixref' are not
%               significant.
%
%               [1,c4] = size(fixref); char = class(fixref)
%
%                  or
%
%               [1,1] = size(fixref); cell = class(fixref)
%
%      abcorr   is the aberration correction to be used in computing
%               the position and orientation of the target body and
%               the location of the illumination source.
%
%               [1,c5] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               For remote sensing applications, where the apparent
%               illumination angles seen by the observer are desired,
%               normally either of the corrections
%
%                  'LT+S'
%                  'CN+S'
%
%               should be used. These and the other supported options
%               are described below. `abcorr' may be any of the
%               following:
%
%                  'NONE'     No aberration correction.
%
%               Let `lt' represent the one-way light time between the
%               observer and the input surface point `spoint' (note: NOT
%               between the observer and the target body's center). The
%               following values of `abcorr' apply to the "reception"
%               case in which photons depart from `spoint' at the
%               light-time corrected epoch et-lt and *arrive* at the
%               observer's location at `et':
%
%                  'LT'       Correct both the position of `spoint' as
%                             seen by the observer, and the position
%                             of the illumination source as seen by
%                             the target, for light time.
%
%                  'LT+S'     Correct both the position of `spoint' as
%                             seen by the observer, and the position
%                             of the illumination source as seen by
%                             the target, for light time and stellar
%                             aberration.
%
%                  'CN'       Converged Newtonian light time
%                             correction. In solving the light time
%                             equations for `spoint' and the
%                             illumination source, the 'CN'
%                             correction iterates until the solution
%                             converges.
%
%                  'CN+S'     Converged Newtonian light time and
%                             stellar aberration corrections. This
%                             option produces a solution that is at
%                             least as accurate at that obtainable
%                             with the 'LT+S' option. Whether the
%                             'CN+S' solution is substantially more
%                             accurate depends on the geometry of the
%                             participating objects and on the
%                             accuracy of the input data. In all
%                             cases this routine will execute more
%                             slowly when a converged solution is
%                             computed.
%
%               Neither case nor white space are significant in
%               `abcorr'. For example, the string
%
%                 'Lt + s'
%
%               is valid.
%
%      obsrvr   is the name of the observing body. The observing body
%               is an ephemeris object: it typically is a spacecraft,
%               the earth, or a surface point on the earth. `obsrvr' is
%               case-insensitive, and leading and trailing blanks in
%               `obsrvr' are not significant. Optionally, you may
%               supply a string containing the integer ID code for
%               the object. For example both 'MOON' and '301' are
%               legitimate strings that indicate the moon is the
%               observer.
%
%               [1,c6] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               `obsrvr' may be not be identical to `target'.
%
%      spoint   is a surface point on the target body, expressed in
%               Cartesian coordinates, relative to the body-fixed
%               target frame designated by `fixref'.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               `spoint' need not be visible from the observer's
%               location at the epoch `et'.
%
%               The components of `spoint' have units of km.
%
%   the call:
%
%      [trgepc, srfvec, phase, solar, emissn] = cspice_illumg( method,      ...
%                                               target, ilusrc, et, fixref, ...
%                                               abcorr, obsrvr, spoint
%
%   returns:
%
%      trgepc   is the "target surface point epoch." `trgepc' is defined
%               as follows: letting `lt' be the one-way light time
%               between the observer and the input surface point
%               `spoint', `trgepc' is either the epoch et-lt or `et'
%               depending on whether the requested aberration correction
%               is, respectively, for received radiation or omitted.
%               `lt' is computed using the method indicated by `abcorr'.
%
%               [1,1] = size(trgepc); double = class(trgepc)
%
%               'trgepc' is expressed as seconds past J2000 TDB.
%
%      srfvec   is the vector from the observer's position at `et' to
%               the aberration-corrected (or optionally, geometric)
%               position of `spoint', where the aberration corrections
%               are specified by `abcorr'. `srfvec' is expressed in the
%               target body-fixed reference frame designated by
%               `fixref', evaluated at `trgepc'.
%
%               [3,1] = size(srfvec); double = class(srfvec)
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
%                     xform  = cspice_pxfrm2 ( fixref, ref, trgepc, et )
%                     refvec = xform * srfvec
%
%      phase    is the phase angle at `spoint', as seen from `obsrvr' at
%               time `et'. This is the angle between the spoint-obsrvr
%               vector and the spoint-illumination source vector. Units
%               are radians. The range of `phase' is [0, pi]. See
%               Particulars below for a detailed discussion of the
%               definition.
%
%               [1,1] = size(incdnc); double = class(incdnc)
%
%      incdnc   is the illumination source incidence angle at `spoint',
%               as seen from `obsrvr' at time `et'. This is the angle
%               between the surface normal vector at `spoint' and the
%               spoint-illumination source vector. Units are radians.
%               The range of `incdnce' is [0, pi]. See Particulars below
%               for a detailed discussion of the definition.
%
%               [1,1] = size(solar); double = class(solar)
%
%      emissn   is the emission angle at `spoint', as seen from `obsrvr'
%               at time `et'. This is the angle between the surface
%               normal vector at `spoint' and the spoint-observer
%               vector. Units are radians. The range of `emissn' is [0,
%               pi]. See Particulars below for a detailed discussion of
%               the definition.
%
%               [1,1] = size(emissn); double = class(emissn)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example(1):
%
%      Find the phase, solar incidence, and emission angles at the
%      sub-solar and sub-spacecraft points on Mars as seen from the Mars
%      Global Surveyor spacecraft at a specified UTC time.
%
%      Use both an ellipsoidal Mars shape model and topographic data
%      provided by a DSK file.  For both surface points, use the 'near
%      point' and 'nadir' definitions for ellipsoidal and DSK shape
%      models, respectively.
%
%      Use converged Newtonian light time and stellar aberration
%      corrections.
%
%      The topographic model is based on data from the MGS MOLA DEM
%      megr90n000cb, which has a resolution of 4 pixels/degree. A
%      triangular plate model was produced by computing a 720 x 1440
%      grid of interpolated heights from this DEM, then tessellating
%      the height grid. The plate model is stored in a type 2 segment
%      in the referenced DSK file.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%         KPL/MK
%
%         File: illum_t2.tm
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
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0011.tls                     Leapseconds
%            mgs_ext12_ipng_mgs95j.bsp        MGS ephemeris
%            megr90n000cb_plate.bds           Plate model based on
%                                             MEGDR DEM, resolution
%                                             4 pixels/degree.
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de430.bsp',
%                                'mar097.bsp',
%                                'pck00010.tpc',
%                                'naif0011.tls',
%                                'mgs_ext12_ipng_mgs95j.bsp',
%                                'megr90n000cb_plate.bds'      )
%         \begintext
%
%
%      function illumg_t2
%
%         %
%         % Load kernel files.
%         %
%         cspice_furnsh( 'illum_t2.tm' )
%
%
%         %
%         % Convert the UTC request time string to seconds past
%         % J2000 TDB.
%         %
%         utc = '2003 OCT 13 06:00:00 UTC';
%
%         et = cspice_str2et( utc );
%
%         %
%         % Assign observer, target, and illumination source
%         % names. The acronym MGS indicates Mars Global
%         % Surveyor. See NAIF_IDS for a list of names
%         % recognized by SPICE.
%         %
%         % Also set the target body-fixed frame and
%         % the aberration correction flag.
%         %
%
%         target = 'Mars';
%         obsrvr = 'MGS';
%         ilusrc = 'Sun';
%         fixref = 'IAU_MARS';
%         abcorr = 'CN+S';
%
%         ilumth  = {'Ellipsoid', 'DSK/Unprioritized' };
%         submth =  {'Near Point/Ellipsoid', 'DSK/Nadir/Unprioritized' };
%
%         for i=1:numel(ilumth)
%
%            %
%            % Find the sub-solar point on the Earth as seen from
%            % the MGS spacecraft at et. Use the 'near point'
%            % style of sub-point definition.
%            %
%            [ssolpt, trgepc, srfvec] = ...
%                           cspice_subslr( submth(i), ...
%                                          target, et, fixref,  ...
%                                          abcorr, obsrvr );
%
%            %
%            % Now find the sub-spacecraft point.
%            %
%            [sscpt, trgepc, srfvec] = ...
%                           cspice_subpnt( submth(i), ...
%                                           target, et, fixref, ...
%                                           abcorr, obsrvr );
%
%            %
%            % Find the phase, solar incidence, and emission
%            % angles at the sub-solar point on the Earth as seen
%            % from MGS at time et.
%            %
%            [ trgepc, srfvec, sslphs, sslsol, sslemi ] = ...
%                           cspice_illumg( ilumth(i),   ...
%                                           target, 'SUN',  et,  fixref, ...
%                                           abcorr,  obsrvr,  ssolpt );
%
%            %
%            % Do the same for the sub-spacecraft point.
%            %
%            [ trgepc, srfvec, sscphs, sscsol, sscemi ] = ...
%                            cspice_illumg( ilumth(i), ...
%                                            target, 'SUN', et, fixref, ...
%                                            abcorr, obsrvr, sscpt );
%
%            %
%            % Convert the angles to degrees and write them out.
%            %
%            sslphs = sslphs * cspice_dpr;
%            sslsol = sslsol * cspice_dpr;
%            sslemi = sslemi * cspice_dpr;
%            sscphs = sscphs * cspice_dpr;
%            sscsol = sscsol * cspice_dpr;
%            sscemi = sscemi * cspice_dpr;
%
%            fprintf( [ '\n'                                            ...
%                       '   cspice_illumg method: %s\n'                 ...
%                       '   cspice_subpnt method: %s\n'                 ...
%                       '   cspice_subslr method: %s\n'                 ...
%                       '\n'                                            ...
%                       '      Illumination angles at the '             ...
%                       'sub-solar point:\n'                            ...
%                       '\n'                                            ...
%                       '      Phase angle            (deg): %15.9f\n'  ...
%                       '      Solar incidence angle  (deg): %15.9f\n'  ...
%                       '      Emission angle         (deg): %15.9f\n'],...
%                              char(ilumth(i)),   ...
%                              char(submth(i)),   ...
%                              char(submth(i)),   ...
%                              sslphs,      ...
%                              sslsol,      ...
%                              sslemi                                    );
%
%            if ( i == 0 )
%
%               fprintf( [ '        The solar incidence angle ' ...
%                          'should be 0.\n'                     ...
%                          '        The emission and phase '    ...
%                          'angles should be equal.\n' ] );
%            end
%
%            fprintf( [ '\n'                                    ...
%                       '      Illumination angles at the '  ...
%                       'sub-s/c point:\n' ...
%                       '\n'...
%                       '      Phase angle            (deg): %15.9f\n'  ...
%                       '      Solar incidence angle  (deg): %15.9f\n'  ...
%                       '      Emission angle         (deg): %15.9f\n'],...
%                       sscphs, ...
%                       sscsol, ...
%                       sscemi                                            );
%
%
%            if ( i == 0 )
%
%              fprintf( [ '        The emission angle '  ...
%                         'should be 0.\n'               ...
%                         '        The solar incidence ' ...
%                         'and phase angles should be equal.\n' ] );
%              end
%
%              fprintf ( '\n' );
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%   MATLAB outputs:
%
%      cspice_illumg method: Ellipsoid
%      cspice_subpnt method: Near Point/Ellipsoid
%      cspice_subslr method: Near Point/Ellipsoid
%
%         Illumination angles at the sub-solar point:
%
%         Phase angle            (deg):   138.370270685
%         Solar incidence angle  (deg):     0.000000000
%         Emission angle         (deg):   138.370270685
%
%         Illumination angles at the sub-s/c point:
%
%         Phase angle            (deg):   101.439331040
%         Solar incidence angle  (deg):   101.439331041
%         Emission angle         (deg):     0.000000002
%
%
%      cspice_illumg method: DSK/Unprioritized
%      cspice_subpnt method: DSK/Nadir/Unprioritized
%      cspice_subslr method: DSK/Nadir/Unprioritized
%
%         Illumination angles at the sub-solar point:
%
%         Phase angle            (deg):   138.387071678
%         Solar incidence angle  (deg):     0.967122745
%         Emission angle         (deg):   137.621480599
%
%         Illumination angles at the sub-s/c point:
%
%         Phase angle            (deg):   101.439331359
%         Solar incidence angle  (deg):   101.555993667
%         Emission angle         (deg):     0.117861156
%
%-Particulars
%
%   Mice contains four routines that compute illumination angles:
%
%      cspice_illumf (same as this routine, except that illumination
%                    and visibility flags are returned)
%
%      cspice_illumg (this routine)
%
%      cspice_ilumin (same as cspice_illumg, except that the sun is fixed
%                    as the illumination source)
%
%      cspice_illum  (deprecated)
%
%   cspice_illumf is the most capable of the set.
%
%
%   Illumination angles
%   ===================
%
%   The term "illumination angles" refers to the following set of
%   angles:
%
%
%      phase angle              Angle between the vectors from the
%                               surface point to the observer and
%                               from the surface point to the
%                               illumination source.
%
%      incidence angle          Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               illumination source.
%
%      emission angle           Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               observer.
%
%   The diagram below illustrates the geometric relationships
%   defining these angles. The labels for the incidence, emission,
%   and phase angles are "inc.", "e.", and "phase".
%
%
%                                                    *
%                                            illumination source
%
%                  surface normal vector
%                            ._                 _.
%                            |\                 /|  illumination
%                              \    phase      /    source vector
%                               \   .    .    /
%                               .            .
%                                 \   ___   /
%                            .     \/     \/
%                                  _\ inc./
%                           .    /   \   /
%                           .   |  e. \ /
%       *             <--------------- *  surface point on
%    viewing            vector            target body
%    location           to viewing
%    (observer)         location
%
%
%   Note that if the target-observer vector, the target normal vector
%   at the surface point, and the target-illumination source vector
%   are coplanar, then phase is the sum of the incidence and emission
%   angles. This rarely occurs; usually
%
%      phase angle  <  incidence angle + emission angle
%
%   All of the above angles can be computed using light time
%   corrections, light time and stellar aberration corrections, or no
%   aberration corrections. In order to describe apparent geometry as
%   observed by a remote sensing instrument, both light time and
%   stellar aberration corrections should be used.
%
%   The way aberration corrections are applied by this routine
%   is described below.
%
%      Light time corrections
%      ======================
%
%         Observer-target surface point vector
%         ------------------------------------
%
%         Let `et' be the epoch at which an observation or remote
%         sensing measurement is made, and let et-lt ("lt" stands
%         for "light time") be the epoch at which the photons
%         received at `et' were emitted from the surface point `spoint'.
%         Note that the light time between the surface point and
%         observer will generally differ from the light time between
%         the target body's center and the observer.
%
%
%         Target body's orientation
%         -------------------------
%
%         Using the definitions of `et' and `lt' above, the target body's
%         orientation at et-lt is used. The surface normal is
%         dependent on the target body's orientation, so the body's
%         orientation model must be evaluated for the correct epoch.
%
%
%         Target body -- illumination source vector
%         -----------------------------------------
%
%         The surface features on the target body near `spoint' will
%         appear in a measurement made at `et' as they were at et-lt.
%         In particular, lighting on the target body is dependent on
%         the apparent location of the illumination source as seen
%         from the target body at et-lt. So, a second light time
%         correction is used to compute the position of the
%         illumination source relative to the surface point.
%
%
%      Stellar aberration corrections
%      ==============================
%
%      Stellar aberration corrections are applied only if
%      light time corrections are applied as well.
%
%         Observer-target surface point body vector
%         -----------------------------------------
%
%         When stellar aberration correction is performed, the
%         direction vector `srfvec' is adjusted so as to point to the
%         apparent position of `spoint': considering `spoint' to be an
%         ephemeris object, `srfvec' points from the observer's
%         position at `et' to the light time and stellar aberration
%         corrected position of `spoint'.
%
%         Target body-illumination source vector
%         --------------------------------------
%
%         The target body-illumination source vector is the apparent
%         position of the illumination source, corrected for light
%         time and stellar aberration, as seen from the target body
%         at time et-lt.
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
%         'SURFACES ='
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
%         "SURFACES = 1"
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
%         'DSK/UNPRIORITIZED/SURFACES = "Mars MEGDR 64 PIXEL/DEG", 3'
%
%
%      Aberration corrections using DSK data
%      -------------------------------------
%
%      For irregularly shaped target bodies, the distance between the
%      observer and the nearest surface intercept need not be a
%      continuous function of time; hence the one-way light time
%      between the intercept and the observer may be discontinuous as
%      well. In such cases, the computed light time, which is found
%      using an iterative algorithm, may converge slowly or not at all.
%      In all cases, the light time computation will terminate, but
%      the result may be less accurate than expected.
%
%   Please refer to the Aberation Corrections Required Reading (ABCORR.REQ)
%   for detailed information describing the nature and calculation of the
%   applied corrections.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the Mice routine cspice_illumg.
%
%   MICE.REQ
%   ABCORR.REQ
%   DSK.REQ
%   FRAMES.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 04-APR-2017, EDW (JPL)
%
%-Index_Entries
%
%   illumination angles general source
%   lighting angles general source
%   phase angle general source
%   incidence angle general source
%   emission angle general source
%
%-&

function [trgepc, srfvec, phase, incdnc, emissn] = cspice_illumg( method, ...
                                              target, ilusrc, et, fixref, ...
                                              abcorr, obsrvr, spoint)

   switch nargin
      case 8

         method = zzmice_str(method);
         target = zzmice_str(target);
         ilusrc = zzmice_str(ilusrc);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         spoint = zzmice_dp(spoint);

      otherwise

         error ( [ 'Usage: [trgepc, srfvec(3), phase, incdnc, emissn] = ' ...
                     'cspice_ilumin( `method`, `target`, `ilusrc`, et, '  ...
                     '`fixref`, `abcorr`, `obsrvr`, spoint(3))' ] )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [ilumin] = mice('illumg_s', method, target, ilusrc, et, ...
                                  fixref, abcorr, obsrvr, spoint);
      trgepc   = reshape( [ilumin(:).trgepc], 1, [] );
      srfvec   = reshape( [ilumin(:).srfvec], 3, [] );
      phase    = reshape( [ilumin(:).phase ], 1, [] );
      incdnc   = reshape( [ilumin(:).incdnc], 1, [] );
      emissn   = reshape( [ilumin(:).emissn], 1, [] );
   catch
      rethrow(lasterror)
   end

