%-Abstract
%
%   CSPICE_LATSRF maps an array of planetocentric longitude/latitude
%   coordinate pairs to surface points on a specified target body.
%
%   The surface of the target body may be represented by a triaxial
%   ellipsoid or by topographic data provided by DSK files.
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
%      method      is a short string providing parameters defining
%                  the computation method to be used. In the syntax
%                  descriptions below, items delimited by brackets
%                  are optional.
%
%                  [1,c1] = size(method); char = class(method)
%
%                     or
%
%                  [1,1] = size(method); cell = class(method)
%
%                  `method' may be assigned the following values:
%
%                     'ELLIPSOID'
%
%                        The surface point computation uses a triaxial
%                        ellipsoid to model the surface of the target
%                        body. The ellipsoid's radii must be available
%                        in the kernel pool.
%
%
%                     'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                        The surface point computation uses topographic
%                        data to model the surface of the target body.
%                        These data must be provided by loaded DSK
%                        files.
%
%                        The surface list specification is optional. The
%                        syntax of the list is
%
%                           <surface 1> [, <surface 2>...]
%
%                        If present, it indicates that data only for the
%                        listed surfaces are to be used; however, data
%                        need not be available for all surfaces in the
%                        list. If absent, loaded DSK data for any surface
%                        associated with the target body are used.
%
%                        The surface list may contain surface names or
%                        surface ID codes. Names containing blanks must
%                        be delimited by double quotes, for example
%
%                           SURFACES = "Mars MEGDR 128 PIXEL/DEG"
%
%                        If multiple surfaces are specified, their names
%                        or IDs must be separated by commas.
%
%                        See the -Particulars section below for details
%                        concerning use of DSK data.
%
%
%                  Neither case nor white space are significant in
%                  `method', except within double-quoted strings. For
%                  example, the string " eLLipsoid " is valid.
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
%      target      is the name of the target body. `target' is
%                  case-insensitive, and leading and trailing blanks in
%                  `target' are not significant. Optionally, you may
%                  supply a string containing the integer ID code for
%                  the object. For example both "MOON" and "301" are
%                  legitimate strings that indicate the Moon is the
%                  target body.
%
%                  When the target body's surface is represented by a
%                  tri-axial ellipsoid, this routine assumes that a
%                  kernel variable representing the ellipsoid's radii is
%                  present in the kernel pool. Normally the kernel
%                  variable would be defined by loading a PCK file.
%
%
%      et          is the epoch for which target surface data will be
%                  selected, if the surface is modeled using DSK data.
%                  In this case, only segments having time coverage that
%                  includes the epoch `et' will be used.
%
%                  `et' is ignored if the target is modeled as an
%                  ellipsoid.
%
%                  `et' is expressed as TDB seconds past J2000 TDB.
%
%      fixref      is the name of a body-fixed reference frame centered
%                  on the target body. `fixref' may be any such frame
%                  supported by the SPICE system, including built-in
%                  frames (documented in the Frames Required Reading)
%                  and frames defined by a loaded frame kernel (FK). The
%                  string `fixref' is case-insensitive, and leading and
%                  trailing blanks in `fixref' are not significant.
%
%                  The output surface points in the array `srfpts' will be
%                  expressed relative to this reference frame.
%
%      lonlat      is an array of pairs of planetocentric longitudes and
%                  latitudes of surface points.
%
%                  [2,n] = size(lonlat); double = class(code)
%
%                  Elements
%
%                     lonlat(1,i)
%                     lonlat(2.i)
%
%                  are, respectively, the planetocentric longitude and
%                  latitude of the Ith surface point, where `i' ranges
%                  from 1 to n.
%
%                  The units of longitude and latitude are radians.
%
%   the call:
%
%      srfpts = cspice_latsrf( method, target, et, fixref, lonlat )
%
%   returns:
%
%      srfpts      is an array of target body surface points
%                  corresponding to the pairs of coordinates in the
%                  input `lonlat' array.
%
%                  [3,n] = size(srfpts); double = class(srfpts)
%
%                  Elements
%
%                     srfpts(1,i)
%                     srfpts(2,i)
%                     srfpts(3,i)
%
%                  are the Cartesian coordinates, expressed in the
%                  reference frame designated by `fixref', of the surface
%                  point corresponding to the Ith pair of input
%                  coordinates, where `i' ranges from 1 to n.
%
%                  If there are multiple solutions for a given input
%                  coordinate pair, this routine will return the point
%                  at those coordinates having the greatest distance
%                  from the origin of the coordinate system.
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
%   1) In the following example program, a DSK file containing a
%      type 2 segment is used to provide a plate model representation
%      of the surface of Phobos.
%
%      Find the surface points on a target body corresponding to a
%      given planetocentric longitude/latitude grid.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: latsrf_ex1.tm
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
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            phobos512.bds                    DSK based on
%                                             Gaskell ICQ Q=512
%                                             plate model
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'pck00010.tpc',
%                                'phobos512.bds' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function latsrf_ex1()
%
%         %
%         % Set target, reference frame, and epoch.
%         %
%         target = 'phobos';
%         fixref = 'iau_phobos';
%         et     = 0.0;
%
%         %
%         % Use both a reference ellipsoid and DSK data
%         % to represent the surface.
%         %
%         method = { 'ELLIPSOID', 'DSK/UNPRIORITIZED' };
%
%         %
%         % Load the meta-kernel.
%         %
%         cspice_furnsh( 'latsrf_ex1.tm' )
%
%         %
%         % Now generate the grid points.  We generate
%         % points along latitude bands, working from
%         % north to south.  The latitude range is selected
%         % to range from +45 to -45 degrees.  Longitude
%         % ranges from 0 to 120 degrees.  The increment
%         % is 90 degrees for latitude and 60 degrees for
%         % longitude.
%         %
%         lat = 45:-90:-45;
%         lon = 0:60:120;
%         n   = 0;
%         grid = eye(2, numel(lat) * numel(lon) );
%
%         for i=1:numel(lat)
%            for j=1:numel(lon)
%
%               n = n+1;
%
%               grid(1,n) = lon(j);
%               grid(2,n) = lat(i);
%
%            end
%         end
%
%         grid = grid * cspice_rpd();
%
%         %
%         % Find the surface points corresponding to the grid points.
%         %
%         % Compute outward normal vectors at the surface points,
%         % using both surface representations.
%         %
%         for i = 1:2
%            srfpts = cspice_latsrf( method(i), target, et, fixref, grid);
%
%            for j=1:n
%
%               %
%               % Use cspice_recrad rather than cspice_reclat to produce
%               % non-negative longitudes.
%               %
%               [ xr, xlon, xlat] = cspice_recrad( srfpts(1:3, j) );
%
%               fprintf( [ '\n%s\n'                         ...
%                          'Intercept for grid point %d:\n' ...
%                          '  Cartesian coordinates: '      ...
%                          '(%11.4e, %11.4e, %11.4e)\n'     ...
%                          '  Latitudinal Coordinates:\n'   ...
%                          '   Longitude (deg): %12.6f\n'   ...
%                          '   Latitude  (deg): %12.6f\n'   ...
%                          '   Radius     (km): %12.6f\n'   ...
%                          '\n'                             ...
%                          '  Original Grid Coordinates:\n' ...
%                          '   Longitude (deg): %12.6f\n'   ...
%                          '   Latitude  (deg): %12.6f\n'   ...
%                          '\n' ],                          ...
%                          char( method(i) ),               ...
%                          j,                               ...
%                          srfpts(1,j),   srfpts(2,j),   srfpts(3,j),    ...
%                          xlon*cspice_dpr(),   xlat*cspice_dpr(),   xr, ...
%                          grid(1,j)*cspice_dpr(), grid(2,j)*cspice_dpr() )
%
%            end
%
%         end
%
%         cspice_kclear();
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      ELLIPSOID
%      Intercept for grid point 1:
%        Cartesian coordinates: ( 7.4550e+00,  0.0000e+00,  7.4550e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):    45.000000
%         Radius     (km):    10.542977
%
%        Original Grid Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):    45.000000
%
%
%      ELLIPSOID
%      Intercept for grid point 2:
%        Cartesian coordinates: ( 3.5966e+00,  6.2296e+00,  7.1933e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):    45.000000
%         Radius     (km):    10.172847
%
%        Original Grid Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):    45.000000
%
%
%      ELLIPSOID
%      Intercept for grid point 3:
%        Cartesian coordinates: (-3.5966e+00,  6.2296e+00,  7.1933e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):    45.000000
%         Radius     (km):    10.172847
%
%        Original Grid Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):    45.000000
%
%
%      ELLIPSOID
%      Intercept for grid point 4:
%        Cartesian coordinates: ( 7.4550e+00,  0.0000e+00, -7.4550e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):   -45.000000
%         Radius     (km):    10.542977
%
%        Original Grid Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):   -45.000000
%
%
%      ELLIPSOID
%      Intercept for grid point 5:
%        Cartesian coordinates: ( 3.5966e+00,  6.2296e+00, -7.1933e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):   -45.000000
%         Radius     (km):    10.172847
%
%        Original Grid Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):   -45.000000
%
%
%      ELLIPSOID
%      Intercept for grid point 6:
%        Cartesian coordinates: (-3.5966e+00,  6.2296e+00, -7.1933e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):   -45.000000
%         Radius     (km):    10.172847
%
%        Original Grid Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):   -45.000000
%
%
%      DSK/UNPRIORITIZED
%      Intercept for grid point 1:
%        Cartesian coordinates: ( 7.1817e+00,  0.0000e+00,  7.1817e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):    45.000000
%         Radius     (km):    10.156402
%
%        Original Grid Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):    45.000000
%
%
%      DSK/UNPRIORITIZED
%      Intercept for grid point 2:
%        Cartesian coordinates: ( 3.5820e+00,  6.2042e+00,  7.1640e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):    45.000000
%         Radius     (km):    10.131412
%
%        Original Grid Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):    45.000000
%
%
%      DSK/UNPRIORITIZED
%      Intercept for grid point 3:
%        Cartesian coordinates: (-3.6854e+00,  6.3832e+00,  7.3707e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):    45.000000
%         Radius     (km):    10.423766
%
%        Original Grid Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):    45.000000
%
%
%      DSK/UNPRIORITIZED
%      Intercept for grid point 4:
%        Cartesian coordinates: ( 8.0269e+00,  0.0000e+00, -8.0269e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):   -45.000000
%         Radius     (km):    11.351730
%
%        Original Grid Coordinates:
%         Longitude (deg):     0.000000
%         Latitude  (deg):   -45.000000
%
%
%      DSK/UNPRIORITIZED
%      Intercept for grid point 5:
%        Cartesian coordinates: ( 3.3336e+00,  5.7739e+00, -6.6672e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):   -45.000000
%         Radius     (km):     9.428818
%
%        Original Grid Coordinates:
%         Longitude (deg):    60.000000
%         Latitude  (deg):   -45.000000
%
%
%      DSK/UNPRIORITIZED
%      Intercept for grid point 6:
%        Cartesian coordinates: (-3.7986e+00,  6.5793e+00, -7.5972e+00)
%        Latitudinal Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):   -45.000000
%         Radius     (km):    10.744021
%
%        Original Grid Coordinates:
%         Longitude (deg):   120.000000
%         Latitude  (deg):   -45.000000
%
%
%-Particulars
%
%   This routine is intended to be used for target body surfaces that
%   have a unique radius for each pair of planetocentric longitude
%   and latitude coordinates.
%
%   If the target surface is represented by topographic data, it is
%   possible for there to be multiple surface points at a given
%   planetocentric longitude and latitude. For example, this can
%   occur if the surface has features such as cliffs, caves, or
%   arches.
%
%   For more complex surfaces, the routine
%
%      DSKSXV {DSK, ray-surface intercept, vectorized}
%
%   may be more suitable. That routine works with rays having vertices
%   anywhere outside of the target body.
%
%
%   Planetocentric coordinates
%   ==========================
%
%   Planetocentric longitude and latitude are defined as follows:
%
%      Longitude of a point P is the angle between the prime meridian
%      and the meridian containing P. The direction of increasing
%      longitude is from the +X axis towards the +Y axis.
%
%      Latitude of a point P is the angle from the XY plane of the
%      ray from the origin through the point.
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
%      Syntax of the METHOD input argument
%      -----------------------------------
%
%      The keywords and surface list in the `method' argument
%      are called "clauses." The clauses may appear in any
%      order, for example
%
%         DSK/<surface list>/UNPRIORITIZED
%         DSK/UNPRIORITIZED/<surface list>
%         UNPRIORITIZED/<surface list>/DSK
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
%         SURFACES =
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
%-Exceptions
%
%   1)  If the target body name input string cannot be converted to an
%       integer ID code, the error SPICE(IDCODENOTFOUND) is signaled
%       by a routine in the call tree of this routine.
%
%   2)  If the input target body-fixed frame `fixref' is not recognized,
%       the error SPICE(NOFRAME) is signaled by a routine in the call
%       tree of this routine. A frame name may fail to be recognized
%       because a required frame specification kernel has not been
%       loaded; another cause is a misspelling of the frame name.
%
%   3)  If the input frame `fixref' is not centered at the target body,
%       the error SPICE(INVALIDFRAME) is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If data are not available to convert between the frame
%       `fixref' and the frame of a DSK segment of interest, an error
%       is signaled by a routine in the call tree of this
%       routine.
%
%   5)  If the input argument `method' cannot be parsed, an error
%       is signaled by either this routine or a routine in
%       the call tree of this routine.
%
%   6)  If the computation method specifies an ellipsoidal target
%       model, and if triaxial radii of the target body have not been
%       loaded into the kernel pool prior to calling cspice_latsrf, an error
%       is signaled by a routine in the call tree of this routine.
%
%   7)  If the computation method specifies an ellipsoidal target
%       model, and if any of the radii of the target body are
%       non-positive, an error is signaled by a routine in the call
%       tree of this routine. The target must be an extended body.
%
%   8)  If `method' specifies that the target surface is represented by
%       DSK data, and no DSK files are loaded for the specified
%       target, an error is signaled by a routine in the call tree
%       of this routine.
%
%   9)  If `method' specifies that the target surface is represented
%       by DSK data, and data representing the portion of the surface
%       corresponding to the coordinates provided in `lonlat' are not
%       available, an error is signaled by a routine in the call
%       tree of this routine.
%
%   10) If a surface point cannot be computed because the ray
%       corresponding to a longitude/latitude pair fails to intersect
%       the target surface as defined by the plate model, the error
%       SPICE(NOINTERCEPT) is signaled by a routine in the call tree
%       of this routine.
%
%   11) If the surface point corresponding to a longitude/latitude
%       pair in `lonlat' does not have matching longitude and latitude
%       (because it is on the opposite side of the origin), the error
%       SPICE(SHAPENOTSUPPORTED) is signaled by a routine in the call
%       tree of this routine.
%
%   12) If the radii are not available in the kernel pool, an error is
%       signaled by a routine in the call tree of this routine.
%
%   13) If the target shape is "ellipsoid" and not all radii of the
%       ellipsoid are strictly positive, the error
%       SPICE(BADAXISLENGTH) is signaled by a routine in the call tree
%       of this routine.
%
%   14) If any of the input arguments, `method', `target', `et',
%       `fixref' or `lonlat', is undefined, an error is signaled by
%       the Matlab error handling system.
%
%   15) If any of the input arguments, `method', `target', `et',
%       `fixref' or `lonlat', is not of the expected type, or it does
%       not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Appropriate kernels must be loaded by the calling program before
%   this routine is called.
%
%   The following data are required:
%
%   -  Shape data for the target body:
%
%         PCK data:
%
%            If the target shape is modeled as an ellipsoid,
%            triaxial radii for the target body must be loaded into
%            the kernel pool. Typically this is done by loading a
%            text PCK file via cspice_furnsh.
%
%         DSK data:
%
%            If the target shape is modeled by DSK data, DSK files
%            containing topographic data for the target body must be
%            loaded. If a surface list is specified, data for at
%            least one of the listed surfaces must be loaded.
%
%   -  Target body orientation data: these may be provided in a
%      text or binary PCK file. In some cases, target body
%      orientation may be provided by one more more CK files. In
%      either case, data are made available by loading the files
%      via cspice_furnsh.
%
%   The following data may be required:
%
%   -  Frame data: if a frame definition is required to convert
%      between the body-fixed frame of the target and the frame of
%      a DSK segment providing topographic data, that definition
%      must be available in the kernel pool. Typically the
%      definition is supplied by loading a frame kernel via cspice_furnsh.
%
%   -  Surface name-ID associations: if surface names are specified
%      in `method', the association of these names with their
%      corresponding surface ID codes must be established by
%      assignments of the kernel variables
%
%         NAIF_SURFACE_NAME
%         NAIF_SURFACE_CODE
%         NAIF_SURFACE_BODY
%
%      Normally these associations are made by loading a text
%      kernel containing the necessary assignments. An example of
%      such a set of assignments is
%
%         NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
%         NAIF_SURFACE_CODE += 1
%         NAIF_SURFACE_BODY += 499
%
%   -  SCLK data: if the target body's orientation is provided by
%      CK files, an associated SCLK kernel must be loaded.
%
%   In all cases, kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   1)  This routine assumes that the origin of the body-fixed
%       reference frame associated with the target body is located in
%       the interior of that body.
%
%   2)  The results returned by this routine may not be meaningful
%       if the target surface has multiple surface points associated
%       with some (longitude, latitude) coordinates.
%
%-Required_Reading
%
%   MICE.REQ
%   FRAMES.REQ
%   PCK.REQ
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
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited -Examples section to comply with NAIF standard. Added example's
%       problem statement. Reduced the number of grid points to compute in
%       code example. Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 03-MAR-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   map latitudinal coordinates to Cartesian surface points
%   map latitudinal coordinates to DSK surface points
%
%-&

function [srfpts] = cspice_latsrf( method, target, et, fixref, lonlat )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         lonlat = zzmice_dp(lonlat);

      otherwise

         error ( ['Usage: [srfpts] = cspice_latsrf( `method`, `target`, ' ...
                                                'et, `fixref`, lonlat )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [srfpts] = mice( 'latsrf_c', method, target, et, fixref, lonlat );
   catch spiceerr
      rethrow(spiceerr)
   end


