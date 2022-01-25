%-Abstract
%
%   CSPICE_SRFNRM maps an array of surface points on a specified target
%   body to the corresponding unit length outward surface normal vectors.
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
%      method   a short string providing parameters defining
%               the computation method to be used.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               In the syntax descriptions below, items delimited by brackets
%               are optional.
%
%               `method' may be assigned the following values:
%
%                  'ELLIPSOID'
%
%                     The normal vector computation uses a triaxial
%                     ellipsoid to model the surface of the target
%                     body. The ellipsoid's radii must be available
%                     in the kernel pool.
%
%                  'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                     The normal vector computation uses topographic
%                     data to model the surface of the target body.
%                     These data must be provided by loaded DSK
%                     files.
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
%                        SURFACES = "Mars MEGDR 128 PIXEL/DEG"
%
%                     If multiple surfaces are specified, their names
%                     or IDs must be separated by commas.
%
%                     See the -Particulars section below for details
%                     concerning use of DSK data.
%
%
%               Neither case nor white space are significant in
%               `method', except within double-quoted strings. For
%               example, the string ' eLLipsoid ' is valid.
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
%      target   the name of the target body.
%
%               [1,c1] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               `target' is case-insensitive, and leading and trailing
%               blanks in `target' are not significant. Optionally, you may
%               supply a string containing the integer ID code for
%               the object. For example both 'MOON' and '301' are
%               legitimate strings that indicate the Moon is the
%               target body.
%
%               When the target body's surface is represented by a
%               tri-axial ellipsoid, this routine assumes that a
%               kernel variable representing the ellipsoid's radii is
%               present in the kernel pool. Normally the kernel
%               variable would be defined by loading a PCK file.
%
%      et       the epoch for which target surface data will be
%               selected, if the surface is modeled using DSK data.
%
%               [1,1] = size(et); double = class(et)
%
%               In this case, only segments having time coverage that
%               includes the epoch `et' will be used.
%
%               `et' is ignored if the target is modeled as an
%               ellipsoid.
%
%               `et' is expressed as TDB seconds past J2000 TDB.
%
%      fixref   the name of a body-fixed reference frame centered
%               on the target body.
%
%               [1,c1] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               `fixref' may be any such frame supported by the SPICE system,
%               including built-in frames (documented in the Frames Required
%               Reading) and frames defined by a loaded frame kernel (FK). The
%               string `fixref' is case-insensitive, and leading and
%               trailing blanks in `fixref' are not significant.
%
%               The input surface points in the array `srfpts' are
%               expressed relative to this reference frame, as are
%               the normal vectors computed by this routine.
%
%      srfpts   an array of target body surface points.
%
%               [3,n] = size(srfpts); double = class(srfpts)
%
%               Elements
%
%                  srfpts(1,i)
%                  srfpts(2,i)
%                  srfpts(3,i)
%
%               are the Cartesian coordinates, expressed in the
%               reference frame designated by `fixref', of the ith
%               surface point in the array. Each surface point
%               represents an offset from the center of that frame.
%
%               All surface points must actually be "on" the surface,
%               that is, the distance of each point from the surface
%               must be less than a small margin. See the -Parameters
%               section below for a description of this margin.
%
%   the call:
%
%      normls = cspice_srfnrm( method, target, et, fixref, srfpts )
%
%   returns:
%
%      normls   an array of unit length, outward normal vectors
%               corresponding to the points in `srfpts'.
%
%               [3,n] = size(normls); double = class(normls)
%
%               Elements
%
%                  normls(1,i)
%                  normls(2,i)
%                  normls(3,i)
%
%               are the Cartesian coordinates, expressed in the
%               reference frame designated by `fixref', of the ith
%               normal vector in the array.
%
%-Parameters
%
%   SPICE_DSK_PTMEMM
%
%               is the default point-surface membership margin. This
%               margin limits the distance an input point can be from
%               a surface and still be considered to lie on that
%               surface.
%
%               The details of the application of
%
%                  SPICE_DSK_PTMEMM
%
%               are implementation-dependent. In the DSK case, roughly
%               speaking, a point-surface distance limit within a DSK
%               segment is set to
%
%                  SPICE_DSK_PTMEMM * MAXR
%
%               where MAXR is the radius of an outer bounding sphere
%               for the segment.
%
%               For shapes modeled as ellipsoids, the expression
%               above is applied to the maximum radius of the
%               ellipsoid.
%
%               See the include file
%
%                  MiceDtl.m
%
%               for the declaration of SPICE_DSK_PTMEMM. This margin
%               can be overridden. See this header file
%               and the routine cspice_dskstl for details.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Compute outward normal vectors at surface points on a target
%      body, where the points correspond to a given planetocentric
%      longitude/latitude grid. Use both ellipsoid and DSK shape
%      models.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: srfnrm_ex1.tm
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
%
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
%      function srfnrm_ex1()
%
%         %
%         % Local constants
%         %
%         META = 'srfnrm_ex1.tm';
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
%         cspice_furnsh( META )
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
%            normls = cspice_srfnrm( method(i), target, et, fixref, srfpts);
%
%            for j=1:n
%
%               %
%               % Use cspice_recrad rather than cspice_reclat to produce
%               % non-negative longitudes.
%               %
%               [ xr, xlon, xlat] = cspice_recrad( srfpts(1:3, j) );
%
%               fprintf( [ '\n%s\n'                                        ...
%                              'Surface point for grid point %d:\n'        ...
%                              '  Latitudinal Coordinates:\n'              ...
%                              '    Longitude           (deg): %12.6f\n'   ...
%                              '    Latitude            (deg): %12.6f\n'   ...
%                              '    Radius              (km): %12.6f\n' ], ...
%                              char( method(i) ),                          ...
%                              j,                                          ...
%                              xlon*cspice_dpr(),                          ...
%                              xlat*cspice_dpr(),                          ...
%                              xr )
%
%               [ nrmrad, nrmlon, nrmlat] = cspice_recrad( normls(1:3, j) );
%
%               fprintf( [ '\n'                                            ...
%                              'Normal vector direction:\n'                ...
%                              '    Longitude (deg):           %12.6f\n'   ...
%                              '    Latitude  (deg):           %12.6f\n'], ...
%                              nrmlon * cspice_dpr(),                      ...
%                              nrmlat * cspice_dpr()  )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      ELLIPSOID
%      Surface point for grid point 1:
%        Latitudinal Coordinates:
%          Longitude           (deg):     0.000000
%          Latitude            (deg):    45.000000
%          Radius              (km):    10.542977
%
%      Normal vector direction:
%          Longitude (deg):               0.000000
%          Latitude  (deg):              63.895146
%
%      ELLIPSOID
%      Surface point for grid point 2:
%        Latitudinal Coordinates:
%          Longitude           (deg):    60.000000
%          Latitude            (deg):    45.000000
%          Radius              (km):    10.172847
%
%      Normal vector direction:
%          Longitude (deg):              66.059787
%          Latitude  (deg):              58.877649
%
%      ELLIPSOID
%      Surface point for grid point 3:
%        Latitudinal Coordinates:
%          Longitude           (deg):   120.000000
%          Latitude            (deg):    45.000000
%          Radius              (km):    10.172847
%
%      Normal vector direction:
%          Longitude (deg):             113.940213
%          Latitude  (deg):              58.877649
%
%      ELLIPSOID
%      Surface point for grid point 4:
%        Latitudinal Coordinates:
%          Longitude           (deg):     0.000000
%          Latitude            (deg):   -45.000000
%          Radius              (km):    10.542977
%
%      Normal vector direction:
%          Longitude (deg):               0.000000
%          Latitude  (deg):             -63.895146
%
%      ELLIPSOID
%      Surface point for grid point 5:
%        Latitudinal Coordinates:
%          Longitude           (deg):    60.000000
%          Latitude            (deg):   -45.000000
%          Radius              (km):    10.172847
%
%      Normal vector direction:
%          Longitude (deg):              66.059787
%          Latitude  (deg):             -58.877649
%
%      ELLIPSOID
%      Surface point for grid point 6:
%        Latitudinal Coordinates:
%          Longitude           (deg):   120.000000
%          Latitude            (deg):   -45.000000
%          Radius              (km):    10.172847
%
%      Normal vector direction:
%          Longitude (deg):             113.940213
%          Latitude  (deg):             -58.877649
%
%      DSK/UNPRIORITIZED
%      Surface point for grid point 1:
%        Latitudinal Coordinates:
%          Longitude           (deg):     0.000000
%          Latitude            (deg):    45.000000
%          Radius              (km):    10.156402
%
%      Normal vector direction:
%          Longitude (deg):             341.337568
%          Latitude  (deg):              62.610726
%
%      DSK/UNPRIORITIZED
%      Surface point for grid point 2:
%        Latitudinal Coordinates:
%          Longitude           (deg):    60.000000
%          Latitude            (deg):    45.000000
%          Radius              (km):    10.131412
%
%      Normal vector direction:
%          Longitude (deg):              48.859884
%          Latitude  (deg):              56.924717
%
%      DSK/UNPRIORITIZED
%      Surface point for grid point 3:
%        Latitudinal Coordinates:
%          Longitude           (deg):   120.000000
%          Latitude            (deg):    45.000000
%          Radius              (km):    10.423766
%
%      Normal vector direction:
%          Longitude (deg):             118.553200
%          Latitude  (deg):              55.906774
%
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 132 lines have been
%      provided.
%
%
%-Particulars
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
%      -------------------------------------
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
%         'Mars MEGDR 128 PIXEL/DEG'                1
%         'Mars MEGDR 64 PIXEL/DEG'                 2
%         'Mars_MRO_HIRISE'                         3
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
%   1)  If the target body name specified in the input string cannot
%       be converted to an integer ID code, the error
%       SPICE(IDCODENOTFOUND) is signaled by a routine in the call
%       tree of this routine.
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
%       loaded into the kernel pool prior to calling cspice_srfnrm, an error
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
%   9)  If `method' specifies that the target surface is represented by
%       DSK data, and data representing the portion of the surface
%       corresponding to the surface points provided in `srfpts' are
%       not available, an error is signaled by a routine in the
%       call tree of this routine.
%
%   10) If an input surface point is not within a small tolerance of
%       the specified surface, the error SPICE(POINTNOTONSURFACE) is
%       signaled by a routine in the call tree of this routine. See
%       the -Parameters section for details.
%
%   11) If the radii are not available in the kernel pool, an error is
%       signaled by a routine in the call tree of this routine.
%
%   12) If the target shape is "ellipsoid" and not all radii of the
%       ellipsoid are strictly positive, the error
%       SPICE(BADAXISLENGTH) is signaled by a routine in the call tree
%       of this routine.
%
%   13) If any of the input arguments, `method', `target', `et',
%       `fixref' or `srfpts', is undefined, an error is signaled by
%       the Matlab error handling system.
%
%   14) If any of the input arguments, `method', `target', `et',
%       `fixref' or `srfpts', is not of the expected type, or it does
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
%        PCK data:
%
%           If the target shape is modeled as an ellipsoid,
%           triaxial radii for the target body must be loaded into
%           the kernel pool. Typically this is done by loading a
%           text PCK file via cspice_furnsh.
%
%        DSK data:
%
%           If the target shape is modeled by DSK data, DSK files
%           containing topographic data for the target body must be
%           loaded. If a surface list is specified, data for at
%           least one of the listed surfaces must be loaded.
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
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   DSK.REQ
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 10-AUG-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement. Reduced the number of grid points
%       to compute in code example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 24-MAR-2016 (EDW) (NJB)
%
%-Index_Entries
%
%   map Cartesian surface points to normal vectors
%   compute normal vectors on topographic surface
%   compute normal vectors on DSK surface
%
%-&

function [normls] = cspice_srfnrm( method, target, et, fixref, srfpts )

   switch nargin
      case 5

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         srfpts = zzmice_dp(srfpts);

      otherwise

         error ( ['Usage: [normls] =  cspice_srfnrm( `method`, `target`, ' ...
                                                'et, `fixref`, srfpts )'] )

   end

   %
   % Call the MEX library.
   %
   try
      [normls] = mice( 'srfnrm_c', method, target, et, fixref, srfpts );
   catch spiceerr
      rethrow(spiceerr)
   end


