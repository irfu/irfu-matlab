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
%                        The normal vector computation uses a triaxial
%                        ellipsoid to model the surface of the target
%                        body. The ellipsoid's radii must be available
%                        in the kernel pool.
%
%
%                     'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                        The normal vector computation uses topographic
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
%                        See the Particulars section below for details
%                        concerning use of DSK data.
%
%
%                  Neither case nor white space are significant in
%                  `method', except within double-quoted strings. For
%                  example, the string ' eLLipsoid ' is valid.
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
%
%      target      is the name of the target body. `target' is
%                  case-insensitive, and leading and trailing blanks in
%                  `target' are not significant. Optionally, you may
%                  supply a string containing the integer ID code for
%                  the object. For example both 'MOON' and '301' are
%                  legitimate strings that indicate the Moon is the
%                  target body.
%
%                  [1,c1] = size(target); char = class(target)
%
%                     or
%
%                  [1,1] = size(target); cell = class(target)
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
%                  [1,1] = size(et); double = class(et)
%
%                  `et' is ignored if the target is modeled as an
%                  ellipsoid.
%
%                  `et' is expressed as TDB seconds past J2000 TDB.
%
%
%      fixref      is the name of a body-fixed reference frame centered
%                  on the target body. `fixref' may be any such frame
%                  supported by the SPICE system, including built-in
%                  frames (documented in the Frames Required Reading)
%                  and frames defined by a loaded frame kernel (FK). The
%                  string `fixref' is case-insensitive, and leading and
%                  trailing blanks in `fixref' are not significant.
%
%                  [1,c1] = size(target); char = class(target)
%
%                     or
%
%                  [1,1] = size(target); cell = class(target)
%
%                  The input surface points in the array `srfpts' are
%                  expressed relative to this reference frame, as are
%                  the normal vectors computed by this routine.
%
%      srfpts      is an array of target body surface points.
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
%                  reference frame designated by `fixref', of the ith
%                  surface point in the array. Each surface point
%                  represents an offset from the center of that frame.
%
%                  All surface points must actually be "on" the surface,
%                  that is, the distance of each point from the surface
%                  must be less than a small margin. See the Parameters
%                  section below for a description of this margin.
%
%   the call:
%
%      normls = cspice_srfnrm( method, target, et, fixref, srfpts )
%
%   returns:
%
%      normls      is an array of unit length, outward normal vectors
%                  corresponding to the points in `srfpts'.
%
%                  [3,n] = size(normls); double = class(normls)
%
%                  Elements
%
%                     normls(1,i)
%                     normls(2,i)
%                     normls(3,i)
%
%                  are the Cartesian coordinates, expressed in the
%                  reference frame designated by `fixref', of the ith
%                  normal vector in the array.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      Compute outward normal vectors at surface points on a target
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
%         File: srfnrm_t.tm
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
%            PATH_SYMBOLS    = 'GEN'
%            PATH_VALUES     = '/ftp/pub/naif/generic_kernels'
%
%            KERNELS_TO_LOAD = ( '$GEN/pck/pck00010.tpc',
%                                '$GEN/dsk/phobos/phobos512.bds' )
%         \begintext
%
%   Example(1):
%
%      function srfnrm_t
%
%         %
%         % Local constants
%         %
%         META = 'srfnrm_t.tm';
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
%         % ranges from 0 to 300 degrees.  The increment
%         % is 45 degrees for latitude and 60 degrees for
%         % longitude.
%         %
%         lat = 45:-45:-45;
%         lon = 0:60:300;
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
%               fprintf( [ '\n%s\n'                                       ...
%                              'Surface point for grid point %d:\n'       ...
%                              '  Latitudinal Coordinates:\n'             ...
%                              '    Longitude           (deg): %12.6f\n'  ...
%                              '    Latitude            (deg): %12.6f\n'  ...
%                              '    Radius              (km): %12.6f\n' ],...
%                              char( method(i) ), ...
%                              j,                 ...
%                              xlon*cspice_dpr(), ...
%                              xlat*cspice_dpr(), ...
%                              xr )
%
%               [ nrmrad, nrmlon, nrmlat] = cspice_recrad( normls(1:3, j) );
%
%               fprintf( [ '\n'                                           ...
%                              'Normal vector direction:\n'               ...
%                              '    Longitude (deg):           %12.6f\n'  ...
%                              '    Latitude  (deg):           %12.6f\n'],...
%                              nrmlon * cspice_dpr(), ...
%                              nrmlat * cspice_dpr()  )
%
%            end
%
%         end
%
%         cspice_kclear();
%
%   Matlab outputs:
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
%      ...
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
%      ...
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
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine srfnrm_c.
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
%   -Mice Version 1.0.0, 24-MAR-2016, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   map Cartesian surface points to normal vectors
%   compute normal vectors on topographic surface
%   compute normal vectors on dsk surface
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
   catch
      rethrow(lasterror)
   end


