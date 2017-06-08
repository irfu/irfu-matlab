%-Abstract
%
%   CSPICE_DSKX02 determines the plate ID and body-fixed coordinates
%   of the intersection of a specified ray with the surface defined by a
%   type 2 DSK plate model.
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
%      handle     the handle of a DSK file containing a type 2
%                 segment from which data are to be fetched.
%
%                 [1,1] = size(handle); int32 = class(handle)
%
%      dladsc     the DLA descriptor associated with the segment
%                 from which data are to be fetched.
%
%                 [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                 int32 = class(dladsc)
%
%      vertex     is the vertex of a ray.  `vertex' is expressed relative
%                 to the body fixed reference frame associated with the
%                 target body.  This reference frame is the same frame
%                 relative to which the vertices of the plate model are
%                 expressed.  Units are km.
%
%                 [3,1] = size(vertex); double = class(vertex)
%
%                 The vertex is required to be outside the target body.
%
%      raydir     is the ray's direction vector.  `raydir' is expressed
%                 relative to the body fixed reference frame associated
%                 with the target body.
%
%                 [3,1] = size(raydir); double = class(raydir)
%
%   the call:
%
%      [ plid, xpt, found] = cspice_dskx02( handle, dladsc, vertex, raydir)
%
%   returns:
%
%      plid       is the ID of the plate closest to the input ray's
%                 vertex at which a ray-surface intercept exists.
%                 If no intercept exists, `plid' is undefined.
%
%                 [1,1] = size(plid); int32 = class(plid)
%
%      xpt        is the ray-target intercept closest to the ray's vertex,
%                 if an intercept exists. `xpt' is expressed relative to
%                 the body-fixed reference frame associated with the target
%                 body.  Units are km.
%
%                 [3,1] = size(xpt); double = class(xpt)
%
%                 If no intercept exists, `xpt' is undefined.
%
%      found      is a logical flag that indicates whether or not the ray
%                 does indeed intersect the target.  If the ray intersects a
%                 plate, `found' is true.  Otherwise `found' is
%                 false.
%
%                 [1,1] = size(found); logical = class(found)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   Example (1):
%
%      Look up all the vertices associated with each plate 
%      of the model contained in a specified type 2 segment. 
%      For this example, we'll show the context of this look-up: 
%      opening the DSK file for read access, traversing a trivial, 
%      one-segment list to obtain the segment of interest. 
%
%      function dskx02_t( dsk )
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see DSKMiceUser.m and
%         % DSKMice02.m.
%         %
%         MiceUser
%
%         NLAT             = 9;
%         NLON             = 9;
%         TOL              = 1.e-12;
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
%            % contains no segments.  This is
%            % unexpected, but we're prepared for it.
%            %
%            txt = sprintf( ['SPICE(NODATA): ' ...
%                  'No segments found in DSK file %s.\n'], dsk );
%            error( txt )
%
%         end
%
%         %
%         % If we made it this far, DLADSC is the
%         % DLA descriptor of the first segment.
%         %
%         % We're going to generate the intercept points
%         % using a set of rays which point toward the
%         % origin and whose vertices are on a
%         % specified lat/lon grid.  To start out we
%         % must pick a reasonable range from the origin
%         % for the vertices:  the range must be large
%         % enough so that the vertices are guaranteed
%         % to be outside the target body but small
%         % enough that we don't lose too much precision
%         % in the surface intercept computation.
%         %
%         % We'll look up the upper bound for the target
%         % radius, then use 2 times this value as the
%         % vertex magnitude.
%         %
%         dskdsc = cspice_dskgd( handle, dladsc );
%
%         maxr =  dskdsc( SPICE_DSK_MX3IDX);
%         r    = 2.0 * maxr;
%
%         %
%         % Now generate the intercept points.  We generate
%         % intercepts along latitude bounds, working from
%         % north to south.  Latitude ranges
%         % from +80 to -80 degrees.  Longitude
%         % ranges from 0 to 320 degrees.  The increment
%         % is 20 degrees for latitude and 40 degrees for
%         % longitude.
%         %
%         for i = 0:(NLAT-1)
%
%            lat = cspice_rpd() * ( 80.0 - 20.0*i );
%
%            for j = 0:(NLON-1)
%
%               lon = cspice_rpd() * 40.0*j;
%
%               %
%               % Produce a ray vertex for the current
%               % lat/lon value.  Negate the vertex to
%               % produce the ray's direction vector.
%               %
%               vertex = cspice_latrec ( r, lon, lat );
%               raydir = -vertex;
%
%               %
%               % Find the surface intercept for this
%               % ray.
%               %
%               [plid, xpt, found] = cspice_dskx02( handle, ...
%                                      dladsc, vertex, raydir );
%
%               %
%               % Since the ray passes through the origin on
%               % the body-fixed frame associated with the
%               % target body, we'd rarely expect to find that
%               % the ray failed to intersect the surface.
%               % For safety, we check the FOUND flag.  (A
%               % "not found" condition could be a sign of
%               % a bug.)
%               %
%               if ~found
%
%                  fprintf ( ['\n'                     ...
%                           'Intercept not found!\n'   ...
%                           '   Ray vertex:\n'         ...
%                           '   Longitude (deg): %f\n' ...
%                           '   Latitude  (deg): %f\n' ...
%                           '   Radius     (km): %e\n' ...
%                           '\n'],                     ...
%                           lon * cspice_dpr(),        ...
%                           lat * cspice_dpr(),        ...
%                           r                           )
%
%               else
%
%                  %
%                  % This is the normal case.  Display the
%                  % intercept plate ID and the intercept
%                  % point in both cartesian and latitudinal
%                  % coordinates.  Show the corresponding ray
%                  % vertex to facilitate validation of results.
%                  %
%                  % Use cspice_recrad rather than cspice_reclat to produce
%                  % non-negative longitudes.
%                  %
%                  [xr, xlon, xlat] = cspice_recrad( xpt );
%
%                  fprintf ( ['\n'                                      ...
%                           'Intercept found:\n'                        ...
%                           '   Plate ID:                 %ld\n'        ...
%                           '   Cartesian Coordinates:    (%e %e %e)\n' ...
%                           '   Latitudinal Coordinates:\n'    ...
%                           '   Longitude (deg): %f\n'         ...
%                           '   Latitude  (deg): %f\n'         ...
%                           '   Radius     (km): %e\n'         ...
%                           '\n'                               ...
%                           '   Ray vertex:\n'                 ...
%                           '   Longitude (deg): %f\n'         ...
%                           '   Latitude  (deg): %f\n'         ...
%                           '   Radius     (km): %e\n'         ...
%                           '\n'],                             ...
%                           plid,                              ...
%                           xpt(1), xpt(2), xpt(3),            ...
%                           xlon * cspice_dpr(),               ...
%                           xlat * cspice_dpr(),               ...
%                           xr,                                ...
%                           lon  * cspice_dpr(),               ...
%                           lat  * cspice_dpr(),               ...
%                           r                                   )
%
%                  %
%                  % Perform sanity checks on the intercept
%                  % coordinates.  Stop the program if any error
%                  % is larger than our tolerance value.
%                  %
%                  if ( abs(xlat-lat) > TOL )
%                     disp( 'Latitude error!' )
%                  end
%
%                  if (  (xlon - lon)  > cspice_pi()  )
%                     xlon = xlon - cspice_twopi();
%                  end
%
%                  if (  (xlon - lon)  > TOL  )
%                     disp( 'Longitude error!' )
%                  end
%
%                  if ( xr  > (1.0+TOL)*maxr  )
%                     disp( 'Radius error!' )
%                  end
%
%               end
%
%               %
%               % End of longitude loop.
%               %
%
%            end
%
%            %
%            % End of latitude loop.
%            %
%
%         end
%
%         %
%         % Close the kernel.  This isn't necessary in a stand-
%         % alone program, but it's good practice in subroutines
%         % because it frees program and system resources.
%         %
%         cspice_dascls( handle );
%
%   Matlab outputs:
%      
%      dskx02_t( '/kernels/gen/dsk/phobos_3_3.bds' )
%
%      Intercept found:
%         Plate ID:                 306238
%         Cartesian Coordinates:    (1.520878e+00 0.000000e+00 8.625327e+00)
%         Latitudinal Coordinates:
%         Longitude (deg): 0.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 8.758387e+00
%      
%         Ray vertex:
%         Longitude (deg): 0.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%      
%      
%      Intercept found:
%         Plate ID:                 317112
%         Cartesian Coordinates:    (1.189704e+00 9.982799e-01 8.807772e+00)
%         Latitudinal Coordinates:
%         Longitude (deg): 40.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 8.943646e+00
%      
%         Ray vertex:
%         Longitude (deg): 40.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%      
%      
%      Intercept found:
%         Plate ID:                 324141
%         Cartesian Coordinates:    (2.777752e-01 1.575341e+00 9.072029e+00)
%         Latitudinal Coordinates:
%         Longitude (deg): 80.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.211980e+00
%      
%         Ray vertex:
%         Longitude (deg): 80.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%      
%         ...
%      
%      Intercept found:
%         Plate ID:                 1893
%         Cartesian Coordinates:    (-7.783952e-01 -1.348220e+00 -8.828997e+00)
%         Latitudinal Coordinates:
%         Longitude (deg): 240.000000
%         Latitude  (deg): -80.000000
%         Radius     (km): 8.965198e+00
%      
%         Ray vertex:
%         Longitude (deg): 240.000000
%         Latitude  (deg): -80.000000
%         Radius     (km): 2.802354e+01
%      
%      
%      Intercept found:
%         Plate ID:                 2467
%         Cartesian Coordinates:    (2.666822e-01 -1.512430e+00 -8.709738e+00)
%         Latitudinal Coordinates:
%         Longitude (deg): 280.000000
%         Latitude  (deg): -80.000000
%         Radius     (km): 8.844099e+00
%      
%         Ray vertex:
%         Longitude (deg): 280.000000
%         Latitude  (deg): -80.000000
%         Radius     (km): 2.802354e+01
%      
%      
%      Intercept found:
%         Plate ID:                 3187
%         Cartesian Coordinates:    (1.132409e+00 -9.502036e-01 -8.383597e+00)
%         Latitudinal Coordinates:
%         Longitude (deg): 320.000000
%         Latitude  (deg): -80.000000
%         Radius     (km): 8.512928e+00
%      
%         Ray vertex:
%         Longitude (deg): 320.000000
%         Latitude  (deg): -80.000000
%         Radius     (km): 2.802354e+01
%      
%-Particulars
%
%   This routine solves the ray-surface intercept problem for a
%   specified ray and a surface represented by triangular plate model.
%   The surface representation is provided by data in a type 2 segment
%   of a DSK file.
%
%   This routine does not assume that the segment from which the surface
%   model data are read represents the entire surface of the target
%   body.  A program could call this routine repeatedly to find the
%   surface intercept of a ray and a shape model partitioned into
%   multiple segments.
%
%   In general, this routine should be expected to run faster when used
%   with smaller shape models.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine dskx02_c.
%
%   MICE.REQ
%   DAS.REQ
%   DSK.REQ
%
%-Version
%
%   -Mice Version 1.0.0, 30-MAR-2014, EDW (JPL), NJB (JPL)
%
%-Index_Entries
%
%   plate and plate model point intersected by ray
%   intersection of ray and surface
%
%-&

function [ plid, xpt, found] = cspice_dskx02( handle, dladsc, vertex, raydir)

   switch nargin
      case 4

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         vertex = zzmice_dp( vertex );
         raydir = zzmice_dp( raydir );

      otherwise

         error ( [ 'Usage: [ plid, xpt, found] = ' ...
                   'cspice_dskx02( handle, dladsc, vertex, raydir)' ] )
   end

   %
   % Call the MEX library.
   %
   try

      [ plid, xpt, found] = mice( 'dskx02_c', ...
                                  handle, dladsc, vertex, raydir );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch
      rethrow(lasterror)
   end


