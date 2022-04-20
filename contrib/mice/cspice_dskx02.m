%-Abstract
%
%   CSPICE_DSKX02 determines the plate ID and body-fixed coordinates of the
%   intersection of a specified ray with the surface defined by a
%   type 2 DSK plate model.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
%   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
%   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
%   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
%   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
%   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
%   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
%   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
%   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
%
%   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
%   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
%   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
%   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
%   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
%   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
%
%   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
%   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
%   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
%   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
%
%-I/O
%
%   Given:
%
%      handle   the file handle of a DSK file containing a shape model for a
%               target body.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               The shape model is stored in a type 2 DSK segment.
%
%      dladsc   the DLA descriptor of a type 2 DSK segment containing plate
%               model data representing the surface of the target body.
%
%               [SPICE_DLA_DSCSIZ,1]  = size(dladsc);
%                               int32 = class(dladsc)
%
%               Normally this descriptor will be obtained by a search
%               through a DSK file using the DLA search routines; see the
%               -Examples header section below for a working code example
%               illustrating a simple search.
%
%      vertex   the vertex of a ray.
%
%               [3,1] = size(vertex); double = class(vertex)
%
%               `vertex' is expressed relative to the body fixed reference
%               frame associated with the target body. This reference frame
%               is the same frame relative to which the vertices of the plate
%               model are expressed. Units are km.
%
%               The vertex is required to be outside the target body.
%
%      raydir   the ray's direction vector.
%
%               [3,1] = size(raydir); double = class(raydir)
%
%               `raydir' is expressed relative to the body fixed reference
%               frame associated with the target body.
%
%   the call:
%
%      [plid, xpt, found] = cspice_dskx02( handle, dladsc, vertex, raydir )
%
%   returns:
%
%      plid     the ID of the plate closest to the input ray's
%               vertex at which a ray-surface intercept exists.
%
%               [1,1] = size(plid); int32 = class(plid)
%
%
%               If no intercept exists, `plid' is undefined.
%
%      xpt      the ray-target intercept closest to the ray's vertex,
%               if an intercept exists.
%
%               [3,1] = size(xpt); double = class(xpt)
%
%               `xpt' is expressed relative to the body-fixed reference frame
%               associated with the target body. Units are km.
%
%               If no intercept exists, `xpt' is undefined.
%
%      found    a logical flag that indicates whether or not the ray
%               does indeed intersect the target.
%
%               [1,1] = size(found); logical = class(found)
%
%               If the ray intersects a plate, `found' is true. Otherwise
%               `found' is false.
%
%-Parameters
%
%   See the parameter definitions file
%
%      MiceDtl.m
%
%   for the values of tolerance parameters used by default by the
%   ray-surface intercept algorithm.
%
%   See the parameter definitions file
%
%      MiceDLA.m
%
%   for declarations of DLA descriptor sizes and documentation of the
%   contents of DLA descriptors.
%
%   See the parameter definitions file
%
%      MiceDSK.m
%
%   for declarations of DSK descriptor sizes and documentation of the
%   contents of DSK descriptors.
%
%   See the parameter definitions file
%
%      MiceDSK.m
%
%   for declarations of DSK data type 2 (plate model) parameters.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Find the surface intercept points corresponding to a latitude/
%      longitude grid of a specified resolution, for a specified
%      target body. This simple program assumes the shape model for
%      the target body is stored in a single type 2 DSK segment, and
%      that this segment is the first one in the DSK file to which it
%      belongs.
%
%      Example code begins here.
%
%
%      function dskx02_ex1()
%
%         %
%         % MiceUser globally defines DSK parameters.
%         % For more information, please see MiceDSK.m.
%         %
%         MiceUser
%
%         NLAT             = 9;
%         NLON             = 9;
%         TOL              = 1.e-12;
%
%         %
%         % Prompt for the name of the file to search.
%         %
%         dsk = input( 'Name of DSK file > ', 's' );
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
%            txt = sprintf( ['SPICE(NODATA): '                             ...
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
%         maxr =  dskdsc(SPICE_DSK_MX3IDX);
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
%               [plid, xpt, found] = cspice_dskx02( handle, dladsc,        ...
%                                                   vertex, raydir );
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
%                  fprintf ( ['\n'                                         ...
%                           'Intercept not found!\n'                       ...
%                           '   Ray vertex:\n'                             ...
%                           '   Longitude (deg): %f\n'                     ...
%                           '   Latitude  (deg): %f\n'                     ...
%                           '   Radius     (km): %e\n'                     ...
%                           '\n'],                                         ...
%                           lon * cspice_dpr(),                            ...
%                           lat * cspice_dpr(),                            ...
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
%                  fprintf ( ['\n'                                         ...
%                           'Intercept found:\n'                           ...
%                           '   Plate ID:                 %ld\n'           ...
%                           '   Cartesian Coordinates: '                   ...
%                           '%12.8f %12.8f %12.8f\n'                       ...
%                           '   Latitudinal Coordinates:\n'                ...
%                           '   Longitude (deg): %f\n'                     ...
%                           '   Latitude  (deg): %f\n'                     ...
%                           '   Radius     (km): %e\n'                     ...
%                           '\n'                                           ...
%                           '   Ray vertex:\n'                             ...
%                           '   Longitude (deg): %f\n'                     ...
%                           '   Latitude  (deg): %f\n'                     ...
%                           '   Radius     (km): %e\n'                     ...
%                           '\n'],                                         ...
%                           plid,                                          ...
%                           xpt(1), xpt(2), xpt(3),                        ...
%                           xlon * cspice_dpr(),                           ...
%                           xlat * cspice_dpr(),                           ...
%                           xr,                                            ...
%                           lon  * cspice_dpr(),                           ...
%                           lat  * cspice_dpr(),                           ...
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
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using the DSK file named phobos_3_3.bds, the output
%      was:
%
%
%      Name of DSK file > phobos_3_3.bds
%
%      Intercept found:
%         Plate ID:                 306238
%         Cartesian Coordinates:   1.52087789   0.00000000   8.62532711
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
%         Cartesian Coordinates:   1.18970365   0.99827989   8.80777185
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
%         Cartesian Coordinates:   0.27777518   1.57534131   9.07202903
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
%
%      Intercept found:
%         Plate ID:                 327994
%         Cartesian Coordinates:  -0.81082405   1.40438846   9.19682344
%         Latitudinal Coordinates:
%         Longitude (deg): 120.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.338699e+00
%
%         Ray vertex:
%         Longitude (deg): 120.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%
%
%      Intercept found:
%         Plate ID:                 329431
%         Cartesian Coordinates:  -1.47820193   0.53802150   8.92132122
%         Latitudinal Coordinates:
%         Longitude (deg): 160.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.058947e+00
%
%         Ray vertex:
%         Longitude (deg): 160.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%
%
%      Intercept found:
%         Plate ID:                 196042
%         Cartesian Coordinates:  -1.49854761  -0.54542673   9.04411256
%         Latitudinal Coordinates:
%         Longitude (deg): 200.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.183633e+00
%
%         Ray vertex:
%         Longitude (deg): 200.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%
%
%      Intercept found:
%         Plate ID:                 235899
%         Cartesian Coordinates:  -0.78240454  -1.35516441   8.87447325
%         Latitudinal Coordinates:
%         Longitude (deg): 240.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.011376e+00
%
%         Ray vertex:
%         Longitude (deg): 240.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 2.802354e+01
%
%
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 1135 lines have
%      been provided.
%
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
%   body. A program could call this routine repeatedly to find the
%   surface intercept of a ray and a shape model partitioned into
%   multiple segments.
%
%   In general, this routine should be expected to run faster when used
%   with smaller shape models.
%
%-Exceptions
%
%   1)  If the input handle is invalid, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If a file read error occurs, the error is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If the input DLA descriptor is invalid, the effect of this
%       routine is undefined. The error *may* be diagnosed by
%       routines in the call tree of this routine, but there are no
%       guarantees.
%
%   4)  If an error occurs while trying to look up any component
%       of the shape model, the error is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If the input ray direction is the zero vector, the error
%       SPICE(ZEROVECTOR) is signaled by a routine in the call tree of
%       this routine.
%
%   6)  If the coarse voxel grid scale of the shape model is less than
%       1, the error SPICE(VALUEOUTOFRANGE) is signaled by a routine
%       in the call tree of this routine.
%
%   7)  If the coarse voxel grid of the shape model contains more than
%       SPICE_DSK02_MAXCGR (see MiceDSK.m) voxels, the error
%       SPICE(GRIDTOOLARGE) is signaled by a routine in the call tree
%       of this routine.
%
%   8)  If the plate list for any intersected voxel is too large
%       for this routine to buffer, the error SPICE(ARRAYTOOSMALL)
%       is signaled by a routine in the call tree of this routine.
%
%   9)  Due to round-off errors, results from this routine may
%       differ across platforms. Results also may differ from
%       those expected---and not necessarily by a small amount.
%       For example, a ray may miss a plate it was expected to
%       hit and instead hit another plate considerably farther
%       from the ray's vertex, or miss the target entirely.
%
%   10) In the event that an intercept point lies on multiple
%       plates (that is, the point is on an edge or vertex),
%       a plate will be selected. Due to round-off error, the
%       selection may vary across platforms.
%
%   11) If any of the input arguments, `handle', `dladsc', `vertex' or
%       `raydir', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   12) If any of the input arguments, `handle', `dladsc', `vertex' or
%       `raydir', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   See the description of the input argument `handle'.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  A. Woo, "Fast Ray-Box Intersection", Graphic Gems I,
%        395-396, Aug. 1990
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
%       Edited the -Examples section to comply with NAIF standard. Updated
%       code example to prompt for the input DSK file. Corrected example's
%       problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 30-MAR-2014 (EDW) (NJB)
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
   catch spiceerr
      rethrow(spiceerr)
   end


