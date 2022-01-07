%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   cspice_latsrf. This routine is supported for purposes of backward
%   compatibility only.
%
%   CSPICE_LLGRID_PL02, given the planetocentric longitude and latitude
%   values of a set of surface points on a specified target body, compute
%   the corresponding rectangular coordinates of those points. The
%   target body's surface is represented by a triangular plate model
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
%               This kernel must contain a type 2 segment that provides a
%               plate model representing the entire surface of the target
%               body.
%
%      dladsc   the DLA descriptor of a DSK segment representing
%               the surface of a target body.
%
%               [SPICE_DLA_DSCSIZ,1] = size(dladsc); int32 = class(dladsc)
%
%      grid     an array of planetocentric longitude/latitude pairs
%               to be mapped to surface points on the target body.
%
%               [2,n] = size(grid); double = class(grid)
%
%               Elements
%
%                  grid(1,i)
%                  grid(2,i)
%
%               are, respectively, the planetocentric longitude and
%               latitude of the ith grid point.
%
%               Units are radians.
%
%   the call:
%
%      [srfpts, pltids] = cspice_llgrid_pl02( handle, dladsc, grid )
%
%   returns:
%
%      srfpts   an array containing the rectangular (Cartesian)
%               coordinates of the surface points on the target body,
%               expressed relative to the body-fixed reference frame of
%               the target body, corresponding to the input grid points.
%
%               [3,n] = size(srfpts); double = class(srfpts)
%
%      pltids   an array of integer ID codes of the plates on which
%               the surface points are located.
%
%               [1,n] = size(pltids); int32 = class(pltids)
%
%               The ith plate ID corresponds to the ith surface point. These
%               ID codes can be use to look up data associated with the
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
%   1) Find the surface points on a target body corresponding to a given
%      planetocentric longitude/latitude grid.
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
%      function llgrid_pl02_ex1( dsknam )
%
%         %
%         % Constants
%         %
%         NLAT     =  9;
%         NLON     =  9;
%         MAXGRID  =  NLAT * NLON;
%         TOL      =  1.d-12;
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
%         % If we made it this far, DLADSC is the
%         % DLA descriptor of the first segment.
%         %
%         % Now generate the grid points.  We generate
%         % points along latitude bands, working from
%         % north to south.  The latitude range is selected
%         % to range from +80 to -80 degrees.  Longitude
%         % ranges from 0 to 320 degrees.  The increment
%         % is 20 degrees for latitude and 40 degrees for
%         % longitude.
%         %
%
%         grid = zeros(2,MAXGRID);
%         n    = 1;
%
%         for  i = 0:(NLAT-1)
%
%            lat = cspice_rpd() * ( 80.0 - 20.0*i );
%
%            for  j = 0:(NLON-1)
%
%               lon = cspice_rpd() * 40.0*j;
%
%               grid(1,n) = lon;
%               grid(2,n) = lat;
%
%               n = n + 1;
%
%            end
%
%         end
%
%         npoints = n - 1;
%
%         %
%         % Find the surface points corresponding to the grid points.
%         %
%         [srfpts, pltids] = cspice_llgrid_pl02( handle, dladsc, grid );
%
%         %
%         % fprintf out the surface points in latitudinal
%         % coordinates and compare the derived lon/lat values
%         % to those of the input grid.
%         %
%         for  i = 1:npoints
%
%            %
%            % Use recrad_c rather than reclat_c to produce
%            % non-negative longitudes.
%            %
%            [ xr, xlon, xlat] = cspice_recrad( srfpts(:,i) );
%
%            fprintf( '\n\nIntercept for grid point  %d\n', i )
%            fprintf( '   Plate ID:              %d\n', pltids(i) )
%            fprintf( '   Cartesian Coordinates: %f  %f  %f\n', srfpts(:,i))
%            fprintf( '   Latitudinal Coordinates:\n')
%            fprintf( '   Longitude (deg): %f\n', xlon * cspice_dpr() )
%            fprintf( '   Latitude  (deg): %f\n', xlat * cspice_dpr() )
%            fprintf( '   Radius     (km): %f\n', xr )
%
%            fprintf( '\nOriginal grid coordinates:\n' )
%            fprintf( '   Longitude (deg): %f\n', grid(1,i) * cspice_dpr() )
%            fprintf( '   Latitude  (deg): %f\n', grid(2,i) * cspice_dpr() )
%
%            %
%            % Perform sanity checks on the intercept
%            % coordinates.  Stop the program if any error
%            % is larger than our tolerance value.
%            %
%            lon = grid(1,i);
%            lat = grid(2,i);
%
%            if ( abs(xlat-lat) > TOL )
%
%               fprintf( 'Latitude error!' )
%               return
%
%            end
%
%            if ( abs(xlon - lon) > cspice_pi() )
%
%               if ( xlon > lon )
%                  xlon = xlon - cspice_twopi()
%               else
%                  xlon = xlon + cspice_twopi()
%               end
%
%            end
%
%            if  ( abs(xlon - lon)  > TOL )
%
%               fprintf( 'Longitude error!\n' )
%               return
%
%            end
%
%         end
%
%         %
%         % Close the kernel.
%         %
%         cspice_dascls( handle )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, with the following variables as inputs
%
%         dsknam = 'phobos_3_3.bds';
%
%      the output was:
%
%
%      Intercept for grid point  1
%         Plate ID:              306238
%         Cartesian Coordinates: 1.520878  0.000000  8.625327
%         Latitudinal Coordinates:
%         Longitude (deg): 0.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 8.758387
%
%      Original grid coordinates:
%         Longitude (deg): 0.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  2
%         Plate ID:              317112
%         Cartesian Coordinates: 1.189704  0.998280  8.807772
%         Latitudinal Coordinates:
%         Longitude (deg): 40.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 8.943646
%
%      Original grid coordinates:
%         Longitude (deg): 40.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  3
%         Plate ID:              324141
%         Cartesian Coordinates: 0.277775  1.575341  9.072029
%         Latitudinal Coordinates:
%         Longitude (deg): 80.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.211980
%
%      Original grid coordinates:
%         Longitude (deg): 80.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  4
%         Plate ID:              327994
%         Cartesian Coordinates: -0.810824  1.404388  9.196823
%         Latitudinal Coordinates:
%         Longitude (deg): 120.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.338699
%
%      Original grid coordinates:
%         Longitude (deg): 120.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  5
%         Plate ID:              329431
%         Cartesian Coordinates: -1.478202  0.538022  8.921321
%         Latitudinal Coordinates:
%         Longitude (deg): 160.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.058947
%
%      Original grid coordinates:
%         Longitude (deg): 160.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  6
%         Plate ID:              196042
%         Cartesian Coordinates: -1.498548  -0.545427  9.044113
%         Latitudinal Coordinates:
%         Longitude (deg): 200.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.183633
%
%      Original grid coordinates:
%         Longitude (deg): 200.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  7
%         Plate ID:              235899
%         Cartesian Coordinates: -0.782405  -1.355164  8.874473
%         Latitudinal Coordinates:
%         Longitude (deg): 240.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 9.011376
%
%      Original grid coordinates:
%         Longitude (deg): 240.000000
%         Latitude  (deg): 80.000000
%
%
%      Intercept for grid point  8
%         Plate ID:              266998
%         Cartesian Coordinates: 0.264512  -1.500123  8.638862
%         Latitudinal Coordinates:
%         Longitude (deg): 280.000000
%         Latitude  (deg): 80.000000
%         Radius     (km): 8.772130
%
%      [...]
%
%
%      Warning: incomplete output. Only 100 out of 1053 lines have
%      been provided.
%
%
%-Particulars
%
%   See the headers of the Mice routines
%
%      cspice_reclat
%      cspice_latrec
%
%   for detailed definitions of Planetocentric coordinates.
%
%-Exceptions
%
%   If any of the listed errors occur, the output arguments are
%   left unchanged.
%
%   1)  If a DSK providing a DSK type 2 plate model has not been
%       loaded prior to calling llgrid_pl02, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the segment associated with the input DLA descriptor is not
%       of data type 2, the error SPICE(WRONGDATATYPE) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If a surface point cannot be computed because the ray corresponding
%       to a longitude/latitude pair fails to intersect the target
%       surface as defined by the plate model, the error
%       SPICE(NOINTERCEPT) is signaled by a routine in the call tree of this
%       routine.
%
%   4)  If any of the input arguments, `handle', `dladsc', or `grid',
%       is undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `handle', `dladsc', or `grid',
%       is not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice interface.
%
%-Files
%
%   The following data are required:
%
%   -  DSK data:  a DSK file containing a plate model representing the
%      target body's surface must be loaded. This kernel must contain
%      a type 2 segment that contains data for the entire surface of
%      the target body.
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
%   2)  The results returned by this routine may not be meaningful
%       if the target surface has multiple surface points associated
%       with some (longitude, latitude) coordinates.
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
%       argument names "spoints" and "plateIDs" to "srfpts" and "pltids"
%       for consistency with other functions.
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
%   DEPRECATED map latitudinal grid to DSK type 2 plate model
%
%-&

function [srfpts, pltids] = cspice_llgrid_pl02( handle, dladsc, grid )

   switch nargin
      case 3

         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         grid   = zzmice_dp( grid );

      otherwise

         error ( [ 'Usage: [srfpts(3), pltids] = '                         ...
            'cspice_llgrid_pl02( handle, dladsc(SPICE_DLA_DSCSIZ), grid )' ] )

   end

   %
   % Call the MEX library.
   %
   try

      [srfpts, pltids] = mice( 'llgrid_pl02', handle, dladsc, grid );
   catch spiceerr
      rethrow(spiceerr)
   end


