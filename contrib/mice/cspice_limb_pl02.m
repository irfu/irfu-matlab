%-Abstract
%
%   Deprecated: This routine has been superseded by the CSPICE routine
%   cspice_limbpt. This routine is supported for purposes of backward
%   compatibility only.
%
%   CSPICE_LIMB_PL02 returns a set of points on the limb of a specified
%   target body, where the target body's surface is represented by a
%   triangular plate model contained in a type 2 DSK segment.
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
%                  surface of a target body.
%
%                  [SPICE_DLA_DSCSIZ,1]  = size(dladsc)
%                                  int32 = class(dladsc)
%
%      target      the name of the target body. `target' is
%                  case-insensitive, and leading and trailing blanks in
%                  `target' are not significant. Optionally, you may supply
%                  a string containing the integer ID code for the object.
%                  For example both 'MOON' and '301' are legitimate strings
%                  that indicate the moon is the target body.
%
%                  This routine assumes that a kernel variable representing
%                  the target's radii is present in the kernel pool.
%                  Normally the kernel variable would be defined by loading
%                  a PCK file.
%
%                  [1,c1] = size(target); char = class(target)
%
%                     or
%
%                  [1,1] = size(target); cell = class(target)
%
%      et          the epoch of participation of the observer,
%                  expressed as ephemeris seconds past J2000 TDB: `et' is
%                  the epoch at which the observer's position is
%                  computed.
%
%                  When aberration corrections are not used, `et' is also
%                  the epoch at which the position and orientation of the
%                  target body are computed.
%
%                  When aberration corrections are used, `et' is the epoch
%                  at which the observer's position relative to the solar
%                  system barycenter is computed; in this case the position
%                  and orientation of the target body are computed at
%                  et-lt, where `lt' is the one-way light time between the
%                  target body's center and the observer. See the
%                  description of `abcorr' below for details.
%
%                  [1,1] = size(et); double = class(et)
%
%      fixfrm      the name of the reference frame relative to which the
%                  output limb points are expressed. This must a
%                  body-centered, body-fixed frame associated with the
%                  target. The frame's axes must be compatible with the
%                  triaxial ellipsoidal shape model associated with the
%                  target body (normally provide via a PCK): this routine
%                  assumes that the first, second, and third ellipsoid radii
%                  correspond, respectively, to the x, y, and z-axes of the
%                  frame designated by `fixfrm'.
%
%                  [1,c2] = size(fixfrm); char = class(fixfrm)
%
%                     or
%
%                  [1,1] = size(fixfrm); cell = class(fixfrm)
%
%                  `fixfrm' may refer to a built-in frame (documented in
%                  the Frames Required Reading) or a frame defined by a
%                  loaded frame kernel (FK).
%
%                  The orientation of the frame designated by `fixfrm' is
%                  evaluated at epoch of participation of the target
%                  body. See the descriptions of `et' and `abcorr' for
%                  details.
%
%      abcorr      indicates the aberration correction to be applied
%                  when computing the observer-target position, the
%                  orientation of the target body, and the target-
%                  source position vector.
%
%                  [1,c3] = size(abcorr); char = class(abcorr)
%
%                     or
%
%                  [1,1] = size(abcorr); cell = class(abcorr)
%
%                  `abcorr' may be any of the following.
%
%                     'NONE'     Apply no correction. Compute the limb
%                                points using the position of the observer
%                                and target, and the orientation of the
%                                target, at `et'.
%
%                  Let `lt' represent the one-way light time between the
%                  observer and the target body's center. The following
%                  values of `abcorr' apply to the "reception" case in
%                  which photons depart from the target body's center at
%                  the light-time corrected epoch et-lt and *arrive* at
%                  the observer's location at `et':
%
%
%                     'LT'       Correct for one-way light time (also
%                                called "planetary aberration") using a
%                                Newtonian formulation. This correction
%                                yields the location of the limb points at
%                                the approximate time they emitted photons
%                                arriving at the observer at `et' (the
%                                difference between light time to the
%                                target center and light time to the limb
%                                points is ignored).
%
%                                The light time correction uses an
%                                iterative solution of the light time
%                                equation. The solution invoked by the
%                                'LT' option uses one iteration.
%
%                                The target position as seen by the
%                                observer and the rotation of the target
%                                body are corrected for light time.
%
%                     'LT+S'     Correct for one-way light time and stellar
%                                aberration using a Newtonian formulation.
%                                This option modifies the position obtained
%                                with the 'LT' option to account for the
%                                observer's velocity relative to the solar
%                                system barycenter. The result is the
%                                apparent limb as seen by the observer.
%
%                     'CN'       Converged Newtonian light time correction.
%                                In solving the light time equation, the
%                                'CN' correction iterates until the
%                                solution converges. The position and
%                                rotation of the target body are corrected
%                                for light time.
%
%                     'CN+S'     Converged Newtonian light time
%                                and stellar aberration corrections.
%
%
%      obsrvr      the name of the observing body. This is typically
%                  a spacecraft, the Earth, or a surface point on the
%                  Earth. `obsrvr' is case-insensitive, and leading and
%                  trailing blanks in `obsrvr' are not significant.
%                  Optionally, you may supply a string containing the
%                  integer ID code for the object. For example both
%                  'EARTH' and '399' are legitimate strings that indicate
%                  the Earth is the observer.
%
%                  [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                     or
%
%                  [1,1] = size(obsrvr); cell = class(obsrvr)
%
%      npoints     the number of limb points to compute.
%
%                  [1,1] = size(npoints); int32 = class(npoints)
%
%                   For values of `npoints' less-than or equal-to zero,
%                   the output arguments return as zeros and empty arrays.
%
%   the call:
%
%      [trgepc, obspos, limbpts, plateids] =                  ...
%                   cspice_limb_pl02( handle, dladsc, target, ...
%                                     et,     fixref, abcorr, ...
%                                     obsrvr, npoints )
%
%   returns:
%
%      trgepc      the "target epoch."  `trgepc' is defined as follows:
%                  letting `lt' be the one-way light time between the
%                  target center and observer, `trgepc' is either the
%                  epoch et-lt or `et' depending on whether the requested
%                  aberration correction is, respectively, for received
%                  radiation or omitted. `lt' is computed using the
%                  method indicated by `abcorr'.
%
%                  [1,1] = size(trgepc); double = class(trgepc)
%
%                  `trgepc' is expressed as seconds past J2000 TDB.
%
%      obspos      the vector from the center of the target body at
%                  epoch `trgepc' to the observer at epoch `et'. `obspos' is
%                  expressed in the target body-fixed reference frame
%                  `fixfrm', which is evaluated at `trgepc'.
%
%                  [3,1] = size(obspos); double = class(obspos)
%
%                  `obspos' is returned to simplify various related
%                  computations that would otherwise be cumbersome. For
%                  example, the vector `xvec' from the observer to the
%                  Ith limb point can be calculated via the expression
%
%                     xvec = imbpts(:,i) - obspos
%
%                  The components of `obspos' are given in units of km.
%
%      limbpts     an array of points on the limb of the target.
%
%                  [3,npoints] = size(limbpts); double = class(limbpts)
%
%                  The ith point is contained in the array elements
%
%                      limbpts(:,i)
%
%                  As described above, each limb point lies on a ray
%                  emanating from the center of the target and passing
%                  through a limb point on the target's reference
%                  ellipsoid. Each limb point *on the reference ellipsoid*
%                  is the point of tangency of a ray that emanates from the
%                  observer. Measured in a cylindrical coordinate system
%                  whose Z-axis is parallel to the observer-target vector,
%                  the magnitude of the separation in longitude between the
%                  limb points is
%
%                     2*Pi / npoints
%
%                  The limb points are expressed in the body-fixed
%                  reference frame designated by `fixfrm'; the
%                  orientation of the frame is evaluated at `trgepc'.
%                  Units are km.
%
%      plateIDs    an array of integer ID codes of the plates on which
%                  the limb points are located. The ith plate ID
%                  corresponds to the ith limb point. These ID codes can
%                  be use to look up data associated with the plates, such
%                  as the plates' vertices or outward normal vectors.
%
%                  [1,npoints] = size(plateIDs); int32 = class(plateIDs)
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   In the following example program, the file
%
%      phobos_3_3.dsk
%
%   is a DSK file containing a type 2 segment that provides a plate model
%   representation of the surface of Phobos.  The file
%
%      mar097.bsp
%
%   is a binary SPK file containing data for Phobos, the Earth, and the
%   Sun for a time interval bracketing the date
%
%      2006 NOV 3 00:00:00 UTC.
%
%   pck00010.tpc is a planetary constants kernel file containing radii
%   and rotation model constants.  naif0010.tls is a leapseconds kernel.
%
%   All of the kernels other than the DSK file should be loaded via
%   a meta-kernel.  An example of the contents of such a kernel is:
%
%      KPL/MK
%
%      \begindata
%
%         KERNELS_TO_LOAD = ( 'naif0010.tls'
%                             'pck00010.tpc'
%                             'mar097.bsp' )
%      \begintext
%
%   Example:
%
%      %
%      % Compute a set of limb points on Phobos as seen from Mars. Perform
%      % a consistency check using the emission angle at each point,
%      % where the emission angle is computed using both a reference
%      % ellipsoid and the actual plate model surface and surface normal.
%      % We expect to see an emission angle of approximately 90 degrees.
%      %
%      %    Version 1.0.0 18-MAR-2014 (NJB)
%      %
%
%      function limb_pl02_t
%
%         %
%         % Constants
%         %
%         NPOINTS     = 3;
%         TIMLEN      = 40;
%         TOL         = 1.d-12;
%         UTCSTR      = '2007 FEB 9 00:00:00 UTC';
%
%         %
%         % Initial values
%         %
%         target      = 'Phobos';
%         abcorr      = 'CN+S';
%         fixfrm      = 'IAU_PHOBOS';
%         obsrvr      = 'Mars';
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
%         dsk = input( 'Enter DSK name         > ','s');
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
%         % Now compute sub-points using both computation
%         % methods. We'll vary the aberration corrections
%         % and the epochs.
%         %
%
%         et = cspice_str2et( UTCSTR );
%
%         timstr = cspice_timout( et,                         ...
%                        'YYYY-MON-DD HR:MN:SC.### ::TDB(TDB)' );
%
%
%         fprintf( '\n\n' )
%         fprintf( '   Observer:                %s\n', obsrvr )
%         fprintf( '   Target:                  %s\n', target )
%         fprintf( '   Observation epoch:       %s\n', timstr )
%         fprintf( '   Aberration correction:   %s\n', abcorr )
%         fprintf( '   Body-fixed frame:        %s\n', fixfrm )
%
%         %
%         % Now compute grid of limb points.
%         %
%         [trgepc, obspos, limbpts, plateIDs] =              ...
%                  cspice_limb_pl02( handle, dladsc, target, ...
%                                    et,     fixfrm, abcorr, ...
%                                    obsrvr, NPOINTS           );
%
%         %
%         % Display the limb points.
%         %
%         for  i = 1:NPOINTS
%
%            [radius, lon, lat] = cspice_reclat( limbpts(:,i) );
%
%            fprintf( '\n' )
%            fprintf( '      Limb point:  %d \n', i )
%            fprintf( '         Radius                     (km):  %f\n', ...
%                                                                  radius)
%            fprintf( '         Planetocentric longitude   (deg): %f\n', ...
%                                                     lon * cspice_dpr() )
%
%            fprintf( '         Planetocentric latitude    (deg): %f\n', ...
%                                                     lat * cspice_dpr() )
%
%            fprintf( '         Plate ID:                         %d\n', ...
%                                                            plateIDs(i) )
%
%            %
%            % Compute the illumination angles using an ellipsoidal
%            % representation of the target's surface. The role of
%            % this representation is to provide an outward surface
%            % normal.
%            %
%            [phase, solar, emissn] = cspice_illum( target, et, abcorr, ...
%                                                   obsrvr, limbpts(:,i) );
%
%            fprintf( '            emission angle derived using\n' )
%            fprintf( '               - an ellipsoidal ' )
%            fprintf( 'reference surface        (deg): %f\n', ...
%                                                  emissn * cspice_dpr() )
%
%
%            %
%            % Compute the illumination angles at the limb point
%            % using the actual plate model surface normal.
%            %
%            [phase, solar, emissn] = cspice_illum_pl02( handle, dladsc, ...
%                                           target, et,  abcorr, obsrvr, ...
%                                           limbpts(:,i)                   );
%
%            fprintf( '               - plate model''s ' )
%            fprintf( 'surface and normal vector (deg): %f\n', ...
%                                                emissn * cspice_dpr() )
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
%   Observer:                Mars
%   Target:                  Phobos
%   Observation epoch:       2007-FEB-09 00:01:05.184 (TDB)
%   Aberration correction:   CN+S
%   Body-fixed frame:        IAU_PHOBOS
%
%      Limb point:  1
%         Radius                     (km):  11.563501
%         Planetocentric longitude   (deg): 91.739066
%         Planetocentric latitude    (deg):  -0.000811
%         Plate ID:                         229468
%            emission angle derived using
%               - an ellipsoidal reference surface        (deg): 90.001006
%               - plate model's surface and normal vector (deg): 110.821665
%
%      Limb point:  2
%         Radius                     (km):  9.537023
%         Planetocentric longitude   (deg): -87.847223
%         Planetocentric latitude    (deg):  59.998792
%         Plate ID:                         235885
%            emission angle derived using
%               - an ellipsoidal reference surface        (deg): 89.999961
%               - plate model's surface and normal vector (deg): 97.681554
%
%      Limb point:  3
%         Radius                     (km):  9.046773
%         Planetocentric longitude   (deg): -88.051726
%         Planetocentric latitude    (deg): -59.997991
%         Plate ID:                         17961
%            emission angle derived using
%               - an ellipsoidal reference surface        (deg): 89.996966
%               - plate model's surface and normal vector (deg): 64.808794
%
%-Particulars
%
%   Boundaries of visible regions on an arbitrary surface are often
%   complicated point sets: boundaries of mountains and craters, if
%   present, may contribute to the overall set. To make the limb
%   computation tractable, we simplify the problem by using a reference
%   ellipsoid for guidance. We compute a set of limb points on the
%   reference ellipsoid for the target body, then use those points to
%   define the latitudes and longitudes of limb points on the surface
%   defined by the specified triangular shape model. As such, the set
%   of limb points found by this routine is just an approximation.
%
%-Required Reading
%
%   For important details concerning this module's function, please
%   refer to the CSPICE routine limb_pl02.
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
%   find limb on a dsk type_2 model
%   find limb on a triangular plate_model
%
%-&

function [trgepc, obspos, limbpts, plateIDs] =              ...
                  cspice_limb_pl02( handle, dladsc, target, ...
                                    et,     fixref, abcorr, ...
                                    obsrvr, npoints           )

   switch nargin
      case 8

         handle  = zzmice_int( handle );
         dladsc  = zzmice_int( dladsc );
         target  = zzmice_str( target );
         et      = zzmice_dp( et );
         fixref  = zzmice_str( fixref );
         abcorr  = zzmice_str( abcorr );
         obsrvr  = zzmice_str( obsrvr );
         npoints = zzmice_int( npoints );

      otherwise

         error ( [ 'Usage: [trgepc, obspos(3), limbpts(3,npoints), ' ...
            'plateIDs(npoints)] = cspice_limb_pl02( handle, '        ...
            'dladsc(SPICE_DLA_DSCSIZ), '                             ...
            '`target`, et, `fixref`, `abcorr`, `obsrvr`, npoints )' ] )

   end

   %
   % Call the MEX library.
   %
   try

      [trgepc, obspos, limbpts, plateIDs] = mice( 'limb_pl02', ...
                                             handle, dladsc, target, et, ...
                                             fixref, abcorr, obsrvr, npoints );
   catch
      rethrow(lasterror)
   end


