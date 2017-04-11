%-Abstract
%
%   CSPICE_SRFXPT computes the surface intercept point of a specified ray
%   on a target body at a specified epoch, optionally corrected for light
%   time and stellar aberration, given an observer and a direction vector
%   defining a ray.
%
%   Deprecated: This routine has been superseded by the routine
%   cspice_sincpt. This routine is supported for purposes of
%   backward compatibility only.
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
%      method   a string providing parameters defining
%               the computation method to use.
%
%               [1,c1] = size(method); char = class(method)
%
%                  or
%
%               [1,1] = size(method); cell = class(method)
%
%               The only currently supported choice:
%
%                  "Ellipsoid"   The intercept computation uses
%                                a triaxial ellipsoid to model
%                                the surface of the target body.
%                                The ellipsoid's radii must be
%                                available in the kernel pool.
%
%               Neither case nor white space are significant in
%               'method'.
%
%      target   the name of the target body. The target body is an
%               ephemeris object (its trajectory is given by SPK data),
%               and is an extended object. Optionally, you may supply the
%               integer ID code for the object as an integer string, i.e.
%               both 'MOON' and '301' are legitimate strings that indicate
%               the Moon is the target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%      et       the ephemeris time(s), expressed as seconds past J2000
%               TDB, at which to compute the surface intercept point on
%               the target body (this epoch represents either the time of
%               signal reception, or transmission, depending on the
%               selected 'abcorr')
%
%               [1,n] = size(et); double = class(et)
%
%      abcorr   the aberration correction to apply when computing the
%               observer-target state and the target body orientation.
%
%               [1,c3] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               For practical purposes, 'CN' (converged Newtonian)
%               represents the best correction choice.
%
%               'abcorr' may be any of the following:
%
%                  'NONE'     Apply no correction. Return the
%                             geometric surface intercept point on the
%                             target body.
%
%               Let 'lt' represent the one-way light time between the
%               observer and the surface intercept point (note: NOT
%               between the observer and the target body's center).
%               The following values of 'abcorr' apply to the
%               "reception" case in which photons depart from the
%               intercept point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at 'et':
%
%
%                  'LT'       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             yields the location of the surface
%                             intercept point at the moment it
%                             emitted photons arriving at the
%                             observer at 'et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             'LT' option uses one iteration.
%
%                             Both the target state as seen by the
%                             observer, and rotation of the target
%                             body, are corrected for light time.
%
%                  'LT+S'     Correct for one-way light time and
%                             stellar aberration using a Newtonian
%                             formulation. This option modifies the
%                             state obtained with the 'LT' option to
%                             account for the observer's velocity
%                             relative to the solar system
%                             barycenter. The result is the apparent
%                             surface intercept point as seen by the
%                             observer.
%
%                  'CN'       Converged Newtonian light time
%                             correction.  In solving the light time
%                             equation, the 'CN' correction iterates
%                             until the solution converges. Both the
%                             state and rotation of the target body
%                             are corrected for light time.
%
%                  'CN+S'     Converged Newtonian light time
%                             and stellar aberration corrections.
%
%               The following values of 'abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at 'et' and arrive at the
%               intercept point at the light-time corrected epoch
%               et+lt:
%
%                  'XLT'      "Transmission" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation. This correction yields the
%                             intercept location at the moment it
%                             receives photons emitted from the
%                             observer's location at 'et'.
%
%                             The light time correction uses an
%                             iterative solution of the light time
%                             equation. The solution invoked by the
%                             'LT' option uses one iteration.
%
%                             Both the target state as seen by the
%                             observer, and rotation of the target
%                             body, are corrected for light time.
%
%                  'XLT+S'    "Transmission" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation  This option modifies the
%                             intercept obtained with the 'XLT'
%                             option to account for the observer's
%                             velocity relative to the solar system
%                             barycenter.
%
%                  'XCN'      Converged Newtonian light time
%                             correction.  This is the same as 'XLT'
%                             correction but with further iterations
%                             to a converged Newtonian light time
%                             solution.
%
%                  'XCN+S'    "Transmission" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%      obsrvr   the name of a observing body. Optionally, you may supply
%               the integer ID code for the object as an integer string,
%               i.e. both 'MOON' and '301' are legitimate strings that
%               indicate the Moon is the observing body.
%
%               [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%      dref     the name of the reference frame relative to which the
%               input direction vector is expressed. This may be any
%               frame supported by the SPICE system, including built-in
%               frames (documented in the Frames Required Reading) and
%               frames defined by a loaded frame kernel (FK).
%
%               [1,c5] = size(dref); char = class(dref)
%
%                  or
%
%               [1,1] = size(dref); cell = class(dref)
%
%               When 'dref' designates a non-inertial frame, the
%               orientation of the frame is evaluated at an epoch
%               dependent on the frame's center and, if the center is
%               not the observer, on the selected aberration
%               correction. See the description of the direction
%               vector 'dvec' for details.
%
%      dvec     Pointing vector emanating from the observer.  The
%               intercept with the target body's surface of the ray
%               defined by the observer and 'dvec' is sought.
%
%               [3,1] = size(dvec); double = class(dvec)
%
%               'dvec' is specified relative to the reference frame
%               designated by 'dref'.
%
%               Non-inertial reference frames are treated as follows:
%               if the center of the frame is at the observer's
%               location, the frame is evaluated at 'et'.  If the
%               frame's center is located elsewhere, then letting
%               'ltcent' be the one-way light time between the observer
%               and the central body associated with the frame, the
%               orientation of the frame is evaluated at et-ltcent,
%               et+ltcent, or 'et' depending on whether the requested
%               aberration correction is, respectively, for received
%               radiation, transmitted radiation, or is omitted.
%               'ltcent' is computed using the method indicated by
%               'abcorr'.
%
%   the call:
%
%      [ spoint, dist, trgepc, obspos, found ] = ...
%            cspice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)
%
%   returns:
%
%      spoint   the surface intercept point on the target body of
%               the ray defined by the observer and the direction
%               vector. If the ray intersects the target body in
%               multiple points, the selected intersection point is
%               the one closest to the observer.  The output
%               argument 'found' (see below) indicates whether an
%               intercept was found.
%
%               [3,1] = size(spoint); double = class(spoint)
%
%               'spoint' is expressed in Cartesian coordinates,
%               relative to the body-fixed frame associated with the
%               target body.  The body-fixed target frame is
%               evaluated at the intercept epoch 'trgepc' (see
%               description below).
%
%               When light time correction is used, the duration of
%               light travel between 'spoint' to the observer is
%               considered to be the one way light time.  When both
%               light time and stellar aberration corrections are
%               used, 'spoint' is selected such that, when 'spoint'
%               is corrected for light time and the vector from the
%               observer to the light-time corrected location of
%               'spoint' is corrected for stellar aberration, the
%               resulting vector is parallel to the ray defined by
%               the observer's location and 'dvec'.
%
%               The components of 'spoint' are given in units of km.
%
%     dist      the distance between the observer and the surface
%               intercept on the target body. 'dist' is given in
%               units of km.
%
%               [1,1] = size(dist); double = class(dist)
%
%     trgepc    the "intercept epoch."  This is the epoch at which
%               the ray defined by 'obsrvr' and 'dvec' intercepts the
%               target surface at 'spoint'.  'trgepc' is defined as
%               follows: letting 'lt' be the one-way light time
%               between the observer and the intercept point,
%               'trgepc' is the epoch et-lt, et+lt, or 'et' depending
%               on whether the requested aberration correction is,
%               respectively, for received radiation, transmitted
%               radiation, or omitted. 'lt' is computed using the
%               method indicated by 'abcorr'.
%
%               [1,1] = size(trgepc); double = class(trgepc)
%
%               'trgepc' is expressed as seconds past J2000 TDB.
%
%     obspos    the vector from the center of the target body at
%               epoch 'trgepc' to the observer at epoch 'et'.
%               'obspos' is expressed in the target body-fixed
%               reference frame evaluated at 'trgepc'.  (This is
%               the frame relative to which 'spoint' is given.)
%
%               [3,1] = size(obspos); double = class(obspos)
%
%               'obspos' is returned to simplify various related
%               computations that would otherwise be cumbersome. For
%               example, the vector 'xvec' from the observer to
%               'spoint' can be calculated via
%
%                  xvec = spoint - obspos
%
%               The components of 'obspos' are given in units of km.
%
%     found     the logical flag indicating whether or not the ray
%               intersects the target.  If an intersection exists
%               'found' will be returned as true. If the ray misses
%               the target, 'found' will return as false.
%
%               [1,1] = size(found); logical = class(found)
%
%      'spoint', 'dist', 'trgepc', 'obspos(3)', and 'found' return with the
%      same vectorization measure (N) as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Assign needed constants.
%      %
%      BUSID   = -94000;
%      MGS     = -94;
%      NCORNR  = 4;
%      ABCORR  = 'LT+S';
%      CAMERA  = 'MGS_MOC_NA';
%      DREF    = 'J2000';
%      METHOD  = 'ELLIPSOID';
%      OBSRVR  = 'MGS';
%      TARGET  = 'MARS';
%      UTC     = '2003 OCT 13 06:00:00 UTC';
%
%      %
%      %    Load kernel files:
%      %
%      %       - Leapseconds kernel
%      %       - MGS SCLK kernel
%      %       - Text PCK file
%      %       - Planetary SPK file
%      %       - MGS I-kernel
%      %       - MGS spacecraft bus C-kernel
%      %       - MGS SPK file
%      %
%      cspice_furnsh( '/kernels/gen/lsk/naif0008.tls' )
%      cspice_furnsh( '/kernels/MGS/sclk/MGS_SCLKSCET.00053.tsc' )
%      cspice_furnsh( '/kernels/MGS/pck/mars_iau2000_v0.tpc' )
%      cspice_furnsh( '/kernels/gen/spk/de405s.bsp' )
%      cspice_furnsh( '/kernels/MGS/ik/moc20.ti' )
%      cspice_furnsh( '/kernels/MGS/ck/mgs_sc_ext12.bc' )
%      cspice_furnsh( '/kernels/MGS/spk/mgs_ext12.bsp' )
%
%      %
%      % Convert the UTC request time to ET (seconds past
%      % J2000, TDB).
%      %
%      et = cspice_str2et( UTC );
%
%      %
%      % Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
%      % ID code.  Then look up the field of view (FOV)
%      % parameters.
%      %
%      [camid, found] = cspice_bodn2c( CAMERA );
%
%      [shape, dref, bsight, bounds] = cspice_getfov( camid, NCORNR);
%
%      disp( ' ' )
%      disp( 'Surface Intercept Locations for Camera' )
%      disp( 'FOV Boundary and Boresight Vectors'     )
%      disp( ' ' )
%
%      txt = sprintf( '   Instrument:             %s', CAMERA);
%      disp( txt )
%
%      txt = sprintf( '   Epoch:                  %s', UTC);
%      disp( txt )
%
%      txt = sprintf( '   Aberration correction:  %s', ABCORR);
%      disp( txt )
%      disp( ' ' )
%
%
%      %
%      % Now compute and display the surface intercepts for the
%      % boresight and all of the FOV boundary vectors.
%      %
%      for i=1:NCORNR+1
%
%         if( i <= NCORNR )
%
%            %
%            % 'bounds' represents a 3 X NCORNR array with each row a bounds
%            % vector. Extract the vectors from 'bounds' using as a vector
%            % segment.
%            %
%            %    corner vector 0: bounds(:,1)
%            %    corner vector 1: bounds(:,2)
%            %    corner vector 2: bounds(:,3)
%            %    corner vector 3: bounds(:,4)
%            %
%            %
%            title = sprintf( 'Corner vector %d', i );
%            dvec = bounds(:,i);
%
%         else
%
%            title = sprintf( 'Boresight vector' );
%            dvec = bsight;
%
%         end
%
%         %
%         % Compute the surface intercept point using
%         % the specified aberration corrections.
%         %
%         [spoint, dist, trgepc, obspos, found] = ...
%               cspice_srfxpt( METHOD, TARGET, et, ABCORR, OBSRVR, dref, dvec );
%
%
%         if( found )
%
%            %
%            % Convert rectangular coordinates to planetocentric
%            % latitude and longitude.  Convert radians to degrees.
%            %
%            [ radius, lon, lat ] = cspice_reclat( spoint );
%
%            lon = lon * cspice_dpr;
%            lat = lat * cspice_dpr;
%
%            %
%            % Display the results.
%            %
%            disp( title )
%            disp( ' ' )
%
%            txt = sprintf( '  Vector in %s frame = ', dref );
%            disp( txt )
%
%            txt = sprintf( '%18.10e%18.10e%18.10e', dvec );
%            disp( txt )
%            disp( ' ' )
%
%            disp( '  Intercept:' )
%            disp( ' ' )
%
%            txt = ...
%              sprintf('     Radius                   (km)  = %18.10e', radius);
%            disp( txt )
%
%            txt = ...
%              sprintf('     Planetocentric Latitude  (deg) = %18.10e', lat);
%            disp( txt )
%
%            txt = ...
%              sprintf('     Planetocentric Longitude (deg) = %18.10e', lon);
%            disp( txt )
%
%            txt = ...
%              sprintf('     Range                    (km)  = %18.10e', dist);
%            disp( txt )
%            disp( ' ' )
%
%         else
%
%            disp( 'Intercept not found.' )
%
%         end
%
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in MATLAB due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Surface Intercept Locations for Camera
%      FOV Boundary and Boresight Vectors
%
%         Instrument:             MGS_MOC_NA
%         Epoch:                  2003 OCT 13 06:00:00 UTC
%         Aberration correction:  LT+S
%
%      Corner vector 1
%
%        Vector in MGS_MOC_NA frame =
%        1.8571383810e-06 -3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849412615e+03
%           Planetocentric Latitude  (deg) =  -4.8477118861e+01
%           Planetocentric Longitude (deg) =  -1.2347365507e+02
%           Range                    (km)  =   3.8898362744e+02
%
%      Corner vector 2
%
%        Vector in MGS_MOC_NA frame =
%        1.8571383810e-06  3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849398244e+03
%           Planetocentric Latitude  (deg) =  -4.8481272936e+01
%           Planetocentric Longitude (deg) =  -1.2339839939e+02
%           Range                    (km)  =   3.8897565851e+02
%
%      Corner vector 3
%
%        Vector in MGS_MOC_NA frame =
%       -1.8571383810e-06  3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849398156e+03
%           Planetocentric Latitude  (deg) =  -4.8481298506e+01
%           Planetocentric Longitude (deg) =  -1.2339840260e+02
%           Range                    (km)  =   3.8897519958e+02
%
%      Corner vector 4
%
%        Vector in MGS_MOC_NA frame =
%       -1.8571383810e-06 -3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849412527e+03
%           Planetocentric Latitude  (deg) =  -4.8477144435e+01
%           Planetocentric Longitude (deg) =  -1.2347365823e+02
%           Range                    (km)  =   3.8898316850e+02
%
%      Boresight vector
%
%        Vector in MGS_MOC_NA frame =
%        0.0000000000e+00  0.0000000000e+00  1.0000000000e+00
%
%        Intercept:
%
%           Radius                   (km)  =   3.3849405358e+03
%           Planetocentric Latitude  (deg) =  -4.8479216591e+01
%           Planetocentric Longitude (deg) =  -1.2343603019e+02
%           Range                    (km)  =   3.8897626607e+02
%
%-Particulars
%
%   A sister version of this routine exists named mice_srfxpt that returns
%   the output arguments as fields in a single structure.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine srfxpt_c.
%
%   MICE.REQ
%   FRAMES.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.3, 12-FEB-2015, EDW (JPL)
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.2, 18-MAY-2010, BVS (JPL)
%
%      Index line now states that this routine is deprecated.
%
%   -Mice Version 1.0.1, 11-NOV-2008, EDW (JPL)
%
%      Edits to header; Abstract now states that this routine is
%      deprecated.
%
%   -Mice Version 1.0.0, 03-JAN-2006, EDW (JPL)
%
%-Index_Entries
%
%   DEPRECATED surface intercept point
%
%-&

function [spoint, dist, trgepc, obspos, found] = ...
          cspice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)

   switch nargin
      case 7

         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         dref   = zzmice_str(dref);
         dvec   = zzmice_dp(dvec);

      otherwise

         error( [ 'Usage: [_spoint(3)_, _dist_, _trgepc_, '   ...
                  '_obspos(3)_, _found_ ] = '                 ...
                  'cspice_srfxpt( `method`, `target`, _et_, ' ...
                  '`abcorr`, `obsrvr`, `dref`, dvec(3))']  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [surf] = mice('srfxpt_s', method, target, et, abcorr, obsrvr, dref, dvec);
      spoint = reshape( [surf.spoint], 3, [] );
      dist   = reshape( [surf.dist]  , 1, [] );
      trgepc = reshape( [surf.trgepc], 1, [] );
      obspos = reshape( [surf.obspos], 3, [] );
      found  = reshape( [surf.found] , 1, [] );
   catch
      rethrow(lasterror)
   end

