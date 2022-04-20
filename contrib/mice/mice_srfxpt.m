%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   cspice_sincpt. This routine is supported for purposes of
%   backward compatibility only.
%
%   MICE_SRFXPT computes the surface intercept point of a specified ray
%   on a target body at a specified epoch, optionally corrected for light
%   time and stellar aberration, given an observer and a direction vector
%   defining a ray.
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
%               `method'.
%
%      target   the name of the target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               The target body is an ephemeris object (its trajectory is
%               given by SPK data), and is an extended object. Optionally,
%               you may supply the integer ID code for the object as an
%               integer string, i.e. both 'MOON' and '301' are legitimate
%               strings that indicate the Moon is the target body.
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
%               `abcorr' may be any of the following:
%
%                  'NONE'     Apply no correction. Return the
%                             geometric surface intercept point on the
%                             target body.
%
%               Let `lt' represent the one-way light time between the
%               observer and the surface intercept point (note: NOT
%               between the observer and the target body's center).
%               The following values of `abcorr' apply to the
%               "reception" case in which photons depart from the
%               intercept point's location at the light-time
%               corrected epoch et-lt and *arrive* at the observer's
%               location at `et':
%
%
%                  'LT'       Correct for one-way light time (also
%                             called "planetary aberration") using a
%                             Newtonian formulation. This correction
%                             yields the location of the surface
%                             intercept point at the moment it
%                             emitted photons arriving at the
%                             observer at `et'.
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
%                             correction. In solving the light time
%                             equation, the 'CN' correction iterates
%                             until the solution converges. Both the
%                             state and rotation of the target body
%                             are corrected for light time.
%
%                  'CN+S'     Converged Newtonian light time
%                             and stellar aberration corrections.
%
%               The following values of `abcorr' apply to the
%               "transmission" case in which photons *depart* from
%               the observer's location at `et' and arrive at the
%               intercept point at the light-time corrected epoch
%               et+lt:
%
%                  'XLT'      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation. This correction yields the
%                             intercept location at the moment it
%                             receives photons emitted from the
%                             observer's location at `et'.
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
%                  'XLT+S'    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation  This option modifies the
%                             intercept obtained with the 'XLT'
%                             option to account for the observer's
%                             velocity relative to the solar system
%                             barycenter.
%
%                  'XCN'      Converged Newtonian light time
%                             correction. This is the same as 'XLT'
%                             correction but with further iterations
%                             to a converged Newtonian light time
%                             solution.
%
%                  'XCN+S'    "Transmission" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%      obsrvr   the name of a observing body.
%
%               [1,c4] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Optionally, you may supply the integer ID code for the object
%               as an integer string, i.e. both 'MOON' and '301' are
%               legitimate strings that indicate the Moon is the observing
%               body.
%
%      dref     the name of the reference frame relative to which the
%               input direction vector is expressed.
%
%               [1,c5] = size(dref); char = class(dref)
%
%                  or
%
%               [1,1] = size(dref); cell = class(dref)
%
%               This may be any frame supported by the SPICE system, including
%               built-in frames (documented in the Frames Required Reading)
%               and frames defined by a loaded frame kernel (FK).
%
%               When `dref' designates a non-inertial frame, the
%               orientation of the frame is evaluated at an epoch
%               dependent on the frame's center and, if the center is
%               not the observer, on the selected aberration
%               correction. See the description of the direction
%               vector `dvec' for details.
%
%      dvec     the pointing vector emanating from the observer.
%
%               [3,1] = size(dvec); double = class(dvec)
%
%               The intercept with the target body's surface of the ray
%               defined by the observer and `dvec' is sought.
%
%               `dvec' is specified relative to the reference frame
%               designated by `dref'.
%
%               Non-inertial reference frames are treated as follows:
%               if the center of the frame is at the observer's
%               location, the frame is evaluated at `et'. If the
%               frame's center is located elsewhere, then letting
%               `ltcent' be the one-way light time between the observer
%               and the central body associated with the frame, the
%               orientation of the frame is evaluated at et-ltcent,
%               et+ltcent, or `et' depending on whether the requested
%               aberration correction is, respectively, for received
%               radiation, transmitted radiation, or is omitted.
%               `ltcent' is computed using the method indicated by
%               `abcorr'.
%
%   the call:
%
%      [surf] = mice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec )
%
%   returns:
%
%      surf     the structure(s) containing the results of the calculation.
%
%               [1,n] = size(surf); struct = class(surf)
%
%               Each structure consists of the fields:
%
%                  spoint   the surface intercept point on the target body of
%                           the ray defined by the observer and the direction
%                           vector.
%
%                           [3,1]  = size(surf(i).spoint);
%                           double = class(surf(i).spoint)
%
%                           If the ray intersects the target body in
%                           multiple points, the selected intersection point
%                           is the one closest to the observer. The output
%                           argument `found' (see below) indicates whether an
%                           intercept was found.
%
%                           `spoint' is expressed in Cartesian coordinates,
%                           relative to the body-fixed frame associated with
%                           the target body. The body-fixed target frame is
%                           evaluated at the intercept epoch `trgepc' (see
%                           description below).
%
%                           When light time correction is used, the duration
%                           of light travel between `spoint' to the observer
%                           is considered to be the one way light time. When
%                           both light time and stellar aberration corrections
%                           are used, `spoint' is selected such that, when
%                           `spoint' is corrected for light time and the
%                           vector from the observer to the light-time
%                           corrected location of `spoint' is corrected for
%                           stellar aberration, the resulting vector is
%                           parallel to the ray defined by the observer's
%                           location and `dvec'.
%
%                           The components of `spoint' are given in units of
%                           km.
%
%                  dist     the distance between the observer and the surface
%                           intercept on the target body.
%
%                           [1,1] = size(surf(i).dist);
%                           double = class(surf(i).dist)
%
%                           `dist' is given in units of km.
%
%                  trgepc   the "intercept epoch."
%
%                           [1,1]  = size(surf(i).trgepc);
%                           double = class(surf(i).trgepc)
%
%                           This is the epoch at which the ray defined by
%                           `obsrvr' and `dvec' intercepts the target
%                           surface at `spoint'. `trgepc' is defined as
%                           follows: letting `lt' be the one-way light time
%                           between the observer and the intercept point,
%                           `trgepc' is the epoch et-lt, et+lt, or `et'
%                           depending on whether the requested aberration
%                           correction is, respectively, for received
%                           radiation, transmitted radiation, or omitted.
%                           `lt' is computed using the method indicated by
%                           `abcorr'.
%
%                           `trgepc' is expressed as seconds past J2000 TDB.
%
%                  obspos   the vector from the center of the target body at
%                           epoch `trgepc' to the observer at epoch `et'.
%
%                           [3,1]  = size(surf(i).obspos);
%                           double = class(surf(i).obspos)
%
%                           `obspos' is expressed in the target body-fixed
%                           reference frame evaluated at `trgepc'. (This is
%                           the frame relative to which `spoint' is given.)
%
%                           `obspos' is returned to simplify various related
%                           computations that would otherwise be cumbersome.
%                           For example, the vector `xvec' from the observer
%                           to `spoint' can be calculated via
%
%                              xvec = spoint - obspos
%
%                           The components of `obspos' are given in units of
%                           km.
%
%                  found    the logical flag indicating whether or not the ray
%                           intersects the target.
%
%                           [1,1]   = size(surf(i).found);
%                           logical = class(surf(i).found)
%
%                           If an intersection exists `found' will be returned
%                           as true. If the ray misses the target, `found'
%                           will return as false.
%
%               `surf' return with the same vectorization measure, N, as `et'.
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
%   1) The following program computes surface intercept points on
%      Mars for the boresight and FOV boundary vectors of the MGS MOC
%      narrow angle camera. The intercepts are computed for a single
%      observation epoch. Light time and stellar aberration
%      corrections are used. For simplicity, camera distortion is
%      ignored.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: srfxpt_ex1.tm
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
%            de405s.bsp                       Planetary ephemeris
%            mars_iau2000_v0.tpc              Planet orientation and
%                                             radii
%            naif0011.tls                     Leapseconds
%            mgs_moc_v20.ti                   MGS MOC instrument
%                                             parameters
%            mgs_sclkscet_00061.tsc           MGS SCLK coefficients
%            mgs_sc_ext12.bc                  MGS s/c bus attitude
%            mgs_ext12_ipng_mgs95j.bsp        MGS ephemeris
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de405s.bsp',
%                                'mars_iau2000_v0.tpc',
%                                'naif0011.tls',
%                                'mgs_moc_v20.ti',
%                                'mgs_sclkscet_00061.tsc',
%                                'mgs_sc_ext12.bc',
%                                'mgs_ext12_ipng_mgs95j.bsp' )
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function srfxpt_ex1()
%
%         %
%         % Assign needed constants.
%         %
%         BUSID   = -94000;
%         MGS     = -94;
%         NCORNR  = 4;
%         ABCORR  = 'LT+S';
%         CAMERA  = 'MGS_MOC_NA';
%         DREF    = 'J2000';
%         METHOD  = 'ELLIPSOID';
%         OBSRVR  = 'MGS';
%         TARGET  = 'MARS';
%         UTC     = '2003 OCT 13 06:00:00 UTC';
%
%         %
%         %    Load kernel files:
%         %
%         %       - Leapseconds kernel
%         %       - MGS SCLK kernel
%         %       - Text PCK file
%         %       - Planetary SPK file
%         %       - MGS I-kernel
%         %       - MGS spacecraft bus C-kernel
%         %       - MGS SPK file
%         %
%         cspice_furnsh( 'srfxpt_ex1.tm' )
%
%         %
%         % Convert the UTC request time to ET (seconds past
%         % J2000, TDB).
%         %
%         et = cspice_str2et( UTC );
%
%         %
%         % Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
%         % ID code.  Then look up the field of view (FOV)
%         % parameters.
%         %
%         [camid, found] = cspice_bodn2c( CAMERA );
%
%         [shape, dref, bsight, bounds] = cspice_getfov( camid, NCORNR );
%
%         disp( ' ' )
%         disp( 'Surface Intercept Locations for Camera' )
%         disp( 'FOV Boundary and Boresight Vectors'     )
%         disp( ' ' )
%
%         txt = sprintf( '   Instrument:             %s', CAMERA);
%         disp( txt )
%
%         txt = sprintf( '   Epoch:                  %s', UTC);
%         disp( txt )
%
%         txt = sprintf( '   Aberration correction:  %s', ABCORR);
%         disp( txt )
%         disp( ' ' )
%
%         %
%         % Now compute and display the surface intercepts for the
%         % boresight and all of the FOV boundary vectors.
%         %
%         for i=1:NCORNR+1
%
%            if( i <= NCORNR )
%
%               %
%               % `bounds' represents a 3 X NCORNR array with each row a
%               % bounds vector. Extract the vectors from 'bounds' using
%               % as a vector segment.
%               %
%               %    corner vector 0: bounds(:,1)
%               %    corner vector 1: bounds(:,2)
%               %    corner vector 2: bounds(:,3)
%               %    corner vector 3: bounds(:,4)
%               %
%               %
%               title = sprintf( 'Corner vector %d', i );
%               dvec = bounds(:,i);
%
%            else
%
%               title = sprintf( 'Boresight vector' );
%               dvec = bsight;
%
%            end
%
%            %
%            % Compute the surface intercept point using
%            % the specified aberration corrections.
%            %
%            [surf] = mice_srfxpt( METHOD, TARGET, et, ABCORR, OBSRVR,     ...
%                                  dref, dvec );
%
%            if( surf.found )
%
%               %
%               % Convert rectangular coordinates to planetocentric
%               % latitude and longitude.  Convert radians to degrees.
%               %
%               [radius, lon, lat] = cspice_reclat( surf.spoint );
%
%               lon = lon * cspice_dpr;
%               lat = lat * cspice_dpr;
%
%               %
%               % Display the results.
%               %
%               disp( title )
%               disp( ' ' )
%
%               txt = sprintf( '  Vector in %s frame = ', dref );
%               disp( txt )
%
%               txt = sprintf( '   %18.10e%18.10e%18.10e', dvec );
%               disp( txt )
%               disp( ' ' )
%
%               disp( '  Intercept:' )
%               disp( ' ' )
%
%               txt = sprintf(['     Radius                    (km)'       ...
%                              ' = %18.10e'], radius);
%               disp( txt )
%
%               txt = sprintf(['     Planetocentric Latitude  (deg)'       ...
%                              ' = %18.10e'], lat);
%               disp( txt )
%
%               txt = sprintf(['     Planetocentric Longitude (deg)'       ...
%                              ' = %18.10e'], lon);
%               disp( txt )
%
%               txt = sprintf(['     Range                     (km)'       ...
%                              ' = %18.10e'], surf.dist );
%               disp( txt )
%               disp( ' ' )
%
%            else
%
%               disp( 'Intercept not found.' )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
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
%           1.8571383810e-06 -3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                    (km) =   3.3849411359e+03
%           Planetocentric Latitude  (deg) =  -4.8477481852e+01
%           Planetocentric Longitude (deg) =  -1.2347407883e+02
%           Range                     (km) =   3.8898310725e+02
%
%      Corner vector 2
%
%        Vector in MGS_MOC_NA frame =
%           1.8571383810e-06  3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                    (km) =   3.3849396988e+03
%           Planetocentric Latitude  (deg) =  -4.8481636267e+01
%           Planetocentric Longitude (deg) =  -1.2339882275e+02
%           Range                     (km) =   3.8897512490e+02
%
%      Corner vector 3
%
%        Vector in MGS_MOC_NA frame =
%          -1.8571383810e-06  3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                    (km) =   3.3849396899e+03
%           Planetocentric Latitude  (deg) =  -4.8481661837e+01
%           Planetocentric Longitude (deg) =  -1.2339882596e+02
%           Range                     (km) =   3.8897466598e+02
%
%      Corner vector 4
%
%        Vector in MGS_MOC_NA frame =
%          -1.8571383810e-06 -3.8015622659e-03  9.9999277403e-01
%
%        Intercept:
%
%           Radius                    (km) =   3.3849411271e+03
%           Planetocentric Latitude  (deg) =  -4.8477507428e+01
%           Planetocentric Longitude (deg) =  -1.2347408199e+02
%           Range                     (km) =   3.8898264817e+02
%
%      Boresight vector
%
%        Vector in MGS_MOC_NA frame =
%           0.0000000000e+00  0.0000000000e+00  1.0000000000e+00
%
%        Intercept:
%
%           Radius                    (km) =   3.3849404102e+03
%           Planetocentric Latitude  (deg) =  -4.8479579751e+01
%           Planetocentric Longitude (deg) =  -1.2343645375e+02
%           Range                     (km) =   3.8897573918e+02
%
%
%-Particulars
%
%   A sister version of this routine exists named cspice_srfxpt that returns
%   the structure field data as separate arguments.
%
%   Given a ray defined by a direction vector and the location of an
%   observer, mice_srfxpt computes the surface intercept point of the ray
%   on a specified target body. mice_srfxpt also determines the distance
%   between the observer and the surface intercept point.
%
%   When aberration corrections are used, this routine finds the
%   value of `spoint' such that, if `spoint' is regarded as an ephemeris
%   object, after the selected aberration corrections are applied to
%   the vector from the observer to `spoint', the resulting vector is
%   parallel to the direction vector `dvec'.
%
%   This routine computes light time corrections using light time
%   between the observer and the surface intercept point, as opposed
%   to the center of the target. Similarly, stellar aberration
%   corrections done by this routine are based on the direction of
%   the vector from the observer to the light-time corrected
%   intercept point, not to the target center. This technique avoids
%   errors due to the differential between aberration corrections
%   across the target body. Therefore it's valid to use aberration
%   corrections with this routine even when the observer is very
%   close to the intercept point, in particular when the
%   observer-intercept point distance is much less than the
%   observer-target center distance. It's also valid to use stellar
%   aberration corrections even when the intercept point is near or
%   on the limb (as may occur in occultation computations using a
%   point target).
%
%   When comparing surface intercept point computations with results
%   from sources other than SPICE, it's essential to make sure the
%   same geometric definitions are used.
%
%-Exceptions
%
%   If any of the listed errors occur, the output arguments are
%   left unchanged.
%
%   1)  If the input argument `method' is not recognized, an error
%       is signaled by a routine in the call tree of this
%       routine.
%
%   2)  If `target' cannot be mapped to an ID code, the error
%       SPICE(IDCODENOTFOUND) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If `obsrvr' cannot be mapped to an ID code, an error is signaled
%       by a routine in the call tree of this routine.
%
%   4)  If the input argument `abcorr' is invalid, an error
%       is signaled by a routine in the call tree of this
%       routine.
%
%   5)  If a body-fixed reference frame associated with the target
%       cannot be found, the error SPICE(NOFRAME) is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If `obsrvr' and `target' map to the same NAIF integer ID codes, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   7)  If frame definition data enabling the evaluation of the state
%       of the target relative to the observer in target body-fixed
%       coordinates have not been loaded prior to calling cspice_srfxpt, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   8)  If the specified aberration correction is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   9)  If insufficient ephemeris data have been loaded prior to
%       calling cspice_srfxpt, an error is signaled by a
%       routine in the call tree of this routine. Note that when
%       light time correction is used, sufficient ephemeris data
%       must be available to propagate the states of both observer
%       and target to the solar system barycenter.
%
%   10) If the computation method has been specified as "Ellipsoid"
%       and triaxial radii of the target body have not been loaded
%       into the kernel pool prior to calling cspice_srfxpt, an error is
%       signaled by a routine in the call tree of this routine.
%
%   11) If any of the radii of the target body are non-positive, an
%       error is signaled by a routine in the call tree of this
%       routine. The target must be an extended body.
%
%   12) If PCK data needed to define the target body-fixed frame have
%       not been loaded prior to calling cspice_srfxpt, an error is signaled
%       by a routine in the call tree of this routine.
%
%   13) If the reference frame designated by `dref' is not recognized
%       by the SPICE frame subsystem, an error is signaled
%       by a routine in the call tree of this routine.
%
%   14) If the direction vector `dvec' is the zero vector, an error
%       is signaled  by a routine in the call tree of this routine.
%
%   15) If any of the input arguments, `method', `target', `et',
%       `abcorr', `obsrvr', `dref' or `dvec', is undefined, an error
%       is signaled by the Matlab error handling system.
%
%   16) If any of the input arguments, `method', `target', `et',
%       `abcorr', `obsrvr', `dref' or `dvec', is not of the expected
%       type, or it does not have the expected dimensions and size, an
%       error is signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPK, PCK, and frame kernels must be loaded by the
%   calling program before this routine is called.  CK, SCLK, and
%   IK kernels may be required as well.
%
%   The following data are required:
%
%   -  SPK data: ephemeris data for target and observer must be
%      loaded. If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%      Typically ephemeris data are made available by loading one
%      or more SPK files via cspice_furnsh.
%
%   -  PCK data: if the computation method is specified as
%      "Ellipsoid," triaxial radii for the target body must be
%      loaded into the kernel pool. Typically this is done by
%      loading a text PCK file via cspice_furnsh.
%
%   -  Further PCK data: rotation data for the target body must
%      be loaded. These may be provided in a text or binary PCK
%      file.
%
%   -  Frame data: if a frame definition is required to convert
%      the observer and target states to the body-fixed frame of
%      the target, that definition must be available in the kernel
%      pool. Similarly, the frame definition required to map
%      between the frame designated by `dref' and the target
%      body-fixed frame must be available. Typically the
%      definitions of frames not already built-in to SPICE are
%      supplied by loading a frame kernel.
%
%   The following data may be required:
%
%   -  CK data: if the frame to which `dref' refers is fixed to a
%      spacecraft instrument or structure, at least one CK file
%      will be needed to permit transformation of vectors between
%      that frame and both J2000 and the target body-fixed frame.
%
%   -  SCLK data: if a CK file is needed, an associated SCLK
%      kernel is required to enable conversion between encoded SCLK
%      (used to time-tag CK data) and barycentric dynamical time
%      (TDB).
%
%   -  IK data: one or more I-kernels may be required to enable
%      transformation of vectors from an instrument-fixed frame to
%      a spacecraft-fixed frame whose attitude is given by a
%      C-kernel.
%
%
%   In all cases, kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   1)  A cautionary note: if aberration corrections are used, and if
%       `dref' is the target body-fixed frame, the epoch at which that
%       frame is evaluated is offset from `et' by the light time between
%       the observer and the *center* of the target body. This light
%       time normally will differ from the light time between the
%       observer and intercept point. Consequently the orientation of
%       the target body-fixed frame at `trgepc' will not match that of
%       the target body-fixed frame at the epoch associated with `dref'.
%       As a result, various derived quantities may not be as
%       expected: for example, `obspos' would not be the inverse of the
%       aberration-corrected position of the target as seen by the
%       observer.
%
%       In many applications the errors arising from this frame
%       discrepancy may be insignificant; however a safe approach is
%       to always use as `dref' a frame other than the target body-fixed
%       frame.
%
%-Required_Reading
%
%   MICE.REQ
%   FRAMES.REQ
%   NAIF_IDS.REQ
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
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 26-OCT-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Modified example
%       to load the required data via a meta-kernel and to map output of
%       example in sister function. Added example's meta-kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Updated -Required_Reading entries.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 12-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Index line now states that this routine is deprecated.
%
%   -Mice Version 1.0.0, 03-JAN-2006 (EDW)
%
%-Index_Entries
%
%   DEPRECATED surface intercept point
%
%-&

function [surf] = mice_srfxpt( method, target, et, abcorr, obsrvr, dref, dvec)

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

         error ( ['Usage: [_surf_ ] = '                                    ...
                  'mice_srfxpt( `method`, `target`,  '                     ...
                  '_et_, `abcorr`, `obsrvr`, '                             ...
                  '`dref`, dvec(6))']  )

   end

   %
   % Call the MEX library. The "_s" suffix indicates a structure type
   % return argument.
   %
   try
      [surf] = mice( 'srfxpt_s', method, target, et,                       ...
                     abcorr,     obsrvr, dref,   dvec );
   catch spiceerr
      rethrow(spiceerr)
   end





