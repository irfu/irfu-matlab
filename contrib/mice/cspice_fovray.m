%-Abstract
%
%   CSPICE_FOVRAY determines if a specified ray is within the field-of-view
%   (FOV) of a specified instrument at a given time.
%
%-Disclaimer
%
%   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
%   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
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
%      inst     indicates the name of an instrument, such as a
%               spacecraft-mounted framing camera.
%
%               [1,c1] = size(inst); char = class(inst)
%
%                  or
%
%               [1,1] = size(inst); cell = class(inst)
%
%               The field of view (FOV) of the instrument will be used to
%               determine if the direction from the observer to a target,
%               represented as a ray, is visible with respect to the
%               instrument.
%
%               The position of the instrument `inst' is considered to
%               coincide with that of the ephemeris object `obsrvr' (see
%               description below).
%
%               The size of the instrument's FOV is constrained by the
%               following: There must be a vector A such that all of
%               the instrument's FOV boundary vectors have an angular
%               separation from A of less than (pi/2)-SPICE_GF_MARGIN radians
%               (see description below). For FOVs that are circular or
%               elliptical, the vector A is the boresight. For FOVs
%               that are rectangular or polygonal, the vector A is
%               calculated.
%
%               See the header of the Mice routine cspice_getfov for a
%               description of the required parameters associated with
%               an instrument.
%
%               Both object names and NAIF IDs are accepted. For
%               example, both 'CASSINI_ISS_NAC' and '-82360' are
%               accepted. Case and leading or trailing blanks are not
%               significant in the string.
%
%      raydir   the direction vector defining a ray of interest.
%
%               [3,1] = size(raydir), double = class(raydir)
%
%               The ray emanates from the location of the ephemeris
%               object designated by the input argument `obsrvr' and
%               is expressed relative to the reference frame designated
%               by `rframe' (see description below).
%
%      rframe   the name of the reference frame associated with
%               the input ray's direction vector `raydir'.
%
%               [1,c2] = size(rframe), char = class(rframe)
%
%                  or
%
%               [1,1] = size(rframe); cell = class(rframe)
%
%               Note: `rframe' does not need to be the instrument's reference
%               frame.
%
%               Since light time corrections are not supported for
%               rays, the orientation of the frame is always evaluated
%               at the epoch associated with the observer, as opposed
%               to the epoch associated with the light-time corrected
%               position of the frame center.
%
%               Case, leading and trailing blanks are not significant
%               in the string.
%
%      abcorr   indicates the aberration corrections to be applied
%               when computing the ray's direction.
%
%               [1,c3] = size(abcorr), char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               The supported aberration correction options are:
%
%                  'NONE'          No correction.
%                  'S'             Stellar aberration correction,
%                                  reception case.
%                  'XS'            Stellar aberration correction,
%                                  transmission case.
%
%               For detailed information, see the geometry finder
%               required reading, gf.req.
%
%               Case, leading and trailing blanks are not significant
%               in the string.
%
%      obsrvr   the name of the body from which the target
%               represented by `raydir' is observed.
%
%               [1,c4] = size(obsrvr), char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               The instrument designated by `inst' is treated as if it were
%               co-located with the observer.
%
%               Both object names and NAIF IDs are accepted. For
%               example, both 'CASSINI' and '-82' are accepted. Case
%               and leading or trailing blanks are not significant in
%               the string.
%
%      et       the observation time(s) in seconds past the J2000
%               epoch.
%
%               [1,n] = size(et), double = class(et)
%
%   the call:
%
%      [visibl] = cspice_fovray( inst, raydir, rframe, abcorr, obsrvr, et )
%
%   returns:
%
%      visibl   true if the ray is "visible", or in the
%               field-of-view, of `inst' at the time `et'.
%
%               [1,n] = size(visibl), logical = class(visibl)
%
%               Otherwise, `visibl' is false.
%
%               `visibl' returns with the same vectorization measure, N,
%               as `et'.
%
%-Parameters
%
%   SPICE_GF_MAXVRT
%
%               is the maximum number of vertices that may be used
%               to define the boundary of the specified instrument's
%               field of view.
%
%   SPICE_GF_MARGIN
%
%               is a small positive number used to constrain the
%               orientation of the boundary vectors of polygonal
%               FOVs. Such FOVs must satisfy the following constraints:
%
%                  1)  The boundary vectors must be contained
%                      within a right circular cone of angular radius
%                      less than than (pi/2) - SPICE_GF_MARGIN radians;
%                      in other words, there must be a vector A such
%                      that all boundary vectors have angular
%                      separation from A of less than
%                      (pi/2)-SPICE_GF_MARGIN radians.
%
%                  2)  There must be a pair of boundary vectors U,
%                      V such that all other boundary vectors lie in
%                      the same half space bounded by the plane
%                      containing U and V. Furthermore, all other
%                      boundary vectors must have orthogonal
%                      projections onto a specific plane normal to this
%                      plane (the normal plane contains the angle
%                      bisector defined by U and V) such that the
%                      projections have angular separation of at least
%                      2*SPICE_GF_MARGIN radians from the plane spanned
%                      by U and V.
%
%               SPICE_GF_MARGIN is currently set to 1.e-12.
%
%   See Mice header file MiceGF.m for declarations and descriptions of
%   parameters used throughout the GF system.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%
%   1) The Cassini Ultraviolet Imaging Spectrograph (UVIS)
%      has been used to measure variations in starlight as
%      rings and moons occult Cassini's view of the stars.
%      One of these events happened at 2008-054T21:31:55.158 UTC.
%      Let's verify that Epsilon CMa (Adhara) was in the
%      Cassini UVIS field-of-view at the observation time.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: fovray_ex1.tm
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
%            File name                      Contents
%            ---------                      --------
%            naif0010.tls                   Leapseconds
%            cpck26Jan2007.tpc              Satellite orientation and
%                                           radii
%            cas00145.tsc                   Cassini SCLK
%            cas_v40.tf                     Cassini frames
%            cas_uvis_v06.ti                Cassini UVIS instrument
%            080428R_SCPSE_08045_08067.bsp  Merged spacecraft,
%                                           planetary, and satellite
%                                           ephemeris
%            08052_08057ra.bc               Orientation for Cassini
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0010.tls'
%                                'cpck26Jan2007.tpc'
%                                'cas00145.tsc'
%                                'cas_v40.tf'
%                                'cas_uvis_v06.ti'
%                                '080428R_SCPSE_08045_08067.bsp'
%                                '08052_08057ra.bc')
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function fovray_ex1()
%
%         %
%         %   Load the meta kernel.
%         %
%         cspice_furnsh ( 'fovray_ex1.tm' );
%
%         %
%         %   Convert the observation time to `et'.
%         %
%         et = cspice_str2et ( '2008-054T21:31:55.158' );
%
%         %
%         %   The variables `right_asc' and `dec' are the right ascension
%         %   and declination of Epsilon CMa in degrees.
%         %
%         right_asc = 104.656;
%         dec       = -28.972;
%
%         %
%         %   Create a unit direction vector pointing from Cassini
%         %   to the specified star. For details on corrections such
%         %   as parallax, please see the example in cspice_gfrfov.
%         %
%
%         raydir = cspice_radrec ( 1, right_asc*cspice_rpd, dec*cspice_rpd );
%
%         %
%         %   Is the star in the field-of-view of Cassini's UVIS?
%         %
%         visible = cspice_fovray ( 'CASSINI_UVIS_FUV_OCC', raydir,        ...
%                                   'J2000', 's', 'cassini', et );
%
%         %
%         %   Put the time in a specified format for output and
%         %   report the result.
%         %
%         time_output = cspice_timout ( et,                                ...
%                                   'YYYY-MON-DD HR:MN:SC.###::TDB (TDB)' );
%
%         if ( visible )
%             fprintf ( 'Epsilon CMa was visible from the Cassini\n' );
%             fprintf ( 'UVIS instrument at %s\n', time_output );
%         end
%
%         %
%         %   Unload kernels.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Epsilon CMa was visible from the Cassini
%      UVIS instrument at 2008-FEB-23 21:33:00.343 (TDB)
%
%
%-Particulars
%
%   To treat the target as an ephemeris object rather than a ray, use
%   the higher-level Mice routine cspice_fovtrg. cspice_fovtrg may be used to
%   determine if ephemeris objects such as Saturn are visible in an
%   instrument's FOV at a given time.
%
%-Exceptions
%
%   1)  If the observer's name cannot be mapped to a NAIF ID code, the
%       error SPICE(IDCODENOTFOUND) is signaled by a routine in the
%       call tree of this routine.
%
%   2)  If the aberration correction flag calls for light time
%       correction, the error SPICE(INVALIDOPTION) is signaled by a
%       routine in the call tree of this routine.
%
%   3)  If the ray's direction vector is zero, the error
%       SPICE(ZEROVECTOR) is signaled by a routine in the call tree of
%       this routine.
%
%   4)  If the instrument name `inst' does not have corresponding NAIF
%       ID code, an error is signaled by a routine in the call
%       tree of this routine.
%
%   5)  If the FOV parameters of the instrument are not present in
%       the kernel pool, an error is signaled by a routine
%       in the call tree of this routine.
%
%   6)  If the FOV boundary has more than SPICE_GF_MAXVRT vertices, an error
%       is signaled by a routine in the call tree of this
%       routine.
%
%   7)  If the instrument FOV shape is a polygon or rectangle, and
%       this routine cannot find a ray R emanating from the FOV vertex
%       such that maximum angular separation of R and any FOV boundary
%       vector is within the limit (pi/2)-SPICE_GF_MARGIN radians, an error
%       is signaled by a routine in the call tree of this routine. If the
%       FOV is any other shape, the same error check will be applied
%       with the instrument boresight vector serving the role of R.
%
%   8)  If the loaded kernels provide insufficient data to compute a
%       requested state vector, an error is signaled by a
%       routine in the call tree of this routine.
%
%   9)  If an error occurs while reading an SPK or other kernel file,
%       the error is signaled by a routine in the call tree
%       of this routine.
%
%   10) If any of the input arguments, `inst', `raydir', `rframe',
%       `abcorr', `obsrvr' or `et', is undefined, an error is signaled
%       by the Matlab error handling system.
%
%   11) If any of the input arguments, `inst', `raydir', `rframe',
%       `abcorr', `obsrvr' or `et', is not of the expected type, or it
%       does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: ephemeris data for the observer at the time
%      `et'. If aberration corrections are used, the state of the
%      observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%
%   -  Data defining the reference frame in which the instrument's
%      FOV is defined must be available in the kernel pool.
%      Additionally the name `inst' must be associated with an ID
%      code.
%
%   -  IK data: the kernel pool must contain data such that
%      the Mice routine cspice_getfov may be called to obtain
%      parameters for `inst'.
%
%   The following data may be required:
%
%   -  CK data: if the frame in which the instrument's FOV is
%      defined is fixed to a spacecraft, at least one CK file will
%      be needed to permit transformation of vectors between that
%      frame and the J2000 frame.
%
%   -  SCLK data: if a CK file is needed, an associated SCLK
%      kernel is required to enable conversion between encoded SCLK
%      (used to time-tag CK data) and barycentric dynamical time
%      (TDB).
%
%   -  Since the input ray direction may be expressed in any
%      frame, additional FKs, CKs, SCLK kernels, PCKs, and SPKs
%      may be required to map the direction to the J2000 frame.
%
%   Kernel data are normally loaded via cspice_furnsh once per program run,
%   NOT every time this routine is called.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   CK.REQ
%   FRAMES.REQ
%   KERNEL.REQ
%   MICE.REQ
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
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Changed the argument names "instrument", "ray_frame", "observer"
%       and "visible" to "inst", "rframe", "obsrvr" and "visibl" for
%       consistency with other routines.
%
%       Added square brackets to output argument.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Corrected
%       the value and changed the name of the parameter "MARGIN" to
%       "SPICE_GF_MARGIN".
%
%       Edited the header to comply with NAIF standard. Completed the
%       list of required readings.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.1, 13-APR-2015 (EDW)
%
%       Edit to correct typos in "Usage" string.
%
%   -Mice Version 1.0.0, 13-NOV-2013 (SCK) (EDW)
%
%-Index_Entries
%
%   Ray in instrument FOV at specified time
%   Ray in instrument field_of_view at specified time
%
%-&

function [visibl] = cspice_fovray( inst,   raydir, rframe,  ...
                                   abcorr, obsrvr, et )

    switch nargin
        case 6

            inst    = zzmice_str(inst);
            raydir  = zzmice_dp (raydir);
            rframe  = zzmice_str(rframe);
            abcorr  = zzmice_str(abcorr);
            obsrvr  = zzmice_str(obsrvr);
            et      = zzmice_dp(et);

        otherwise

            error ( ['Usage: [_visibl_] = '              ...
                  'cspice_xfmsta( `inst`, _raydir[6]_, ' ...
                  '`rframe`, `abcorr`, `obsrvr`, _et_]' ] )

   end


   %
   % Call the MEX library.
   %
   try
      [visibl] = mice('fovray_c', inst, raydir,   rframe, ...
                                  abcorr,     obsrvr, et );
      [visibl] = zzmice_logical( visibl );
   catch spiceerr
      rethrow(spiceerr)
   end

