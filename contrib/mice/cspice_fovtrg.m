%-Abstract
%
%   CSPICE_FOVTRG determines if a specified ephemeris object is within
%   the field-of-view (FOV) of a specified instrument at a given time.
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
%               determine if the target is visible with respect to the
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
%      target   the name of the target body.
%
%               [1,c2] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               This routine determines if the target body appears in the
%               instrument's field of view.
%
%               Both object names and NAIF IDs are accepted. For
%               example, both 'Moon' and '301' are accepted. Case and
%               leading or trailing blanks are not significant in the
%               string.
%
%      tshape   a string indicating the geometric model used to
%               represent the shape of the target body.
%
%               [1,c3] = size(tshape); char = class(tshape)
%
%                  or
%
%               [1,1] = size(tshape); cell = class(tshape)
%
%               The supported options are:
%
%                  'ELLIPSOID'     Use a triaxial ellipsoid model,
%                                  with radius values provided via the
%                                  kernel pool. A kernel variable
%                                  having a name of the form
%
%                                     'BODYnnn_RADII'
%
%                                  where nnn represents the NAIF
%                                  integer code associated with the
%                                  body, must be present in the kernel
%                                  pool. This variable must be
%                                  associated with three numeric
%                                  values giving the lengths of the
%                                  ellipsoid's X, Y, and Z semi-axes.
%
%                  'POINT'         Treat the body as a single point.
%
%               Case and leading or trailing blanks are not
%               significant in the string.
%
%      tframe   the name of the body-fixed, body-centered reference
%               frame associated with the target body.
%
%               [1,c4] = size(tframe); char = class(tframe)
%
%                  or
%
%               [1,1] = size(tframe); cell = class(tframe)
%
%               Examples of such names are 'IAU_SATURN' (for Saturn) and
%               'ITRF93' (for Earth).
%
%               If the target body is modeled as a point, `tframe'
%               is ignored and should be left blank. (Ex: ' ').
%
%               Case and leading or trailing blanks bracketing a
%               non-blank frame name are not significant in the string.
%
%      abcorr   indicates the aberration corrections to be applied
%               when computing the target's position and orientation.
%
%               [1,c5] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               For remote sensing applications, where the apparent
%               position and orientation of the target seen by the
%               observer are desired, normally either of the
%               corrections:
%
%                  'LT+S'
%                  'CN+S'
%
%               should be used. These and the other supported options
%               are described below.
%
%               Supported aberration correction options for
%               observation (the case where radiation is received by
%               observer at `et') are:
%
%                  'NONE'         No correction.
%                  'LT'           Light time only
%                  'LT+S'         Light time and stellar aberration.
%                  'CN'           Converged Newtonian (CN) light time.
%                  'CN+S'         CN light time and stellar aberration.
%
%               Supported aberration correction options for
%               transmission (the case where radiation is emitted from
%               observer at `et') are:
%
%                  'XLT'          Light time only.
%                  'XLT+S'        Light time and stellar aberration.
%                  'XCN'          Converged Newtonian (CN) light time.
%                  'XCN+S'        CN light time and stellar aberration.
%
%               Case, leading and trailing blanks are not significant
%               in the string.
%
%      obsrvr   the name of the body from which the target is
%               observed.
%
%               [1,c6] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               The instrument `inst' is treated as if it were co-located
%               with the observer.
%
%               Both object names and NAIF IDs are accepted. For
%               example, both 'CASSINI' and '-82' are accepted. Case and
%               leading or trailing blanks are not significant in the
%               string.
%
%      et       the observation time(s) in seconds past the J2000
%               epoch.
%
%               [1,n] = size(et); double = class(et)
%
%   the call:
%
%      [visibl] = cspice_fovtrg( inst,   target, tshape, ...
%                                tframe, abcorr, obsrvr, et )
%
%   returns:
%
%      visibl   true if `target' is fully or partially in the
%               field-of-view of `inst' at the time `et'.
%
%               [1,n] = size(visibl); logical = class(visibl)
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
%   1) A spectacular picture was taken by Cassini's
%      narrow-angle camera on Oct. 6, 2010 that shows
%      six of Saturn's moons. Let's verify that the moons
%      in the picture are Epimetheus, Atlas, Daphnis, Pan,
%      Janus, and Enceladus.
%
%      To see this picture, visit:
%      http://photojournal.jpl.nasa.gov/catalog/PIA12741
%      or go to the PDS Image Node's Image Atlas at
%      http://pds-imaging.jpl.nasa.gov/search/search.html.
%      Select Cassini as the mission, ISS as the instrument,
%      and enter 1_N1665078907.122 as the Product ID in the
%      Product tab. Note: these directions may change as the
%      PDS Imaging Node changes.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: fovtrg_ex1.tm
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
%            File name                     Contents
%            ---------                     --------
%            naif0010.tls                  Leapseconds
%            cpck*.tpc                     Satellite orientation and
%                                          radii
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            cas_rocks_v18.tf              FK for small satellites
%                                          around Saturn
%            cas_v40.tf                    Cassini FK
%            cas_iss_v10.ti                Cassini ISS IK
%            cas00149.tsc                  Cassini SCLK
%            *.bsp                         Ephemeris for Cassini,
%                                          planets, and satellites
%            10279_10284ra.bc              Orientation for Cassini
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0010.tls'
%                                'cpck14Oct2010.tpc'
%                                'cpck_rock_21Jan2011_merged.tpc'
%                                'pck00010.tpc'
%                                'cas_rocks_v18.tf'
%                                'cas_v40.tf'
%                                'cas_iss_v10.ti'
%                                'cas00149.tsc'
%                                '110317AP_RE_90165_18018.bsp'
%                                '110120BP_IRRE_00256_25017.bsp'
%                                '101210R_SCPSE_10256_10302.bsp'
%                                '10279_10284ra.bc'              )
%
%         \begintext
%
%         For project meta-kernels similar to the one shown
%         here, please see the CASSINI SPICE PDS archive.
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function fovtrg_ex1()
%
%         %
%         % Load the meta kernel.
%         %
%         cspice_furnsh( 'fovtrg_ex1.tm' );
%
%         %
%         % Retrieve Cassini's NAIF ID.
%         %
%         [cassini_id, found] = cspice_bodn2c( 'CASSINI' );
%
%         if (~found)
%            fprintf( 'Could not find ID code for Cassini.\n' )
%            return
%         end
%
%         %
%         % Convert the image tag SCLK to ET.
%         %
%         et = cspice_scs2e( cassini_id, '1665078907.122' );
%
%         %
%         % Convert `et' to a string format for the output.
%         %
%         time_format = 'YYYY-MON-DD HR:MN:SC.###::TDB (TDB)';
%         time = cspice_timout( et, time_format );
%
%         %
%         % Search through all of Saturn's moons to see if each
%         % satellite was in the ISS NAC's field-of-view at
%         % the image time. We're going to take advantage of the
%         % fact that all Saturn's moons have a NAIF ID of 6xx.
%         %
%         fprintf( 'At time %s the following were\n', time     )
%         fprintf( 'in the field of view of CASSINI_ISS_NAC\n' )
%         for body_id = 600:699
%            %
%            % Check to see if the `body_id' has a translation.
%            %
%            [body, found] = cspice_bodc2n( body_id );
%
%            if (found)
%               %
%               % Check to see if a body-fixed frame for this ID exists.
%               % If the frame is not in the kernel pool, we cannot
%               % perform the visibility test. The main cause of a
%               % failure is a missing kernel.
%               %
%               [frame_code, frame_name, found] = cspice_cidfrm( body_id );
%
%               if (found)
%                  %
%                  % Is this body in the field-of-view of Cassini's
%                  % ISS narrow-angle camera?
%                  %
%                  visibl = cspice_fovtrg( 'CASSINI_ISS_NAC', body,        ...
%                                          'Ellipsoid',       frame_name,  ...
%                                          'CN+S',            'CASSINI', et );
%                  if ( visibl )
%                     fprintf( '  %s\n', body )
%                  end
%               end
%            end
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      At time 2010-OCT-06 17:09:45.346 (TDB) the following were
%      in the field of view of CASSINI_ISS_NAC
%        ENCELADUS
%        JANUS
%        EPIMETHEUS
%        ATLAS
%        PAN
%        DAPHNIS
%        ANTHE
%
%
%      Note: there were actually 7 of Saturn's satellites in the
%      field-of-view of Cassini's narrow-angle camera. However, Anthe
%      is very small and was probably obscured by other objects or
%      shadow.
%
%-Particulars
%
%   To treat the target as a ray rather than as an ephemeris object,
%   use the higher-level Mice routine cspice_fovray. cspice_fovray may be used
%   to determine if distant target objects such as stars are visible
%   in an instrument's FOV at a given time, as long as the direction
%   from the observer to the target can be modeled as a ray.
%
%-Exceptions
%
%   1)  If the name of either the target or observer cannot be
%       translated to a NAIF ID code, an error is signaled by
%       a routine in the call tree of this routine.
%
%   2)  If the specified aberration correction is an unrecognized
%       value, an error is signaled by a routine
%       in the call tree of this routine.
%
%   3)  If the radii of a target body modeled as an ellipsoid cannot
%       be determined by searching the kernel pool for a kernel
%       variable having a name of the form
%
%          'BODYnnn_RADII'
%
%       where nnn represents the NAIF integer code associated with
%       the body, an error is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If the target and observer bodies are the same, an error is
%       signaled by a routine in the call tree of this routine.
%
%   5)  If the body model specifier `tshape' is invalid, an error is
%       signaled by either this routine or a routine in the call tree
%       of this routine.
%
%   6)  If a target body-fixed reference frame associated with a
%       non-point target is not recognized, an error is signaled by a
%       routine in the call tree of this routine.
%
%   7)  If a target body-fixed reference frame is not centered at the
%       corresponding target body, an error is signaled by a routine
%       in the call tree of this routine.
%
%   8)  If the instrument name `inst' does not have a corresponding NAIF
%       ID code, an error is signaled by a routine in the call tree of
%       this routine.
%
%   9)  If the FOV parameters of the instrument are not present in
%       the kernel pool, an error is signaled by a routine
%       in the call tree of this routine.
%
%   10) If the FOV boundary has more than SPICE_GF_MAXVRT vertices, an error
%       is signaled by a routine in the call tree of this
%       routine.
%
%   11) If the instrument FOV shape is a polygon or rectangle, and
%       this routine cannot find a ray R emanating from the FOV vertex
%       such that maximum angular separation of R and any FOV boundary
%       vector is within the limit (pi/2)-SPICE_GF_MARGIN radians, an
%       error is signaled by a routine in the call tree of this routine.
%       If the FOV is any other shape, the same error check will be
%       applied with the instrument boresight vector serving the role
%       of R.
%
%   12) If the loaded kernels provide insufficient data to compute a
%       requested state vector, an error is signaled by a
%       routine in the call tree of this routine.
%
%   13) If an error occurs while reading an SPK or other kernel file,
%       the error is signaled by a routine in the call tree
%       of this routine.
%
%   14) If any of the input arguments, `inst', `target', `tshape',
%       `tframe', `abcorr', `obsrvr' or `et', is undefined, an error
%       is signaled by the Matlab error handling system.
%
%   15) If any of the input arguments, `inst', `target', `tshape',
%       `tframe', `abcorr', `obsrvr' or `et', is not of the expected
%       type, or it does not have the expected dimensions and size, an
%       error is signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: ephemeris data for target and observer that
%      describe the ephemerides of these objects at the time `et'.
%      If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%
%   -  Frame data: if a frame definition is required to convert
%      the observer and target states to the body-fixed frame of
%      the target, that definition must be available in the kernel
%      pool. Typically the definitions of frames not already
%      built-in to SPICE are supplied by loading a frame kernel.
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
%   -  PCK data: bodies modeled as triaxial ellipsoids must have
%      orientation data provided by variables in the kernel pool.
%
%      Bodies modeled as triaxial ellipsoids must have radii
%      lengths provided by variables in the kernel pool.
%
%   -  CK data: if the frame in which the instrument's FOV is
%      defined is fixed to a spacecraft, at least one CK file will
%      be needed to permit transformation of vectors between that
%      frame and both J2000 and the target body-fixed frame.
%
%   -  SCLK data: if a CK file is needed, an associated SCLK
%      kernel is required to enable conversion between encoded SCLK
%      (used to time-tag CK data) and barycentric dynamical time
%      (TDB).
%
%   Kernel data are normally loaded via cspice_furnsh once per program run,
%   NOT every time this routine is called.
%
%-Restrictions
%
%   1)  The reference frame associated with `inst' must be centered at
%       the observer or must be inertial. No check is done to ensure
%       this.
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
%       Changed argument names "instrument", "target_shape",
%       "target_frame", "observer" and "visible" to "inst", "tshape",
%       "tframe", "obsrvr" and "visibl" for consistency with other
%       routines.
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
%   -Mice Version 1.0.0, 16-FEB-2012 (SCK)
%
%-Index_Entries
%
%   Target in instrument FOV at specified time
%   Target in instrument field_of_view at specified time
%
%-&

function [visibl] = cspice_fovtrg( inst,   target, tshape, ...
                                   tframe, abcorr, obsrvr, et )

    switch nargin
        case 7

            inst   = zzmice_str(inst);
            target = zzmice_str(target);
            tshape = zzmice_str(tshape);
            tframe = zzmice_str(tframe);
            abcorr = zzmice_str(abcorr);
            obsrvr = zzmice_str(obsrvr);
            et     = zzmice_dp(et);

        otherwise

            error ( ['Usage: [_visibl_] = ' ...
                  'cspice_fovtrg( `inst`, `target`, '    ...
                  '`tshape`, `tframe`, `abcorr`, ' ...
                  '`obsrvr`, _et_]' ] )

   end


   %
   % Call the MEX library. An "_s" suffix indicates a structure type
   % return argument (not present in this case).
   %
   try
      [visibl] = mice('fovtrg_c', inst,   target, tshape, ...
                                  tframe, abcorr, obsrvr, et );
      [visibl] = zzmice_logical( visibl );
   catch spiceerr
      rethrow(spiceerr)
   end

