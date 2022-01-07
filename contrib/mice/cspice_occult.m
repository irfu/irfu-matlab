%-Abstract
%
%   CSPICE_OCCULT determines the occultation condition (not occulted,
%   partially occulted, etc.) of one target relative to another target as
%   seen by an observer at a given time.
%
%   The surfaces of the target bodies may be represented by triaxial
%   ellipsoids, points, or by topographic data provided by DSK files.
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
%      targ1    the name of the first target body.
%
%               [1,c1] = size(targ1); char = class(targ1)
%
%                  or
%
%               [1,1] = size(targ1); cell = class(targ1)
%
%               Both object names and NAIF IDs are accepted. For example,
%               both 'Moon' and '301' are accepted.
%
%      shape1   a string indicating the geometric model used to represent the
%               shape of the first target body.
%
%               [1,c2] = size(shape1); char = class(shape1)
%
%                  or
%
%               [1,1] = size(shape1); cell = class(shape1)
%
%               The supported options are:
%
%                  'ELLIPSOID'
%
%                      Use a triaxial ellipsoid model with radius
%                      values provided via the kernel pool. A kernel
%                      variable having a name of the form
%
%                         'BODYnnn_RADII'
%
%                      where nnn represents the NAIF integer code
%                      associated with the body, must be present in
%                      the kernel pool. This variable must be
%                      associated with three numeric values giving the
%                      lengths of the ellipsoid's X, Y, and Z
%                      semi-axes.
%
%                  'POINT'
%
%                      Treat the body as a single point. When a point
%                      target is specified, the occultation conditions
%                      can only be total, annular, or none.
%
%                  'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
%
%                      Use topographic data provided by DSK files to
%                      model the body's shape. These data must be
%                      provided by loaded DSK files.
%
%                      The surface list specification is optional. The
%                      syntax of the list is
%
%                         <surface 1> [, <surface 2>...]
%
%                      If present, it indicates that data only for the
%                      listed surfaces are to be used; however, data
%                      need not be available for all surfaces in the
%                      list. If absent, loaded DSK data for any surface
%                      associated with the target body are used.
%
%                      The surface list may contain surface names or
%                      surface ID codes. Names containing blanks must
%                      be delimited by double quotes, for example
%
%                         SURFACES = "Mars MEGDR 128 PIXEL/DEG"
%
%                      If multiple surfaces are specified, their names
%                      or IDs must be separated by commas.
%
%                      See the -Particulars section below for details
%                      concerning use of DSK data.
%
%               The combinations of the shapes of the target bodies
%               `targ1' and `targ2' must be one of:
%
%                  One ELLIPSOID, one POINT
%                  Two ELLIPSOIDs
%                  One DSK, one POINT
%
%               Case and leading or trailing blanks are not
%               significant in the string `shape1'.
%
%      frame1   the name of the body-fixed, body-centered reference frame
%               associated with the first target body.
%
%               [1,c3] = size(frame1); char = class(frame1)
%
%                  or
%
%               [1,1] = size(frame1); cell = class(frame1)
%
%               Examples of such names are 'IAU_SATURN' (for Saturn) and
%               'ITRF93' (for the Earth).
%
%               If the first target body is modeled as a point, `frame1'
%               should be left blank (Ex: ' ').
%
%               Case and leading or trailing blanks bracketing a
%               non-blank frame name are not significant in the string.
%
%      targ2    the name of the second target body.
%
%               [1,c4] = size(targ2); char = class(targ2)
%
%                  or
%
%               [1,1] = size(targ2); cell = class(targ2)
%
%               See the description of `targ1' above for more details.
%
%      shape2   the shape specification for the body designated by `targ2'.
%
%               [1,c5] = size(shape2); char = class(shape2)
%
%                  or
%
%               [1,1] = size(shape2); cell = class(shape2)
%
%               See the description of `shape1' above for details.
%
%      frame2   the name of the body-fixed, body-centered reference frame
%               associated with the second target body.
%
%               [1,c6] = size(frame2); char = class(frame2)
%
%                  or
%
%               [1,1] = size(frame2); cell = class(frame2)
%
%               See the description of `frame1' above for more details.
%
%      abcorr   indicates the aberration corrections to be applied to the
%               state of each target body to account for one-way light time.
%
%               [1,c7] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               Stellar aberration corrections are ignored if specified,
%               since these corrections don't improve the accuracy of the
%               occultation determination.
%
%               See the header of the Mice routine cspice_spkezr for a
%               detailed description of the aberration correction
%               options. For convenience, the options supported by
%               this routine are listed below:
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       "Reception" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'CN'       "Reception" case: converged
%                             Newtonian light time correction.
%
%                  'XLT'      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'XCN'      "Transmission" case: converged
%                             Newtonian light time correction.
%
%               Case and blanks are not significant in the string
%               `abcorr'.
%
%      obsrvr   the name of the body from which the occultation is observed.
%
%               [1,c8] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               See the description of `targ1' above for more details.
%
%      et       the observation time(s) in seconds past the J2000 epoch.
%
%               [1,n] = size(et); double = class(et)
%
%   the call:
%
%      [ocltid] = cspice_occult( targ1,  shape1, frame1,                   ...
%                                targ2,  shape2, frame2,                   ...
%                                abcorr, obsrvr, et      )
%
%   returns:
%
%      ocltid   an integer occultation code(s) indicating the geometric
%               relationship of the three bodies.
%
%               [1,n] = size(ocltid); int32 = class(ocltid)
%
%               The meaning of the sign of `ocltid' is given below.
%
%                   Code sign          Meaning
%                   ---------          ------------------------------
%                      > 0             The second target is
%                                      partially or fully occulted
%                                      by the first.
%
%                      < 0             The first target is
%                                      partially of fully
%                                      occulted by the second.
%
%                      = 0             No occultation.
%
%               Possible `ocltid' values and meanings are given below. The
%               variable names indicate the type of occultation and which
%               target is in the back. For example, SPICE_OCCULT_TOTAL1
%               represents a total occultation in which the first target is in
%               the back of (or is occulted by) the second target.
%
%               When the target shapes are DSK and POINT, the only
%               possible occultation conditions are total, annular,
%               or none.
%
%                  Name                 Code   Meaning
%                  -------------------  -----  -----------------------
%                  SPICE_OCCULT_TOTAL1   -3    Total occultation of
%                                              first target by second.
%
%                  SPICE_OCCULT_ANNLR1   -2    Annular occultation of
%                                              first target by second.
%                                              The second target does
%                                              not block the limb of
%                                              the first.
%
%                  SPICE_OCCULT_PARTL1   -1    Partial occultation of
%                                              first target by second
%                                              target.
%
%                  SPICE_OCCULT_NOOCC     0    No occultation or
%                                              transit: both objects
%                                              are completely visible
%                                              to the observer.
%
%                  SPICE_OCCULT_PARTL2    1    Partial occultation of
%                                              second target by first
%                                              target.
%
%                  SPICE_OCCULT_ANNLR2    2    Annular occultation of
%                                              second target by first.
%
%                  SPICE_OCCULT_TOTAL2    3    Total occultation of
%                                              second target by first.
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
%   1) Find whether MRO is occulted by Mars as seen by
%      the DSS-13 ground station at a few specific
%      times.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: occult_ex1.tm
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
%            File name                       Contents
%            ---------                       --------
%            de421.bsp                       Planetary ephemeris
%            pck00010.tpc                    Planet orientation and
%                                            radii
%            naif0010.tls                    Leapseconds
%            earth_latest_high_prec.bpc      Earth latest binary PCK
%            earthstns_itrf93_050714.bsp     DSN station SPK
%            mro_psp22.bsp                   MRO ephemeris
%            earth_topo_050714.tf            Topocentric reference
%                                            frames for DSN
%                                            stations
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'mro_psp22.bsp',
%                                'earthstns_itrf93_050714.bsp',
%                                'earth_latest_high_prec.bpc',
%                                'pck00010.tpc',
%                                'naif0010.tls',
%                                'earth_topo_050714.tf' )
%         \begintext
%
%         End of meta-kernel.
%
%
%      Example code begins here.
%
%
%      function occult_ex1()
%
%         %
%         %   MiceUser is a file that makes occultation-specific
%         %   variables global. You must call MiceUser to have access to
%         %   these variable names. These variables are defined so you don't
%         %   need to remember what the integer codes mean that the
%         %   occultation routine returns. For more information, please see
%         %   MiceUser.m and MiceOccult.m.
%         %
%         MiceUser;
%
%         targ1    = 'MRO';
%         shape1   = 'point';
%         targ2    = 'Mars';
%         shape2   = 'ellipsoid';
%         observer = 'DSS-13';
%         dt = 1000;
%         out_form = '%s %s %s %s wrt %s\n';
%         out_char = ['totally occulted by  ';
%                     'transited by         ';
%                     'partially occulted by';
%                     'not occulted by      '];
%
%         %
%         %   Load the meta kernel.
%         %
%         cspice_furnsh ( 'occult_ex1.tm' );
%
%         et_start = cspice_str2et ( '2012-jan-5 1:15:00 UTC' );
%         et_stop  = cspice_str2et ( '2012-jan-5 2:50:00 UTC' );
%
%         for et = et_start : dt : et_stop
%
%             %
%             %   Calculate the type of occultation that
%             %   corresponds to time ET.
%             %
%             ocltid = cspice_occult ( targ1, shape1, ' ',               ...
%                                      targ2, shape2, 'iau_mars',        ...
%                                      'cn', observer, et );
%             %
%             %   Output the results.
%             %
%             time = cspice_timout ( et, 'YYYY-MM-DD HR:MN ::UTC-8');
%
%             %
%             %   Remember: You must call 'MiceUser' before
%             %   using the parameters like 'SPICE_OCCULT_TOTAL1' in
%             %   the case statements below.
%             %
%             switch ocltid
%                 case SPICE_OCCULT_TOTAL1
%                     fprintf (out_form, time, targ1, out_char(1,:),     ...
%                              targ2,  observer )
%                 case SPICE_OCCULT_ANNLR1
%                     fprintf (out_form, time, targ1, out_char(2,:),     ...
%                              targ2,  observer )
%                 case SPICE_OCCULT_PARTL1
%                     fprintf (out_form, time, targ1, out_char(3,:),     ...
%                              targ2,  observer )
%                 case SPICE_OCCULT_NOOCC
%                     fprintf (out_form, time, targ1, out_char(4,:),     ...
%                              targ2,  observer )
%                 case SPICE_OCCULT_PARTL2
%                     fprintf (out_form, time, targ2, out_char(3,:),     ...
%                              targ1,  observer )
%                 case SPICE_OCCULT_ANNLR2
%                     fprintf (out_form, time, targ2, out_char(2,:),     ...
%                              targ1,  observer )
%                 case SPICE_OCCULT_TOTAL2
%                     fprintf (out_form, time, targ2, out_char(1,:),     ...
%                              targ1,  observer )
%                 otherwise
%                     fprintf ( 'Bad occultation code: %d\n', ocltid )
%             end
%
%         end
%
%         %
%         %   Unload kernels
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      2012-01-04 17:15 Mars transited by          MRO wrt DSS-13
%      2012-01-04 17:31 MRO not occulted by       Mars wrt DSS-13
%      2012-01-04 17:48 MRO totally occulted by   Mars wrt DSS-13
%      2012-01-04 18:04 MRO totally occulted by   Mars wrt DSS-13
%      2012-01-04 18:21 MRO not occulted by       Mars wrt DSS-13
%      2012-01-04 18:38 Mars transited by          MRO wrt DSS-13
%
%
%-Particulars
%
%   For many purposes, modeling extended bodies as triaxial
%   ellipsoids is adequate for determining whether one body is
%   occulted by another as seen from a specified observer.
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
%      both represent a surface as a set of triangular plates, the
%      union of those sets of plates is considered to represent the
%      surface.
%
%      Currently only unprioritized data selection is supported.
%      Because prioritized data selection may be the default behavior
%      in a later version of the routine, the UNPRIORITIZED keyword is
%      required in the `shape1' and `shape2' arguments.
%
%
%      Syntax of the shape input arguments for the DSK case
%      ----------------------------------------------------
%
%      The keywords and surface list in the target shape arguments
%      `shape1' and `shape2' are called "clauses." The clauses may
%      appear in any order, for example
%
%         'DSK/<surface list>/UNPRIORITIZED'
%         'DSK/UNPRIORITIZED/<surface list>'
%         'UNPRIORITIZED/<surface list>/DSK'
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
%         'SURFACES = '
%
%      followed by a comma-separated list of one or more surface
%      identifiers. The identifiers may be names or integer codes in
%      string format. For example, suppose we have the surface
%      names and corresponding ID codes shown below:
%
%         Surface Name                              ID code
%         ------------                              -------
%         "Mars MEGDR 128 PIXEL/DEG"                1
%         "Mars MEGDR 64 PIXEL/DEG"                 2
%         "Mars_MRO_HIRISE"                         3
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
%      An example of a shape argument that could be constructed
%      using one of the surface lists above is
%
%         'DSK/UNPRIORITIZED/SURFACES = "Mars MEGDR 64 PIXEL/DEG", 3'
%
%-Exceptions
%
%   1)  If the target or observer body names input by the user are
%       not recognized, an error is signaled by a routine in
%       the call tree of this routine.
%
%   2)  If the input shapes are not accepted, an error is signaled by
%       a routine in the call tree of this routine.
%
%   3)  If both input shapes are points, an error is signaled by a
%       routine in the call tree of this routine.
%
%   4)  If the radii of a target body modeled as an ellipsoid cannot
%       be determined by searching the kernel pool for a kernel
%       variable having a name of the form
%
%          'BODYnnn_RADII'
%
%       where nnn represents the NAIF integer code associated with
%       the body, an error is signaled by a routine in the
%       call tree of this routine.
%
%   5)  If any of the target or observer bodies (targ1, targ2, or
%       obsrvr) are the same, an error is signaled
%       by a routine in the call tree of this routine.
%
%   6)  If the loaded kernels provide insufficient data to compute any
%       required state vector, an error is signaled by a routine in
%       the call tree of this routine.
%
%   7)  If an error occurs while reading an SPK or other kernel,
%       the error is signaled by a routine in the call tree
%       of this routine.
%
%   8)  If the aberration correction specification `abcorr' is invalid,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   9)  If either `shape1' or `shape2' specifies that the target surface
%       is represented by DSK data, and no DSK files are loaded for
%       the specified target, an error is signaled by a routine in
%       the call tree of this routine.
%
%   10) If either `shape1' or `shape2' specifies that the target surface
%       is represented by DSK data, but the shape specification is
%       invalid, an error is signaled by a routine in the call tree
%       of this routine.
%
%   11) If any of the input arguments, `targ1', `shape1', `frame1',
%       `targ2', `shape2', `frame2', `abcorr', `obsrvr' or `et', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   12) If any of the input arguments, `targ1', `shape1', `frame1',
%       `targ2', `shape2', `frame2', `abcorr', `obsrvr' or `et', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   Appropriate SPICE kernels must be loaded by the calling program
%   before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: the calling application must load ephemeris data
%      for the target, source and observer that cover the time
%      instant specified by the argument `et'. If aberration
%      corrections are used, the states of the target bodies and of
%      the observer relative to the solar system barycenter must be
%      calculable from the available ephemeris data. Typically
%      ephemeris data
%      are made available by loading one or more SPK files via
%      cspice_furnsh.
%
%   -  PCK data: bodies modeled as triaxial ellipsoids must have
%      semi-axis lengths provided by variables in the kernel pool.
%      Typically these data are made available by loading a text
%      PCK file via cspice_furnsh.
%
%   -  FK data: if either of the reference frames designated by
%      `frame1' or `frame2' are not built in to the SPICE system,
%      one or more FKs specifying these frames must be loaded.
%
%   The following data may be required:
%
%   -  DSK data: if either `shape1' or `shape2' indicates that DSK
%      data are to be used, DSK files containing topographic data
%      for the target body must be loaded. If a surface list is
%      specified, data for at least one of the listed surfaces must
%      be loaded.
%
%   -  Surface name-ID associations: if surface names are specified
%      in `shape1' or `shape2', the association of these names with
%      their corresponding surface ID codes must be established by
%      assignments of the kernel variables
%
%         NAIF_SURFACE_NAME
%         NAIF_SURFACE_CODE
%         NAIF_SURFACE_BODY
%
%      Normally these associations are made by loading a text
%      kernel containing the necessary assignments. An example
%      of such a set of assignments is
%
%         NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
%         NAIF_SURFACE_CODE += 1
%         NAIF_SURFACE_BODY += 499
%
%   -  CK data: either of the body-fixed frames to which `frame1' or
%      `frame2' refer might be a CK frame. If so, at least one CK
%      file will be needed to permit transformation of vectors
%      between that frame and the J2000 frame.
%
%   -  SCLK data: if a CK file is needed, an associated SCLK
%      kernel is required to enable conversion between encoded SCLK
%      (used to time-tag CK data) and barycentric dynamical time
%      (TDB).
%
%   Kernel data are normally loaded once per program run, NOT every
%   time this routine is called.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   DSK.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 2.1.0, 03-NOV-2021 (EDW) (JDR)
%
%       Changed argument names "target1", "target2", "observer", "time" and
%       "occult_code" to "targ1", "targ2", "obsrvr", "et" and "ocltid" for
%       consistency with other functions.
%
%       Edited the header to comply with NAIF standard. Added square brackets
%       to output argument in function declaration.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 2.0.0, 04-APR-2017 (EDW) (NJB)
%
%       Header update to reflect support for use of DSKs.
%
%   -Mice Version 1.0.0, 14-NOV-2013 (SCK)
%
%-Index_Entries
%
%   occultation type at a specified time
%
%-&

function [ocltid] = cspice_occult ( targ1, shape1,   frame1, ...
                                    targ2, shape2,   frame2, ...
                                    abcorr,  obsrvr, et )

   switch nargin
      case 9

         targ1  = zzmice_str(targ1);
         shape1 = zzmice_str(shape1);
         frame1 = zzmice_str(frame1);
         targ2  = zzmice_str(targ2);
         shape2 = zzmice_str(shape2);
         frame2 = zzmice_str(frame2);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         et     = zzmice_dp(et);

      otherwise

         error ( ['Usage: [_output_state_] = ' ...
                  'cspice_occult( `targ1`, `shape1`, `frame1`,' ...
                  '`targ2`, `shape2`, `frame2`, `abcorr`, '...
                  '`obsrvr`, _et_)'] )

   end

   %
   % Call the MEX library. An "_s" suffix indicates a structure type
   % return argument (not present in this case).
   %
   try
      [ocltid] = mice('occult_c', targ1, shape1,   frame1, ...
                                  targ2, shape2,   frame2, ...
                                  abcorr,  obsrvr, et );

   catch spiceerr
      rethrow(spiceerr)
   end
