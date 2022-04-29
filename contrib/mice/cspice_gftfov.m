%-Abstract
%
%   CSPICE_GFTFOV determines time intervals when a specified ephemeris
%   object intersects the space bounded by the field-of-view (FOV) of a
%   specified instrument.
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
%      inst     the string naming the instrument, such as a
%               spacecraft-mounted framing camera, the field of view
%               (FOV) of which is to be used for a target intersection
%               search: times when the specified target intersects the
%               region of space corresponding to the FOV are sought.
%
%               [1,c1] = size(inst); char = class(inst)
%
%               The position of the instrument designated by `inst' is
%               considered to coincide with that of the ephemeris
%               object designated by the input argument `obsrvr' (see
%               description below).
%
%               `inst' must have a corresponding NAIF ID and a frame
%               defined, as is normally done in a frame kernel. It
%               must also have an associated reference frame and a FOV
%               shape, boresight and boundary vertices (or reference
%               vector and reference angles) defined, as is usually
%               done in an instrument kernel.
%
%               See the header of the Mice routine cspice_getfov for a
%               description of the required parameters associated with
%               an instrument.
%
%      target   the string naming the `target' body, the appearances
%               of which in the specified instrument's field of view are
%               sought.
%
%               [1,c2] = size(target); char = class(target)
%
%               The body must be an ephemeris object.
%
%               Optionally, you may supply the integer NAIF ID code
%               for the body as a string. For example both 'MOON' and
%               '301' are legitimate strings that designate the Moon.
%
%               The `target' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%      tshape   the string naming the geometric model used to
%               represent the shape of the `target' body.
%
%               [1,c3] = size(tshape); char = class(tshape)
%
%               The supported options are:
%
%                  'ELLIPSOID'   Use a triaxial ellipsoid model,
%                                with radius values provided via the
%                                kernel pool. A kernel variable
%                                having a name of the form
%
%                                   'BODYnnn_RADII'
%
%                                where nnn represents the NAIF
%                                integer code associated with the
%                                body, must be present in the kernel
%                                pool. This variable must be
%                                associated with three numeric
%                                values giving the lengths of the
%                                ellipsoid's X, Y, and Z semi-axes.
%
%                  'POINT'       Treat the body as a single point.
%
%               The `tshape' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      tframe   the string naming the body-fixed, body-centered
%               reference frame associated with the target body.
%
%               [1,c4] = size(tframe); char = class(tframe)
%
%               Examples of such names are 'IAU_SATURN' (for Saturn) and
%               'ITRF93' (for the Earth).
%
%               If the target body is modeled as a point, `tframe'
%               is ignored and should be left blank.
%
%               The `tframe' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%      abcorr   the string indicating the aberration corrections to apply
%               to the state evaluations to account for one-way light time and
%               stellar aberration.
%
%               [1,c5] = size(abcorr); char = class(abcorr)
%
%               For remote sensing applications, where the apparent
%               position and orientation of the target seen by the
%               observer are desired, normally either of the
%               corrections
%
%                  'LT+S'
%                  'CN+S'
%
%               should be used. These and the other supported options
%               are described below.
%
%                 'NONE'      Apply no correction.
%
%               Supported aberration correction options for reception case
%               (radiation is received by observer at `et') are:
%
%                  'LT'       Correct for one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     Correct for one-way light time and stellar
%                             aberration using a Newtonian formulation.
%
%                  'CN'       Correct for one-way light time using a converged
%                             Newtonian light time correction.
%
%                  'CN+S'     Correct for one-way light time and stellar
%                             aberration using a converged Newtonian light
%                             time and stellar aberration corrections.
%
%               Supported aberration correction options for transmission case
%               (radiation is emitted from observer at ET) are:
%
%                  'XLT'      Correct for one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    Correct for one-way light time and stellar
%                             aberration using a Newtonian formulation.
%
%                  'XCN'      Correct for one-way light time using a converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    Correct for one-way light time and stellar
%                             aberration using a converged Newtonian light
%                             time and stellar aberration corrections.
%
%               For detailed information, see the geometry finder
%               required reading, gf.req.
%
%               The `abcorr' string lacks sensitivity to case, and to leading
%               and trailing blanks.
%
%      obsrvr   the string naming the body from which the target is
%               observed.
%
%               [1,c6] = size(abcorr); char = class(abcorr)
%
%               The instrument designated by `inst' is treated as if it were
%               co-located with the observer.
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the observer
%               is Earth.
%
%      step     the step size to use in the search to use in the search.
%
%               [1,1] = size(step); double = class(step)
%
%               `step' must be short enough for a search using step
%               to locate the time intervals where the specified
%               angular separation function is monotone increasing or
%               decreasing. However, `step' must not be *too* short, or
%               the search will take an unreasonable amount of time.
%
%               The choice of `step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               `step' has units of TDB seconds.
%
%      cnfine   the SPICE window that confines the time
%               period over which the specified search is conducted.
%
%               [2m,1] = size(cnfine); double = class(cnfine)
%
%               `cnfine' may consist of a single interval or a collection
%               of intervals.
%
%               In some cases the confinement window can be used to
%               greatly reduce the time period that must be searched
%               for the desired solution. See the -Particulars section
%               below for further discussion.
%
%      nintvls  the maximum number of intervals to return in `result'.
%
%               [1,1] = size(nintvls); int32 = class(nintvls)
%
%               Note: this value should equal at least the number of expected
%               intervals. Recall two double precision values define
%               an interval.
%
%   the call:
%
%      result = cspice_gftfov( inst,   target, tshape, tframe, abcorr,     ...
%                              obsrvr, step,   cnfine, nintvls )
%
%   returns:
%
%      result   the SPICE window of intervals, contained within the
%               confinement window `cnfine', on which the specified
%               constraint is satisfied.
%
%               [2n,1] = size(result); double = class(result)
%
%               If the search is for local extrema, or for absolute
%               extrema with `adjust' set to zero, then normally each
%               interval of `result' will be a singleton: the left and
%               right endpoints of each interval will be identical.
%
%               If no times within the confinement window satisfy the
%               constraint, `result' will return with cardinality zero.
%
%-Parameters
%
%   All parameters described here are declared in the Mice include file
%   MiceGF.m. See that file for parameter values.
%
%   SPICE_GF_CNVTOL
%
%               is the convergence tolerance used for finding endpoints
%               of the intervals comprising the result window.
%               SPICE_GF_CNVTOL is used to determine when binary searches
%               for roots should terminate: when a root is bracketed
%               within an interval of length SPICE_GF_CNVTOL, the root is
%               considered to have been found.
%
%               The accuracy, as opposed to precision, of roots found
%               by this routine depends on the accuracy of the input
%               data. In most cases, the accuracy of solutions will be
%               inferior to their precision.
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
%                  1)  The boundary vectors must be contained within
%                      a right circular cone of angular radius less
%                      than (pi/2) - SPICE_GF_MARGIN radians; in other
%                      words, there must be a vector A such that all
%                      boundary vectors have angular separation from
%                      A of less than (pi/2)-SPICE_GF_MARGIN radians.
%
%                  2)  There must be a pair of boundary vectors U, V
%                      such that all other boundary vectors lie in
%                      the same half space bounded by the plane
%                      containing U and V. Furthermore, all other
%                      boundary vectors must have orthogonal
%                      projections onto a plane normal to this plane
%                      such that the projections have angular
%                      separation of at least 2*SPICE_GF_MARGIN radians
%                      from the plane spanned by U and V.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Search for times when Saturn's satellite Phoebe is within
%      the FOV of the Cassini narrow angle camera (CASSINI_ISS_NAC).
%      To simplify the problem, restrict the search to a short time
%      period where continuous Cassini bus attitude data are
%      available.
%
%      Use a step size of 10 seconds to reduce chances of missing
%      short visibility events.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gftfov_ex1.tm
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
%            -----------------------------   ----------------------
%            naif0012.tls                    Leapseconds
%            pck00010.tpc                    Satellite orientation
%                                            and radii
%            041014R_SCPSE_01066_04199.bsp   CASSINI, planetary and
%                                            Saturn satellite
%                                            ephemeris
%            cas_v40.tf                      Cassini FK
%            04161_04164ra.bc                Cassini bus CK
%            cas00071.tsc                    Cassini SCLK kernel
%            cas_iss_v10.ti                  Cassini IK
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'pck00010.tpc',
%                                '041014R_SCPSE_01066_04199.bsp',
%                                'cas_v40.tf',
%                                '04161_04164ra.bc',
%                                'cas00071.tsc',
%                                'cas_iss_v10.ti'            )
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function gftfov_ex1()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gftfov_ex1.tm' )
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2004 JUN 11 06:30:00 TDB',                ...
%                               '2004 JUN 11 12:00:00 TDB' } );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         %Initialize inputs for the search.
%         %
%         inst    = 'CASSINI_ISS_NAC';
%         target  = 'PHOEBE';
%         tshape  = 'ELLIPSOID';
%         tframe  = 'IAU_PHOEBE';
%         abcorr  = 'LT+S';
%         obsrvr  = 'CASSINI';
%         step    = 10.;
%         nintvls = MAXWIN;
%
%         result = cspice_gftfov( inst,   target, tshape, tframe, abcorr,  ...
%                                 obsrvr, step,   cnfine, nintvls );
%
%
%         %
%         % List the beginning and ending times in each interval
%         % if result contains data.
%         %
%         for i=1:numel(result)/2
%
%            [left, right] = cspice_wnfetd( result, i );
%
%            output = cspice_timout( [left,right], TIMFMT );
%
%            if( isequal( left, right) )
%
%               disp( ['Event time: ' output(1,:)] )
%
%            else
%
%               disp( ['From : ' output(1,:)] )
%               disp( ['To   : ' output(2,:)] )
%               disp( ' ' )
%
%            end
%
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      From : 2004-JUN-11 07:35:27.066980 (TDB)
%      To   : 2004-JUN-11 08:48:03.954696 (TDB)
%
%      From : 2004-JUN-11 09:02:56.580046 (TDB)
%      To   : 2004-JUN-11 09:35:04.038509 (TDB)
%
%      From : 2004-JUN-11 09:49:56.476397 (TDB)
%      To   : 2004-JUN-11 10:22:04.242879 (TDB)
%
%      From : 2004-JUN-11 10:36:56.283772 (TDB)
%      To   : 2004-JUN-11 11:09:04.397165 (TDB)
%
%      From : 2004-JUN-11 11:23:56.020645 (TDB)
%      To   : 2004-JUN-11 11:56:04.733536 (TDB)
%
%
%-Particulars
%
%   This routine determines a set of one or more time intervals
%   within the confinement window when any portion of a specified
%   target body appears within the field of view of a specified
%   instrument. We'll use the term "visibility event" to designate
%   such an appearance. The set of time intervals resulting from the
%   search is returned as a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient use
%   of this routine in user applications.
%
%   The Search Process
%   ==================
%
%   The search for visibility events is treated as a search for state
%   transitions: times are sought when the state of the target body
%   changes from "not visible" to "visible" or vice versa.
%
%   Step Size
%   =========
%
%   Each interval of the confinement window is searched as follows:
%   first, the input step size is used to determine the time
%   separation at which the visibility state will be sampled.
%   Starting at the left endpoint of an interval, samples will be
%   taken at each step. If a state change is detected, a root has
%   been bracketed; at that point, the "root"--the time at which the
%   state change occurs---is found by a refinement process, for
%   example, via binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the visibility state is constant:
%   the step size should be shorter than the shortest visibility event
%   duration and the shortest period between visibility events, within
%   the confinement window.
%
%   Having some knowledge of the relative geometry of the target and
%   observer can be a valuable aid in picking a reasonable step size.
%   In general, the user can compensate for lack of such knowledge by
%   picking a very short step size; the cost is increased computation
%   time.
%
%   Note that the step size is not related to the precision with which
%   the endpoints of the intervals of the result window are computed.
%   That precision level is controlled by the convergence tolerance.
%
%   Convergence Tolerance
%   =====================
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   "convergence tolerance." The convergence tolerance used by this
%   routine is set by the parameter SPICE_GF_CNVTOL.
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't become the limiting factor in the accuracy of
%   solutions found by this routine. In general the accuracy of input
%   data will be the limiting factor.
%
%   The user may change the convergence tolerance from the default
%   SPICE_GF_CNVTOL value by calling the routine cspice_gfstol, e.g.
%
%      cspice_gfstol( tolerance value in seconds )
%
%   Call cspice_gfstol prior to calling this routine. All subsequent
%   searches will use the updated tolerance value.
%
%   Setting the tolerance tighter than SPICE_GF_CNVTOL is unlikely to be
%   useful, since the results are unlikely to be more accurate.
%   Making the tolerance looser will speed up searches somewhat,
%   since a few convergence steps will be omitted. However, in most
%   cases, the step size is likely to have a much greater affect on
%   processing time than would the convergence tolerance.
%
%   The Confinement Window
%   ======================
%
%   The simplest use of the confinement window is to specify a time
%   interval within which a solution is sought. However, the
%   confinement window can, in some cases, be used to make searches
%   more efficient. Sometimes it's possible to do an efficient search
%   to reduce the size of the time period over which a relatively
%   slow search of interest must be performed.
%
%-Exceptions
%
%   1)  In order for this routine to produce correct results,
%       the step size must be appropriate for the problem at hand.
%       Step sizes that are too large may cause this routine to miss
%       roots; step sizes that are too small may cause this routine
%       to run unacceptably slowly and in some cases, find spurious
%       roots.
%
%       This routine does not diagnose invalid step sizes, except that
%       if the step size is non-positive, an error is signaled by a
%       routine in the call tree of this routine.
%
%   2)  Due to numerical errors, in particular,
%
%          - Truncation error in time values
%          - Finite tolerance value
%          - Errors in computed geometric quantities
%
%       it is *normal* for the condition of interest to not always be
%       satisfied near the endpoints of the intervals comprising the
%       result window.
%
%       The result window may need to be contracted slightly by the
%       caller to achieve desired results. The SPICE window routine
%       cspice_wncond can be used to contract the result window.
%
%   3)  If the name of either the target or observer cannot be
%       translated to a NAIF ID code, an error is signaled by
%       a routine in the call tree of this routine.
%
%   4)  If the specified aberration correction is an unrecognized
%       value, an error is signaled by a routine
%       in the call tree of this routine.
%
%   5)  If the radii of a target body modeled as an ellipsoid cannot
%       be determined by searching the kernel pool for a kernel
%       variable having a name of the form
%
%          'BODYnnn_RADII'
%
%       where nnn represents the NAIF integer code associated with
%       the body, an error is signaled by a routine in the
%       call tree of this routine.
%
%   6)  If the target body coincides with the observer body `obsrvr', an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   7)  If the body model specifier `tshape' is invalid, an error is
%       signaled by either this routine or a routine in the call tree
%       of this routine.
%
%   8)  If a target body-fixed reference frame associated with a
%       non-point target is not recognized, an error is signaled by a
%       routine in the call tree of this routine.
%
%   9)  If a target body-fixed reference frame is not centered at the
%       corresponding target body, an error is signaled by a routine
%       in the call tree of this routine.
%
%   10) If the instrument name `inst' does not have corresponding NAIF
%       ID code, an error is signaled by a routine in the call
%       tree of this routine.
%
%   11) If the FOV parameters of the instrument are not present in
%       the kernel pool, an error is signaled by a routine
%       in the call tree of this routine.
%
%   12) If the FOV boundary has more than SPICE_GF_MAXVRT vertices, an error
%       is signaled by a routine in the call tree of this
%       routine.
%
%   13) If the instrument FOV is polygonal, and this routine cannot
%       find a ray R emanating from the FOV vertex such that maximum
%       angular separation of R and any FOV boundary vector is within
%       the limit (pi/2)-margin radians, an error is signaled
%       by a routine in the call tree of this routine. If the FOV
%       is any other shape, the same error check will be applied with
%       the instrument boresight vector serving the role of R.
%
%   14) If the loaded kernels provide insufficient data to compute a
%       requested state vector, an error is signaled by a
%       routine in the call tree of this routine.
%
%   15) If an error occurs while reading an SPK or other kernel file,
%       the error is signaled by a routine in the call tree
%       of this routine.
%
%   16) If the output SPICE window `result' has insufficient capacity
%       to contain the number of intervals on which the specified
%       visibility condition is met, an error is signaled
%       by a routine in the call tree of this routine.
%
%   17) If any of the input arguments, `inst', `target', `tshape',
%       `tframe', `abcorr', `obsrvr', `step', `cnfine' or `nintvls', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   18) If any of the input arguments, `inst', `target', `tshape',
%       `tframe', `abcorr', `obsrvr', `step', `cnfine' or `nintvls', is
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
%   -  SPK data: ephemeris data for target and observer that
%      describes the ephemeris of these objects for the period
%      defined by the confinement window, `cnfine', must be
%      loaded. If aberration corrections are used, the states of
%      target and observer relative to the solar system barycenter
%      must be calculable from the available ephemeris data.
%      Typically ephemeris data are made available by loading one
%      or more SPK files via cspice_furnsh.
%
%   -  Frame data: if a frame definition is required to convert
%      the observer and target states to the body-fixed frame of
%      the target, that definition must be available in the kernel
%      pool. Typically the definitions of frames not already
%      built-in to SPICE are supplied by loading a frame kernel.
%
%      Data defining the reference frame associated with the
%      instrument designated by `inst' must be available in the
%      kernel pool. Additionally the name `inst' must be associated
%      with an ID code. Normally these data are  made available by
%      loading a frame kernel via cspice_furnsh.
%
%   -  IK data: the kernel pool must contain data such that
%      the Mice routine cspice_getfov may be called to obtain
%      parameters for `inst'. Normally such data are provided by
%      an IK via cspice_furnsh.
%
%   The following data may be required:
%
%   -  PCK data: bodies modeled as triaxial ellipsoids must have
%      orientation data provided by variables in the kernel pool.
%      Typically these data are made available by loading a text
%      PCK file via cspice_furnsh.
%
%      Bodies modeled as triaxial ellipsoids must have semi-axis
%      lengths provided by variables in the kernel pool. Typically
%      these data are made available by loading a text PCK file via
%      cspice_furnsh.
%
%   -  CK data: if the instrument frame is fixed to a spacecraft,
%      at least one CK file will be needed to permit transformation
%      of vectors between that frame and both J2000 and the target
%      body-fixed frame.
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
%   1)  The reference frame associated with `inst' must be
%       centered at the observer or must be inertial. No check is done
%       to ensure this.
%
%   2)  The kernel files to be used by cspice_gftfov must be loaded (normally
%       via the Mice routine cspice_furnsh) before cspice_gftfov is called.
%
%-Required_Reading
%
%   MICE.REQ
%   CK.REQ
%   FRAMES.REQ
%   GF.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%   SPK.REQ
%   TIME.REQ
%   WINDOWS.REQ
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
%   -Mice Version 1.2.0, 11-AUG-2021 (EDW) (JDR)
%
%       Changed input argument name "room" to "nintvls" for consistency
%       with other routines.
%
%       Edited the header to comply with NAIF standard. Updated
%       Example's kernels set to use PDS archived data.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.0, 12-MAY-2012 (EDW)
%
%       Corrected minor typo in header.
%
%       Renamed the argument 'size' to 'room'. "size" is a Matlab function
%       name and it's seriously dumb to use a function name word as an
%       argument name.
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.0, 15-APR-2009 (EDW)
%
%-Index_Entries
%
%   GF target in instrument FOV search
%
%-&

function [result] = cspice_gftfov( inst,   target, tshape, tframe,         ...
                                   abcorr, obsrvr, step,   cnfine, nintvls )

   switch nargin

      case 9

         inst    = zzmice_str(inst);
         target  = zzmice_str(target);
         tshape  = zzmice_str(tshape);
         tframe  = zzmice_str(tframe);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         step    = zzmice_dp(step);
         cnfine  = zzmice_win(cnfine);
         nintvls    = zzmice_int(nintvls,    [1, int32(inf)/2] );

      otherwise

         error ( [ 'Usage: [result] = cspice_gftfov( `inst`, '             ...
                                   '`target`, `tshape`, `tframe` '         ...
                                   '`abcorr`, `obsrvr`, step, '            ...
                                   'cnfine, nintvls )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gftfov_c', inst,  target,  tshape,  tframe,         ...
                                  abcorr, obsrvr, step,                    ...
                                 [zeros(6,1); cnfine], nintvls);

   catch spiceerr
      rethrow(spiceerr)
   end
