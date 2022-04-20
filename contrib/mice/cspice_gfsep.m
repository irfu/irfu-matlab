%-Abstract
%
%   CSPICE_GFSEP determines the time intervals when the angular separation
%   between the position vectors of two target bodies relative to an
%   observer satisfies a numerical relationship.
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
%      targ1    the name of the first body of interest.
%
%               [1,c1] = size(targ1); char = class(targ1)
%
%                  or
%
%               [1,1] = size(targ1); cell = class(targ1)
%
%               You can also supply the integer ID code for the object
%               as an integer string. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               target body.
%
%      shape1   the name of the geometric model used to represent the shape
%               of the `targ1' body.
%
%               [1,c2] = size(shape1); char = class(shape1)
%
%                  or
%
%               [1,1] = size(shape1); cell = class(shape1)
%
%               Models supported by this routine:
%
%                  'SPHERE'        Treat the body as a sphere with radius
%                                  equal to the maximum value of
%                                  BODYnnn_RADII
%
%                  'POINT'         Treat the body as a point;
%                                  radius has value zero.
%
%               The `shape1' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      frame1   the name of the body-fixed reference frame corresponding
%               to `targ1'.
%
%               [1,c3] = size(frame1); char = class(frame1)
%
%                  or
%
%               [1,1] = size(frame1); cell = class(frame1)
%
%               cspice_gfsep does not currently use this argument's value,
%               its use is reserved for future shape models. The value 'NULL'
%               will suffice for 'POINT' and 'SPHERE' shaped bodies.
%
%      targ2    the name of the second body of interest.
%
%               [1,c4] = size(targ2); char = class(targ2)
%
%                  or
%
%               [1,1] = size(targ2); cell = class(targ2)
%
%               You can also supply the integer ID code for the object
%               as an integer string. For example both 'MOON' and '301'
%               are legitimate strings that indicate the Moon is the
%               target body.
%
%      shape2   the name of the geometric model used to represent
%               the shape of the `targ2'.
%
%               [1,c5] = size(shape2); char = class(shape2)
%
%                  or
%
%               [1,1] = size(shape2); cell = class(shape2)
%
%               Models supported by this routine:
%
%                 'SPHERE'        Treat the body as a sphere with radius
%                                 equal to the maximum value of
%                                 BODYnnn_RADII
%
%                 'POINT'         Treat the body as a single point;
%                                 radius has value zero.
%
%               The `shape2' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      frame2   the name of the body-fixed reference frame corresponding
%               to `targ2'.
%
%               [1,c6] = size(frame2); char = class(frame2)
%
%                  or
%
%               [1,1] = size(frame2); cell = class(frame2)
%
%               cspice_gfsep does not currently use this argument's value,
%               its use is reserved for future shape models. The value 'NULL'
%               will suffice for 'POINT' and 'SPHERE' shaped bodies.
%
%      abcorr   describes the aberration corrections to apply to the state
%               evaluations to account for one-way light time and stellar
%               aberration.
%
%               [1,c7] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               This routine accepts the same aberration corrections as does
%               the routine cspice_spkezr. See the header of cspice_spkezr
%               for a detailed description of the aberration correction
%               options. For convenience, the options are listed below:
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       "Reception" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     "Reception" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'CN'       "Reception" case: converged
%                             Newtonian light time correction.
%
%                  'CN+S'     "Reception" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%                  'XLT'      "Transmission" case: correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    "Transmission" case: correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'XCN'      "Transmission" case: converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    "Transmission" case: converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               The `abcorr' string lacks sensitivity to case, and to
%               embedded, leading and trailing blanks.
%
%      obsrvr   the name of the observing body.
%
%               [1,c8] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Optionally, you may supply the ID code of the object as an
%               integer string. For example, both 'EARTH' and '399' are
%               legitimate strings to supply to indicate the
%               observer is Earth.
%
%      relate   the constraint relational operator on the angular separation.
%
%               [1,c9] = size(relate); char = class(relate)
%
%                  or
%
%               [1,1] = size(relate); cell = class(relate)
%
%               The result window found by this routine indicates the time
%               intervals where the constraint is satisfied.
%
%               Supported values of relate and corresponding meanings are
%               shown below:
%
%                  '>'       Separation is greater than the reference
%                            value `refval'.
%
%                  '='       Separation is equal to the reference
%                            value `refval'.
%
%                  '<'       Separation is less than the reference
%                            value `refval'.
%
%                  'ABSMAX'  Separation is at an absolute maximum.
%
%                  'ABSMIN'  Separation is at an absolute  minimum.
%
%                  'LOCMAX'  Separation is at a local maximum.
%
%                  'LOCMIN'  Separation is at a local minimum.
%
%               The caller may indicate that the region of interest
%               is the set of time intervals where the quantity is
%               within a specified angular separation of an absolute extremum.
%               The argument adjust (described below) is used to
%               specify this angular separation.
%
%               Local extrema are considered to exist only in the
%               interiors of the intervals comprising the confinement
%               window:  a local extremum cannot exist at a boundary
%               point of the confinement window.
%
%               Negative Angular Separation
%
%                  For those searches using a SPHERE shape identifier for
%                  either target body, the angular separation function
%                  returns a negative value when the bodies overlap (occult).
%
%               The `relate' string lacks sensitivity to case, leading
%               and trailing blanks.
%
%      refval   reference value used together with `relate' argument to define
%               an equality or inequality to be satisfied by the angular
%               separation between the specified target and observer.
%
%               [1,1] = size(refval); double = class(refval)
%
%               See the discussion of `relate' above for further information.
%
%               The units of 'refval' are radians.
%
%      adjust   value used to modify searches for absolute extrema.
%
%               [1,1] = size(adjust); double = class(adjust)
%
%               When `relate' is set to 'ABSMAX' or 'ABSMIN' and `adjust' is
%               set to a positive value, cspice_gfsep finds times when the
%               angular separation between the bodies is within `adjust'
%               radians of the specified extreme value.
%
%               For `relate' set to 'ABSMAX', the result window contains
%               time intervals when the angular separation has
%               values between ABSMAX - adjust and ABSMAX.
%
%               For `relate' set to 'ABSMIN', the result window contains
%               time intervals when the angular separation has
%               values between ABSMIN and ABSMIN + adjust.
%
%               `adjust' is not used for searches for local extrema,
%               equality or inequality conditions.
%
%      step     time step size to use in the search.
%
%               [1,1] = size(step); double = class(step)
%
%               `step' must be short enough for a search using this step size
%               to locate the time intervals where coordinate function of the
%               observer-target vector is monotone increasing or decreasing.
%               However, `step' must not be *too* short, or the search will
%               take an unreasonable amount of time.
%
%               The choice of `step' affects the completeness but not
%               the precision of solutions found by this routine; the
%               precision is controlled by the convergence tolerance.
%               See the discussion of the parameter SPICE_GF_CNVTOL for
%               details.
%
%               `step' has units of TDB seconds.
%
%      nintvls  value specifying the number of intervals in the internal
%               workspace array used by this routine.
%
%               [1,1] = size(nintvls); int32 = class(nintvls)
%
%               `nintvls' should be at least as large as the number of
%               intervals within the search region on which the specified
%               observer-target vector coordinate function is monotone
%               increasing or decreasing. It does no harm to pick a value of
%               `nintvls' larger than the minimum required to execute the
%               specified search, but if chosen too small, the search will
%               fail.
%
%      cnfine   a SPICE window that confines the time period over which the
%               specified search is conducted.
%
%               [2m,1] = size(cnfine); double = class(cnfine)
%
%               `cnfine' may consist of a single interval or a collection of
%               intervals.
%
%               In some cases the confinement window can be used to
%               greatly reduce the time period that must be searched
%               for the desired solution. See the -Particulars section
%               below for further discussion.
%
%               See the -Examples section below for a code example
%               that shows how to create a confinement window.
%
%               In some cases the observer's state may be computed at
%               times outside of `cnfine' by as much as 2 seconds. See
%               -Particulars for details.
%
%   the call:
%
%      [result] = cspice_gfsep( targ1,  shape1, frame1,  targ2,  shape2,   ...
%                               frame2, abcorr, obsrvr,  relate, refval,   ...
%                               adjust, step,   nintvls, cnfine )
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
%               SPICE_GF_CNVTOL is used to determine when binary
%               searches for roots should terminate: when a root is
%               bracketed within an interval of length SPICE_GF_CNVTOL,
%               the root is considered to have been found.
%
%               The accuracy, as opposed to precision, of roots found
%               by this routine depends on the accuracy of the input
%               data. In most cases, the accuracy of solutions will be
%               inferior to their precision.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Determine the times of local maxima of the angular separation
%      between the moon and sun as observed from earth from
%      Jan 1, 2007 to Jan 1, 2008.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: gfsep_ex1.tm
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
%            de421.bsp                     Planetary ephemeris
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00009.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function gfsep_ex1()
%
%         MAXWIN  =  1000;
%         TIMFMT  = 'YYYY-MON-DD HR:MN:SC.###### (TDB) ::TDB ::RND';
%
%         %
%         % Load kernels.
%         %
%         cspice_furnsh( 'gfsep_ex1.tm' );
%
%         %
%         % Store the time bounds of our search interval in
%         % the cnfine confinement window.
%         %
%         et = cspice_str2et( { '2007 JAN 01', '2008 JAN 01'} );
%
%         cnfine = cspice_wninsd( et(1), et(2) );
%
%         %
%         % Prompt for the inputs.
%         %
%         targ1  = input( 'First body     > ', 's' );
%         targ2  = input( 'Second body    > ', 's' );
%         obsrvr = input( 'Observing body > ', 's' );
%
%         %
%         % Search using a step size of 6 days (in units of seconds).
%         %
%         step   = 6.*cspice_spd;
%         adjust = 0.;
%         refval = 0;
%
%         shape1 = 'SPHERE';
%         frame1 = 'NULL';
%         shape2 = 'SPHERE';
%         frame2 = 'NULL';
%         abcorr = 'NONE';
%         relate = 'LOCMAX';
%         nintvls = MAXWIN;
%
%         result = cspice_gfsep( targ1,   shape1, frame1,                  ...
%                                targ2,   shape2, frame2,                  ...
%                                abcorr,  obsrvr, relate,                  ...
%                                refval,  adjust, step,                    ...
%                                nintvls, cnfine         );
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
%      platform, using 'MOON' as first body, 'EARTH' as second body
%      and 'SUN' as observing body, the output was:
%
%
%      First body     > MOON
%      Second body    > EARTH
%      Observing body > SUN
%      Event time: 2007-JAN-11 11:21:20.214305 (TDB)
%      Event time: 2007-JAN-26 01:43:41.027309 (TDB)
%      Event time: 2007-FEB-10 04:49:53.431964 (TDB)
%      Event time: 2007-FEB-24 13:18:18.953256 (TDB)
%      Event time: 2007-MAR-11 20:41:59.571964 (TDB)
%      Event time: 2007-MAR-26 01:20:26.860201 (TDB)
%      Event time: 2007-APR-10 10:24:39.017514 (TDB)
%      Event time: 2007-APR-24 14:00:49.422728 (TDB)
%      Event time: 2007-MAY-09 21:53:25.643532 (TDB)
%      Event time: 2007-MAY-24 03:14:05.873982 (TDB)
%      Event time: 2007-JUN-08 07:24:13.686616 (TDB)
%      Event time: 2007-JUN-22 16:45:56.506850 (TDB)
%      Event time: 2007-JUL-07 15:30:03.706532 (TDB)
%      Event time: 2007-JUL-22 06:26:17.397353 (TDB)
%      Event time: 2007-AUG-05 23:03:21.625229 (TDB)
%      Event time: 2007-AUG-20 20:14:56.801678 (TDB)
%      Event time: 2007-SEP-04 07:13:25.162360 (TDB)
%      Event time: 2007-SEP-19 10:16:42.721117 (TDB)
%      Event time: 2007-OCT-03 17:11:17.188939 (TDB)
%      Event time: 2007-OCT-19 00:30:31.300060 (TDB)
%      Event time: 2007-NOV-02 05:43:48.902220 (TDB)
%      Event time: 2007-NOV-17 14:38:21.314771 (TDB)
%      Event time: 2007-DEC-01 20:50:27.562519 (TDB)
%      Event time: 2007-DEC-17 04:04:46.933247 (TDB)
%      Event time: 2007-DEC-31 13:43:52.558812 (TDB)
%
%
%   2) Determine the time of local maxima elongation of the
%      Moon as seen from Earth for the same time interval
%      as the previous example, i.e. find the local maxima of
%      the angular separation between the Moon and the Sun as
%      seen from the Earth, by running the code in example #1.
%
%
%      When Example #1 was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, using 'MOON' as first body, 'SUN' as second body
%      and 'EARTH' as observing body, the output was:
%
%
%      First body     > MOON
%      Second body    > SUN
%      Observing body > EARTH
%      Event time: 2007-JAN-03 14:20:24.617627 (TDB)
%      Event time: 2007-FEB-02 06:16:24.101517 (TDB)
%      Event time: 2007-MAR-03 23:22:41.994972 (TDB)
%      Event time: 2007-APR-02 16:49:16.135505 (TDB)
%      Event time: 2007-MAY-02 09:41:43.830081 (TDB)
%      Event time: 2007-JUN-01 01:03:44.527470 (TDB)
%      Event time: 2007-JUN-30 14:15:26.576292 (TDB)
%      Event time: 2007-JUL-30 01:14:49.000963 (TDB)
%      Event time: 2007-AUG-28 10:39:01.388249 (TDB)
%      Event time: 2007-SEP-26 19:25:51.509426 (TDB)
%      Event time: 2007-OCT-26 04:30:56.625105 (TDB)
%      Event time: 2007-NOV-24 14:31:04.331185 (TDB)
%      Event time: 2007-DEC-24 01:40:12.235392 (TDB)
%
%
%-Particulars
%
%   This routine determines a set of one or more time intervals
%   within the confinement window for which the angular separation
%   between the two bodies satisfies some defined relationship.
%   The resulting set of intervals is returned as a SPICE window.
%
%   Below we discuss in greater detail aspects of this routine's
%   solution process that are relevant to correct and efficient
%   use of this routine in user applications.
%
%   The Search Process
%   ==================
%
%   Regardless of the type of constraint selected by the caller, this
%   routine starts the search for solutions by determining the time
%   periods, within the confinement window, over which the specified
%   angular separation function is monotone increasing and monotone
%   decreasing. Each of these time periods is represented by a SPICE
%   window. Having found these windows, all of the angular separation
%   function's local extrema within the confinement window are known.
%   Absolute extrema then can be found very easily.
%
%   Within any interval of these "monotone" windows, there will be at
%   most one solution of any equality constraint. Since the boundary
%   of the solution set for any inequality constraint is contained in
%   the union of
%
%   -  the set of points where an equality constraint is met
%
%   -  the boundary points of the confinement window
%
%   the solutions of both equality and inequality constraints can be
%   found easily once the monotone windows have been found.
%
%
%   Step Size
%   =========
%
%   The monotone windows (described above) are found using a two-step
%   search process. Each interval of the confinement window is
%   searched as follows: first, the input step size is used to
%   determine the time separation at which the sign of the rate of
%   change of angular separation (angular separation rate) will be
%   sampled. Starting at the left endpoint of an interval, samples
%   will be taken at each step. If a change of sign is found, a
%   root has been bracketed; at that point, the time at which the
%   angular separation rate is zero can be found by a refinement
%   process, for example, using a binary search.
%
%   Note that the optimal choice of step size depends on the lengths
%   of the intervals over which the distance function is monotone:
%   the step size should be shorter than the shortest of these
%   intervals (within the confinement window).
%
%   The optimal step size is *not* necessarily related to the lengths
%   of the intervals comprising the result window. For example, if
%   the shortest monotone interval has length 10 days, and if the
%   shortest result window interval has length 5 minutes, a step size
%   of 9.9 days is still adequate to find all of the intervals in the
%   result window. In situations like this, the technique of using
%   monotone windows yields a dramatic efficiency improvement over a
%   state-based search that simply tests at each step whether the
%   specified constraint is satisfied. The latter type of search can
%   miss solution intervals if the step size is longer than the
%   shortest solution interval.
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
%
%   Convergence Tolerance
%   =====================
%
%   As described above, the root-finding process used by this routine
%   involves first bracketing roots and then using a search process
%   to locate them. "Roots" are both times when local extrema are
%   attained and times when the distance function is equal to a
%   reference value. All endpoints of the intervals comprising the
%   result window are either endpoints of intervals of the
%   confinement window or roots.
%
%   Once a root has been bracketed, a refinement process is used to
%   narrow down the time interval within which the root must lie.
%   This refinement process terminates when the location of the root
%   has been determined to within an error margin called the
%   "convergence tolerance." The default convergence tolerance
%   used by this routine is set by the parameter SPICE_GF_CNVTOL (defined
%   in MiceGF.m).
%
%   The value of SPICE_GF_CNVTOL is set to a "tight" value so that the
%   tolerance doesn't become the limiting factor in the accuracy of
%   solutions found by this routine. In general the accuracy of input
%   data will be the limiting factor.
%
%   The user may change the convergence tolerance from the default
%   SPICE_GF_CNVTOL value by calling the routine cspice_gfstol, e.g.
%
%      cspice_gfstol( tolerance value );
%
%   Call cspice_gfstol prior to calling this routine. All subsequent
%   searches will use the updated tolerance value.
%
%   Setting the tolerance tighter than SPICE_GF_CNVTOL is unlikely to be
%   useful, since the results are unlikely to be more accurate.
%   Making the tolerance looser will speed up searches somewhat,
%   since a few convergence steps will be omitted. However, in most
%   cases, the step size is likely to have a much greater effect
%   on processing time than would the convergence tolerance.
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
%   Certain types of searches require the state of the observer,
%   relative to the solar system barycenter, to be computed at times
%   slightly outside the confinement window `cnfine'. The time window
%   that is actually used is the result of "expanding" `cnfine' by a
%   specified amount "T": each time interval of `cnfine' is expanded by
%   shifting the interval's left endpoint to the left and the right
%   endpoint to the right by T seconds. Any overlapping intervals are
%   merged. (The input argument `cnfine' is not modified.)
%
%   The window expansions listed below are additive: if both
%   conditions apply, the window expansion amount is the sum of the
%   individual amounts.
%
%   -  If a search uses an equality constraint, the time window
%      over which the state of the observer is computed is expanded
%      by 1 second at both ends of all of the time intervals
%      comprising the window over which the search is conducted.
%
%   -  If a search uses stellar aberration corrections, the time
%      window over which the state of the observer is computed is
%      expanded as described above.
%
%   When light time corrections are used, expansion of the search
%   window also affects the set of times at which the light time-
%   corrected state of the target is computed.
%
%   In addition to the possible 2 second expansion of the search
%   window that occurs when both an equality constraint and stellar
%   aberration corrections are used, round-off error should be taken
%   into account when the need for data availability is analyzed.
%
%   Negative Angular Separation
%   ===========================
%
%   For those searches using a SPHERE shape identifier for both
%   target bodies, the angular separation function returns a
%   negative value when the bodies overlap (occult), e.g.
%   a search for an ABSMIN of angular separation in a
%   confinement window covering an occultation event will
%   return the time when the apparent center of the
%   occulting body passes closest to the apparent center of
%   the occulted body.
%
%
%   Elongation
%   ===========================
%
%   The angular separation of two targets as seen from an observer
%   where one of those targets is the sun is known as elongation.
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
%       This routine does not diagnose invalid step sizes, except
%       that if the step size is non-positive, an error is signaled
%       by a routine in the call tree of this routine.
%
%   2)  Due to numerical errors, in particular,
%
%          - truncation error in time values
%          - finite tolerance value
%          - errors in computed geometric quantities
%
%       it is *normal* for the condition of interest to not always be
%       satisfied near the endpoints of the intervals comprising the
%       `result' window. One technique to handle such a situation,
%       slightly contract `result' using the window routine cspice_wncond.
%
%   3)  If `result' has insufficient capacity to contain the
%       number of intervals on which the specified distance condition
%       is met, an error is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If an error (typically cell overflow) occurs during
%       window arithmetic, the error is signaled by a routine
%       in the call tree of this routine.
%
%   5)  If the relational operator `relate' is not recognized, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   6)  If the aberration correction specifier contains an
%       unrecognized value, an error is signaled by a routine in the
%       call tree of this routine.
%
%   7)  If `adjust' is negative, an error is signaled by a routine in
%       the call tree of this routine.
%
%   8)  If either of the input body names, `targ1', `targ2' do not map
%       to NAIF ID codes, an error is signaled by a routine in the
%       call tree of this routine.
%
%   9)  If either of the input body shape names, `shape1', `shape2',
%       are not recognized by the GF subsystem, an error is signaled
%       by a routine in the call tree of this routine.
%
%   10) If either of the input body frame names, `frame1', `frame2',
%       are not recognized by the frame subsystem, an error is
%       signaled by a routine in the call tree of this routine.
%
%   11) If either of the input body frames, `frame1', `frame2',
%       are not centered on the corresponding body (`frame1' on `targ1',
%       `frame2' on `targ2'), an error is signaled by a routine in the
%       call tree of this routine.
%
%   12) If required ephemerides or other kernel data are not
%       available, an error is signaled by a routine in the call tree
%       of this routine.
%
%   13) If any of the input arguments, `targ1', `shape1', `frame1',
%       `targ2', `shape2', `frame2', `abcorr', `obsrvr', `relate',
%       `refval', `adjust', `step', `nintvls' or `cnfine', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   14) If any of the input arguments, `targ1', `shape1', `frame1',
%       `targ2', `shape2', `frame2', `abcorr', `obsrvr', `relate',
%       `refval', `adjust', `step', `nintvls' or `cnfine', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
%
%-Files
%
%   Appropriate SPK and PCK kernels must be loaded by the
%   calling program before this routine is called.
%
%   The following data are required:
%
%   -  SPK data: the calling application must load ephemeris data
%      for the targets, observer, and any intermediate objects in
%      a chain connecting the targets and observer that cover the
%      time period specified by the window `cnfine'. If aberration
%      corrections are used, the states of target and observer
%      relative to the solar system barycenter must be calculable
%      from the available ephemeris data. Typically ephemeris data
%      are made available by loading one or more SPK files using
%      cspice_furnsh.
%
%   -  PCK data: bodies modeled as triaxial ellipsoids must have
%      semi-axis lengths provided by variables in the kernel pool.
%      Typically these data are made available by loading a text
%      PCK file using cspice_furnsh.
%
%   -  If non-inertial reference frames are used, then PCK
%      files, frame kernels, C-kernels, and SCLK kernels may be
%      needed.
%
%   -  In some cases the observer's state may be computed at times
%      outside of `cnfine' by as much as 2 seconds; data required to
%      compute this state must be provided by loaded kernels. See
%      -Particulars for details.
%
%   Such kernel data are normally loaded once per program
%   run, NOT every time this routine is called.
%
%-Restrictions
%
%   1)  The kernel files to be used by this routine must be loaded
%       (normally using the Mice routine cspice_furnsh) before this
%       routine is called.
%
%   2)  This routine has the side effect of re-initializing the
%       angular separation quantity utility package. Callers may
%       need to re-initialize the package after calling this routine.
%
%   3)  Due to the current logic implemented in SPICE, a direct
%       search for zero angular separation of two point targets will
%       always fails, i.e.,
%
%          relate = '='
%          refval = 0.0
%
%       Use `relate' values of 'ABSMIN' or 'LOCMIN' to detect such an
%       event(s).
%
%-Required_Reading
%
%   MICE.REQ
%   GF.REQ
%   SPK.REQ
%   CK.REQ
%   TIME.REQ
%   WINDOWS.REQ
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
%   -Mice Version 1.1.0, 03-NOV-2021 (EDW) (JDR)
%
%       Updated header to describe use of expanded confinement window.
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel, modified example code to prompt for
%       the required inputs and added a second example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 17-MAR-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Typo correction in version IDs in -Version section.
%
%   -Mice Version 1.0.2, 05-SEP-2012 (EDW)
%
%       Edit to comments to correct search description.
%
%       Header updated to describe use of cspice_gfstol.
%
%   -Mice Version 1.0.1, 29-DEC-2009 (EDW)
%
%       Edited argument descriptions. Removed mention of "ELLIPSOID"
%       shape from 'shape1' and 'shape2' as that option is not yet
%       implemented.
%
%   -Mice Version 1.0.0, 15-APR-2009 (NJB) (EDW)
%
%-Index_Entries
%
%   GF angular separation search
%
%-&

function [result] = cspice_gfsep( targ1, shape1, frame1,          ...
                                  targ2, shape2, frame2,          ...
                                  abcorr, obsrvr, relate, refval, ...
                                  adjust, step, nintvls, cnfine )

   switch nargin

      case 14

         targ1   = zzmice_str(targ1);
         shape1  = zzmice_str(shape1);
         frame1  = zzmice_str(frame1);
         targ2   = zzmice_str(targ2);
         shape2  = zzmice_str(shape2);
         frame2  = zzmice_str(frame2);
         abcorr  = zzmice_str(abcorr);
         obsrvr  = zzmice_str(obsrvr);
         relate  = zzmice_str(relate);
         refval  = zzmice_dp(refval);
         adjust  = zzmice_dp(adjust);
         step    = zzmice_dp(step);
         nintvls = zzmice_int(nintvls, [1, int32(inf)/2] );
         cnfine  = zzmice_win(cnfine);

      otherwise

         error ( [ 'Usage: [result] = cspice_gfsep( `targ1`, `shape1`, ' ...
                               '`frame1`, `targ2`, `shape2`, `frame2`, ' ...
                               '`abcorr`, `obsrvr`, `relate`, refval, '  ...
                               'adjust, step, nintvls, cnfine )' ] )

   end

   %
   % Call the GF routine, add to 'cnfine' the space needed for
   % the control segment.
   %
   try

      [result] = mice('gfsep_c',  targ1, shape1, frame1,          ...
                                  targ2, shape2, frame2,          ...
                                  abcorr, obsrvr, relate, refval, ...
                                  adjust, step, nintvls,          ...
                                  [zeros(6,1); cnfine] );
   catch spiceerr
      rethrow(spiceerr)
   end
