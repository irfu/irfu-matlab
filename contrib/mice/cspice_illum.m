%-Abstract
%
%   Deprecated: This routine has been superseded by the Mice routine
%   cspice_ilumin. This routine is supported for purposes of
%   backward compatibility only.
%
%   CSPICE_ILLUM calculates the illumination angles at a specified
%   surface point of a target body.
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
%      target   the name of the target body.
%
%               [1,c1] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               `target' is case-insensitive, and leading and trailing
%               blanks in `target' are not significant. Optionally, you may
%               supply a string containing the integer ID code for the
%               object. For example both 'MOON' and '301' are legitimate
%               strings that indicate the moon is the target body.
%
%      et       the epoch(s), specified in ephemeris seconds past J2000, at
%               which the apparent illumination angles at the specified
%               surface point on the target body, as seen from the observing
%               body, are to be computed.
%
%               [1,n] = size(et); double = class(et)
%
%      abcorr   the aberration correction to be used in computing the
%               location and orientation of the target body and the location
%               of the Sun.
%
%               [1,c2] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               Possible values are:
%
%                  'NONE'        No aberration correction.
%
%                  'LT'          Correct the position and
%                                orientation of target body for
%                                light time, and correct the
%                                position of the Sun for light
%                                time.
%
%                  'LT+S'        Correct the observer-target vector
%                                for light time and stellar
%                                aberration, correct the
%                                orientation of the target body
%                                for light time, and correct the
%                                target-Sun vector for light time
%                                and stellar aberration.
%
%                  'CN'          Converged Newtonian light time
%                                correction. In solving the light
%                                time equation, the 'CN'
%                                correction iterates until the
%                                solution converges (three
%                                iterations on all supported
%                                platforms). Whether the 'CN+S'
%                                solution is substantially more
%                                accurate than the 'LT' solution
%                                depends on the geometry of the
%                                participating objects and on the
%                                accuracy of the input data. In
%                                all cases this routine will
%                                execute more slowly when a
%                                converged solution is computed.
%                                See the -Particulars section of
%                                cspice_spkezr for a discussion of
%                                precision of light time
%                                corrections.
%
%                                Both the state and rotation of
%                                the target body are corrected for
%                                light time.
%
%                     'CN+S'     Converged Newtonian light time
%                                correction and stellar aberration
%                                correction.
%
%                                Both the state and rotation of
%                                the target body are corrected for
%                                light time.
%
%      obsrvr   the name of the observing body, typically a spacecraft, the
%               earth, or a surface point on the earth.
%
%               [1,c3] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               `obsrvr' is case-insensitive, and leading and trailing
%               blanks in `obsrvr' are not significant. Optionally, you may
%               supply a string containing the integer ID code for the
%               object. For example both 'EARTH' and '399' are legitimate
%               strings that indicate the earth is the observer.
%
%               `obsrvr' may be not be identical to `target'.
%
%      spoint   an array representing a surface point or points on the target
%               body, expressed in rectangular body-fixed (body equator and
%               prime meridian) coordinates.
%
%               [3,n] = size(spoint); double = class(spoint)
%
%               Each `spoint' element (spoint(:,i)) corresponds to the same
%               element index in `et' (et(i)) and need not be visible from the
%               observer's location at time `et'.
%
%               Note: The design of cspice_illum supposes the input 'spoint'
%               originates as the output of another Mice routine. Still, in
%               the event the user requires an 'spoint' constant over a vector
%               of 'et', such as a constant station location at (x,y,z),
%               construct 'spoint' with the MATLAB code:
%
%                  N            = numel(et);
%                  spoint       = eye(3, N);
%                  spoint(1,:)  = x;
%                  spoint(2,:)  = y;
%                  spoint(3,:)  = z;
%
%   the call:
%
%      [phase, solar, emissn] = cspice_illum( target, et,    abcorr,       ...
%                                             obsrvr, spoint        )
%
%   returns:
%
%      phase    the phase angle(s) at `spoint', as seen from `obsrvr' at time
%               `et'.
%
%               [1,n] = size(phase); double = class(phase)
%
%               This is the angle between the spoint-obsrvr vector and the
%               spoint-Sun vector. Units are radians. The range of `phase' is
%               [0, pi]. See -Particulars below for a detailed discussion of
%               the definition.
%
%      solar    the solar incidence angle(s) at `spoint', as seen from
%               `obsrvr' at time `et'.
%
%               [1,n] = size(solar); double = class(solar)
%
%               This is the angle between the surface normal vector at
%               `spoint' and the spoint-Sun vector. Units are radians. The
%               range of `solar' is [0, pi]. See -Particulars below for a
%               detailed discussion of the definition.
%
%      emissn   the emission angle(s) at `spoint', as seen from `obsrvr' at
%               time `et'.
%
%               [1,n] = size(emissn); double = class(emissn)
%
%               This is the angle between the surface normal vector at
%               `spoint' and the spoint-observer vector. Units are radians.
%               The range of `emissn' is [0, pi]. See -Particulars below for
%               a detailed discussion of the definition.
%
%               `phase', `solar' and `emissn' return with the same
%               vectorization measure, N, as `et' and `spoint'.
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
%   1) Compute the time evolution of the phase, solar, and
%      emission angles for the intercept sub-point of the
%      MGS orbiter from Aug 1, 2003 to Aug 3, 2003.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: illum_ex1.tm
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
%            mar097.bsp                       Mars satellite ephemeris
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0011.tls                     Leapseconds
%            mgs_ext12_ipng_mgs95j.bsp        MGS ephemeris
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'mar097.bsp',
%                                'pck00010.tpc',
%                                'naif0011.tls',
%                                'mgs_ext12_ipng_mgs95j.bsp' )
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function illum_ex1()
%
%         TARGET   = 'MARS';
%         OBSERVER = 'MGS';
%         CORRECT  = 'LT+S';
%
%         %
%         % Define the start and stop time for the computations.
%         %
%         START_TIME = '1 Aug 2003';
%         STOP_TIME  = '3 Aug 2003';
%
%         %
%         % Number of steps?
%         %
%         STEP = 6;
%
%         %
%         % Load the standard leapseconds and PCK kernels,
%         % and the Mars and MGS SPK kernels.
%         %
%         cspice_furnsh( 'illum_ex1.tm' )
%
%         %
%         % Convert the strings to ephemeris time J2000.
%         %
%         et_start = cspice_str2et( START_TIME );
%         et_stop = cspice_str2et( STOP_TIME );
%
%         %
%         % Length of a step in seconds for STEP steps.
%         %
%         space = (et_stop - et_start)/STEP;
%
%         %
%         % Create a vector of ephemeris times.
%         %
%         et = [0:(STEP-1)]*space + et_start;
%
%         %
%         % Start at 'et_start', take STEP steps
%         % of space 'length'. At each time, calculate the
%         % intercept sub-point of the observer, then calculate
%         % the illumination angles at the sub-point.
%         %
%         [pos, alt] = cspice_subpt( 'Intercept', TARGET, et, ...
%                                     CORRECT,    OBSERVER        );
%
%         [ phase, solar, emissn] = cspice_illum( TARGET, et, ...
%                                        CORRECT, OBSERVER, pos   );
%
%         %
%         % Convert the et value to UTC for human comprehension.
%         %
%         utc    = cspice_et2utc( et, 'C', 3 );
%         phase  = phase  * cspice_dpr;
%         solar  = solar  * cspice_dpr;
%         emissn = emissn * cspice_dpr;
%
%         for i = 1:STEP
%
%            %
%            % Output the times and lighting angles in degrees.
%            %
%            txt = sprintf( 'UTC           : %s', utc(i,:) );
%            disp( txt )
%
%            txt = sprintf( 'Emission angle: %14.6f', emissn(i) );
%            disp( txt )
%
%            txt = sprintf( 'Solar angle   : %14.6f', solar(i)  );
%            disp( txt )
%
%            txt = sprintf( 'Phase angle   : %14.6f', phase(i)  );
%            disp( txt )
%
%            disp( ' ' )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      UTC           : 2003 AUG 01 00:00:00.000
%      Emission angle:       0.195943
%      Solar angle   :     141.888500
%      Phase angle   :     141.855904
%
%      UTC           : 2003 AUG 01 08:00:00.000
%      Emission angle:       0.124756
%      Solar angle   :     131.066188
%      Phase angle   :     131.143966
%
%      UTC           : 2003 AUG 01 16:00:00.000
%      Emission angle:       0.330871
%      Solar angle   :     111.685115
%      Phase angle   :     111.947851
%
%      UTC           : 2003 AUG 02 00:00:00.000
%      Emission angle:       0.238157
%      Solar angle   :      89.263211
%      Phase angle   :      89.469408
%
%      UTC           : 2003 AUG 02 08:00:00.000
%      Emission angle:       0.082844
%      Solar angle   :      67.144481
%      Phase angle   :      67.111755
%
%      UTC           : 2003 AUG 02 16:00:00.000
%      Emission angle:       0.316717
%      Solar angle   :      48.044076
%      Phase angle   :      47.887931
%
%
%-Particulars
%
%   The term "illumination angles" refers to following set of
%   angles:
%
%
%      solar incidence angle    Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               `sun'.
%
%      emission angle           Angle between the surface normal at
%                               the specified surface point and the
%                               vector from the surface point to the
%                               observer.
%
%      phase angle              Angle between the vectors from the
%                               surface point to the observing body's
%                               location and from the surface point
%                               to the `sun'.
%
%
%   The diagram below illustrates the geometrical relationships
%   defining these angles. The labels for the solar incidence,
%   emission, and phase angles are "s.i.", "e.", and "phase".
%
%
%                                                    *
%                                                   sun
%
%                  surface normal vector
%                            ._                 _.
%                            |\                 /|  sun vector
%                              \    phase      /
%                               \   .    .    /
%                               .            .
%                                 \   ___   /
%                            .     \/     \/
%                                  _\ s.i./
%                           .    /   \   /
%                           .   |  e. \ /
%       *             <--------------- *  surface point on
%    viewing            vector            target body
%    location           to viewing
%    (observer)         location
%
%
%   Note that if the target-observer vector, the target normal vector
%   at the surface point, and the target-sun vector are coplanar,
%   then phase is the sum of incidence and emission. This is rarely
%   true; usually
%
%      phase angle  <  solar incidence angle + emission angle
%
%   All of the above angles can be computed using light time
%   corrections, light time and stellar aberration corrections, or
%   no aberration corrections. The way aberration corrections
%   are used is described below.
%
%   Care must be used in computing light time corrections. The
%   guiding principle used here is "describe what appears in
%   an image." We ignore differential light time; the light times
%   from all points on the target to the observer are presumed to be
%   equal.
%
%
%      Observer-target body vector
%      ---------------------------
%
%      Let `et' be the epoch at which an observation or remote
%      sensing measurement is made, and let et - lt ("LT" stands
%      for "light time") be the epoch at which the photons received
%      at `et' were emitted from the body (we use the term "emitted"
%      loosely here).
%
%      The correct observer-target vector points from the observer's
%      location at `et' to the target body's location at et - lt.
%      The target-observer vector points in the opposite direction.
%
%      Since light time corrections are not symmetric, the correct
%      target-observer vector CANNOT be found by computing the light
%      time corrected position of the observer as seen from the
%      target body.
%
%
%      Target body's orientation
%      -------------------------
%
%      Using the definitions of `et' and `lt' above, the target
%      body's orientation at et - lt is used. The surface
%      normal is dependent on the target body's orientation, so
%      the body's orientation model must be evaluated for the correct
%      epoch.
%
%
%      Target body -- `sun' vector
%      ---------------------------
%
%      All surface features on the target body will appear in
%      a measurement made at `et' as they were at et-lt. In
%      particular, lighting on the target body is dependent on
%      the apparent location of the `sun' as seen from the target
%      body at et-lt. So, a second light time correction is used
%      in finding the apparent location of the `sun'.
%
%
%   Stellar aberration corrections, when used, are applied as follows:
%
%
%      Observer-target body vector
%      ---------------------------
%
%      In addition to light time correction, stellar aberration is
%      used in computing the apparent target body position as seen
%      from the observer's location at time `et'. This apparent
%      position defines the observer-target body vector.
%
%
%      Target body-Sun vector
%      ----------------------
%
%      The target body-Sun vector is the apparent position of the `sun',
%      corrected for light time and stellar aberration, as seen from
%      the target body at time et-lt. Note that the target body's
%      position is not affected by the stellar aberration correction
%      applied in finding its apparent position as seen by the
%      observer.
%
%
%   Once all of the vectors, as well as the target body's
%   orientation, have been computed with the proper aberration
%   corrections, the element of time is eliminated from the
%   computation. The problem becomes a purely geometrical one,
%   and is described by the diagram above.
%
%-Exceptions
%
%   1)  If `target' and `obsrvr' are not distinct, the error
%       SPICE(BODIESNOTDISTINCT) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If no SPK (ephemeris) data are available for the observer,
%       target, and Sun at the time specified by `et', an error is
%       signaled by a routine in the call tree of this routine. If
%       light time corrections are used, SPK data for the target body
%       must be available at the time et - lt, where `lt' is the one-way
%       light time from the target to the observer at `et'.
%       Additionally, SPK data must be available for the Sun at the
%       time et - lt - lt2, where `lt2' is the light time from the Sun
%       to the target body at time et - lt.
%
%   3)  If PCK data defining the orientation or shape of the target
%       body are unavailable, an error is signaled by a routine in the
%       call tree of this routine.
%
%   4)  If no body-fixed frame is associated with the target body, the
%       error SPICE(NOFRAME) is signaled by a routine in the call tree
%       of this routine.
%
%   5)  If name of target or observer cannot be translated to its NAIF
%       ID code, the error SPICE(IDCODENOTFOUND) is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If radii for `target' are not found in the kernel pool, an error
%       is signaled by a routine in the call tree of this routine.
%
%   7)  If the size of the `target' body radii kernel variable is not
%       three, an error is signaled by a routine in the call tree of
%       this routine.
%
%   8)  If any of the three `target' body radii is less-than or equal to
%       zero, an error is signaled by a routine in the call tree of
%       this routine.
%
%   9)  If any of the input arguments, `target', `et', `abcorr',
%       `obsrvr' or `spoint', is undefined, an error is signaled by
%       the Matlab error handling system.
%
%   10) If any of the input arguments, `target', `et', `abcorr',
%       `obsrvr' or `spoint', is not of the expected type, or it does
%       not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%   11) If the input vectorizable arguments `et' and `spoint' do not
%       have the same measure of vectorization (N), an error is
%       signaled by the Mice interface.
%
%-Files
%
%   No files are input to this routine. However, cspice_illum expects
%   that the appropriate SPK and PCK files have been loaded via
%   cspice_furnsh.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
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
%   B.V. Semenov        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel. Reduced the time window and the number of
%       steps used in the code example.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.3, 18-NOV-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.0.2, 18-MAY-2010 (BVS)
%
%       Index lines now state that this routine is deprecated.
%
%   -Mice Version 1.0.1, 30-DEC-2008 (EDW)
%
%       Edits to header; -Abstract now states that this routine is
%       deprecated.
%
%       Corrected misspellings.
%
%   -Mice Version 1.0.0, 15-DEC-2005 (EDW)
%
%-Index_Entries
%
%   DEPRECATED illumination angles
%   DEPRECATED lighting angles
%   DEPRECATED phase angle
%   DEPRECATED emission angle
%   DEPRECATED solar incidence angle
%
%-&

function [phase, solar, emissn] = cspice_illum( target, et, abcorr, ...
                                                obsrvr, spoint )

   switch nargin
      case 5

         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
         spoint = zzmice_dp(spoint);

      otherwise

         error ( ['Usage: [_phase_, _solar_, _emissn_] = '           ...
                          'cspice_illum( `target`, _et_, `abcorr`, ' ...
                          '`obsrvr`, _spoint(3)_)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [phase, solar, emissn] = mice('illum_c', target, et, ...
                                               abcorr, obsrvr, spoint );
   catch spiceerr
      rethrow(spiceerr)
   end




