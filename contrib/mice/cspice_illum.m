%-Abstract
%
%   CSPICE_ILLUM calculates the illumination angles at a specified
%   surface point of a target body.
%
%   Deprecated: This routine has been superseded by the routine
%   cspice_ilumin. This routine is supported for purposes of
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
%      target   the name of the target body. Optionally, you may supply
%               a string containing the integer ID code for the object.
%               For example both 'MOON' and '301' are legitimate strings
%               that indicate the moon is the target body.
%
%               [1,c1] = size(target); char = class(target)
%
%                  or
%
%               [1,1] = size(target); cell = class(target)
%
%               The target and observer define a state vector whose
%               position component points from the observer to the target.
%
%               Case and leading or trailing blanks are not significant
%               in the string 'target'.
%
%      et       the  ephemeris time(s) at which to compute the
%               apparent illumination angles at the specified surface
%               point on the target body, as seen from the observing
%               body.
%
%               [1,n] = size(et); double = class(et)
%
%      abcorr   describes the aberration corrections to apply to the state
%               evaluations to account for one-way light time and stellar
%               aberration.
%
%               [1,c2] = size(abcorr); char = class(abcorr)
%
%                  or
%
%               [1,1] = size(abcorr); cell = class(abcorr)
%
%               This routine accepts the same aberration corrections as does
%               the routine spkezr_c. See the header of spkezr_c for a
%               detailed description of the aberration correction options.
%               For convenience, the options are listed below:
%
%                  'NONE'     Apply no correction.
%
%                  'LT'       "Reception" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'LT+S'     "Reception" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'CN'       "Reception" case:  converged
%                             Newtonian light time correction.
%
%                  'CN+S'     "Reception" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%                  'XLT'      "Transmission" case:  correct for
%                             one-way light time using a Newtonian
%                             formulation.
%
%                  'XLT+S'    "Transmission" case:  correct for
%                             one-way light time and stellar
%                             aberration using a Newtonian
%                             formulation.
%
%                  'XCN'      "Transmission" case:  converged
%                             Newtonian light time correction.
%
%                  'XCN+S'    "Transmission" case:  converged
%                             Newtonian light time and stellar
%                             aberration corrections.
%
%               The 'abcorr' string lacks sensitivity to case, and to embedded,
%               leading and trailing blanks.
%
%      obsrvr   the name of the observing body, typically a spacecraft,
%               the earth, or a surface point on the earth. Optionally, 
%               you may supply a string containing the integer ID code for
%               the object.  For example both "EARTH" and "399" are legitimate
%               strings that indicate the earth is the observer.
%
%               [1,c3] = size(obsrvr); char = class(obsrvr)
%
%                  or
%
%               [1,1] = size(obsrvr); cell = class(obsrvr)
%
%               Case and leading or trailing blanks are not significant
%               in the string 'obsrvr'.
%
%      spoint   an array representing a surface point or points on the target
%               body, expressed in rectangular body-fixed (body equator and
%               prime meridian) coordinates. Each 'spoint' element
%               (spoint(:,i)) corresponds to the same element index in 'et'
%               (et(i)) and need not be visible from the observer's location
%               at time 'et'.
%
%               [3,n] = size(spoint); double = class(spoint)
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
%      [phase, solar, emissn] = cspice_illum( target, et, abcorr, ...
%                                             obsrvr, spoint)
%
%   returns:
%
%      phase    the values(s) of phase angle at 'spoint', as seen from
%               'obsrvr' at time 'et'. This is the angle between the
%               'spoint'-'obsrvr' vector and the 'spoint'-sun vector.
%               Units are radians. The range of 'phase' is [0, pi].
%
%               [1,n] = size(phase); double = class(phase)
%
%      solar    the values(s) of incidence angle at `spoint', as seen
%               from 'obsrvr' at time 'et'. This is the angle between
%               the surface normal vector at 'spoint' and the
%               'spoint'-sun vector. Units are radians. The range of
%               'solar' is [0, pi].
%
%               [1,n] = size(solar); double = class(solar)
%
%      emissn   the values(s) of the emission angle at 'spoint', as seen
%               from  'obsrvr' at time 'et'. This is the angle between
%               the surface normal vector at 'spoint' and the
%               'spoint'-observer vector.  Units are radians.  The range
%               of 'emissn' is [0, pi].
%
%               [1,n] = size(emissn); double = class(emissn)
%
%               'phase', 'solar', 'emissn' return with the same
%               vectorization measure, N, as 'et'.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      % Compute the time evolution of the phase, solar, and
%      % emission angles for the intercept sub-point of the
%      % MGS orbiter from Feb 1, 2003 to April 1, 2003.
%      %
%      TARGET   = 'MARS';
%      OBSERVER = 'MGS';
%      CORRECT  = 'LT+S';
%
%      %
%      % Assign the MGS SPK kernel path-name to a string variable.
%      %
%      MGS = '/kernels/MGS/spk/spk_m_030102-030403_021004.bsp';
%
%      %
%      % Define the start and stop time for the computations.
%      %
%      START_TIME = '1 Feb 2003';
%      STOP_TIME  = '1 APR 2003';
%
%      %
%      % Number of steps?
%      %
%      STEP = 75;
%
%      %
%      % Load the standard leapseconds and PCK kernels, and the MGS SPK
%      % kernel.
%      %
%      cspice_furnsh( 'standard.tm' )
%      cspice_furnsh( MGS )
%
%      %
%      % Convert the strings to ephemeris time J2000.
%      %
%      et_start = cspice_str2et( START_TIME );
%      et_stop = cspice_str2et( STOP_TIME );
%
%      %
%      % Length of a step in seconds for STEP steps.
%      %
%      space = (et_stop - et_start)/STEP;
%
%      %
%      % Create a vector of ephemeris times.
%      %
%      et = [0:(STEP-1)]*space + et_start;
%
%      %
%      % Start at 'et_start', take STEP steps
%      % of space 'length'. At each time, calculate the
%      % intercept sub-point of the observer, then calculate
%      % the illumination angles at the sub-point.
%      %
%      [pos, alt] = cspice_subpt( 'Intercept', TARGET, et, CORRECT, OBSERVER );
%
%      [ phase, solar, emissn] = cspice_illum( TARGET, et, CORRECT, ...
%                                               OBSERVER, pos );
%
%      %
%      % Convert the et value to UTC for human comprehension.
%      %
%      utc    = cspice_et2utc( et, 'C', 3 );
%      phase  = phase  * cspice_dpr;
%      solar  = solar  * cspice_dpr;
%      emissn = emissn * cspice_dpr;
%
%      for i = 1:STEP
%
%         %
%         % Output the times and lighting angles in degrees.
%         %
%         txt = sprintf( 'UTC           : %s', utc(i,:) );
%         disp( txt )
%
%         txt = sprintf( 'Emission angle: %14.6f', emissn(i) );
%         disp( txt )
%
%         txt = sprintf( 'Solar angle   : %14.6f', solar(i)  );
%         disp( txt )
%
%         txt = sprintf( 'Phase angle   : %14.6f', phase(i)  );
%         disp( txt )
%
%         disp( ' ' )
%      end
%
%   MATLAB outputs:
%
%      UTC           : 2003 FEB 01 00:00:00.000
%      Emission angle:     0.128603
%      Solar angle   :    49.308223
%      Phase angle   :    49.230858
%
%      UTC           : 2003 FEB 01 18:52:48.000
%      Emission angle:     0.301783
%      Solar angle   :   143.954297
%      Phase angle   :   144.017336
%
%                   ...
%
%      UTC           : 2003 MAR 30 10:14:24.000
%      Emission angle:       0.160300
%      Solar angle   :     139.575312
%      Phase angle   :     139.585791
%
%      UTC           : 2003 MAR 31 05:07:12.000
%      Emission angle:       0.291396
%      Solar angle   :      58.922414
%      Phase angle   :      58.711270
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine illum_c.
%
%   MICE.REQ
%   KERNEL.REQ
%   NAIF_IDS.REQ
%   SPK.REQ
%   TIME.REQ
%
%-Version
%
%   -Mice Version 1.0.3, 18-NOV-2014, EDW (JPL)
%
%       Edited I/O section to conform to NAIF standard for Mice documentation.
%
%   -Mice Version 1.0.2, 18-MAY-2010, BVS (JPL)
%
%      Index lines now state that this routine is deprecated.
%
%   -Mice Version 1.0.1, 30-DEC-2008, EDW (JPL)
%
%      Edits to header; Abstract now states that this routine is
%      deprecated.
%
%      Corrected misspellings.
%
%    -Mice Version 1.0.0, 15-DEC-2005, EDW (JPL)
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
   catch
      rethrow(lasterror)
   end




