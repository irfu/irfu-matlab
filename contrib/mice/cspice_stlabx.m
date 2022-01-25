%-Abstract
%
%   CSPICE_STLABX corrects the position of a target for the stellar
%   aberration effect on radiation transmitted from a specified observer to
%   the target.
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
%      pobj     the cartesian position vector of an object with respect to
%               the observer, possibly corrected for light time.
%
%               [3,1] = size(pobj); double = class(pobj)
%
%               Units are km.
%
%      vobs     the cartesian velocity vector of the observer with respect to
%               the Solar System barycenter.
%
%               [3,1] = size(vobs); double = class(vobs)
%
%               Units are km/s.
%
%   the call:
%
%      [corpos] = cspice_stlabx( pobj, vobs )
%
%   returns:
%
%      corpos   the position of the object relative to the observer,
%               corrected for the stellar aberration effect on radiation
%               directed toward the target.
%
%               [3,1] = size(corpos); double = class(corpos)
%
%               This correction is the inverse of the usual stellar
%               aberration correction: the corrected vector indicates the
%               direction in which radiation must be emitted from the
%               observer, as seen in an inertial reference frame having
%               velocity equal to that of the observer, in order to reach the
%               position indicated by the input vector `pobj'.
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
%   1) Compute the apparent position of the Moon relative to the
%      Earth, corrected for one way light-time and stellar aberration
%      effect on radiation transmitted from the Earth to the Moon,
%      given the geometric state of the Earth relative to the Solar
%      System Barycenter, and the difference between the stelar
%      aberration corrected and uncorrected position vectors, taking
%      several steps.
%
%      First, compute the light-time corrected state of the Moon body
%      as seen by the Earth, using its geometric state. Then apply
%      the correction for stellar aberration to the light-time
%      corrected state of the target body, both for the transmission
%      case.
%
%      The code in this example could be replaced by a single call
%      to cspice_spkpos:
%
%          [pos, lt] = cspice_spkpos( 'MOON',  et,      ...
%                                     'J2000', 'XLT+S', ...
%                                     'EARTH'           );
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: stlabx_ex1.tm
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
%            de418.bsp                     Planetary ephemeris
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de418.bsp',
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
%      function stlabx_ex1()
%
%         %
%         % Assign an observer, Earth, target, Moon, time of interest
%         % and reference frame for returned vectors.
%         %
%         idobs  = 399;
%         idtarg = 301;
%         utcstr = 'July 4 2004';
%         reffrm = 'J2000';
%
%         %
%         % Load the needed kernels.
%         %
%         cspice_furnsh( 'stlabx_ex1.tm' );
%
%         %
%         % Convert the time string to ephemeris time.
%         %
%         [et] = cspice_str2et( utcstr );
%
%         %
%         % Get the state of the observer with respect to the solar
%         % system barycenter.
%         %
%         [sobs] = cspice_spkssb( idobs, et, reffrm );
%
%         %
%         % Get the light-time corrected position `pos' of the target
%         % body `idtarg' as seen by the observer. Normally we would
%         % call cspice_spkpos to obtain this vector, but we already have
%         % the state of the observer relative to the solar system
%         % barycenter, so we can avoid looking up that state twice
%         % by calling cspice_spkapo.
%         %
%         [pos, lt] = cspice_spkapo( idtarg, et, reffrm, sobs, 'XLT' );
%
%         %
%         % Output the uncorrected vector.
%         %
%         fprintf( 'Uncorrected position vector\n' )
%         fprintf( '    %18.6f %18.6f %18.6f\n', pos(1), pos(2), pos(3) )
%
%         %
%         % Apply the correction for stellar aberration to the
%         % light-time corrected position of the target body.
%         %
%         [pcorr] = cspice_stlabx( pos, sobs(4:6) );
%
%         %
%         % Output the corrected position vector and the apparent
%         % difference from the uncorrected vector.
%         %
%         fprintf( '\n' )
%         fprintf( 'Corrected position vector\n' )
%         fprintf( '    %18.6f %18.6f %18.6f\n', ...
%                   pcorr(1), pcorr(2), pcorr(3) )
%
%         %
%         % Apparente difference.
%         %
%         appdif = pos - pcorr;
%         fprintf( '\n' )
%         fprintf( 'Apparent difference\n' )
%         fprintf( '    %18.6f %18.6f %18.6f\n',   ...
%                  appdif(1), appdif(2), appdif(3) )
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
%      Uncorrected position vector
%               201809.933536     -260878.049826     -147716.077987
%
%      Corrected position vector
%               201782.730972     -260894.375627     -147724.405897
%
%      Apparent difference
%                   27.202563          16.325802           8.327911
%
%
%-Particulars
%
%   In order to transmit radiation from an observer to a specified
%   target, the emission direction must be corrected for one way
%   light time and for the motion of the observer relative to the
%   solar system barycenter. The correction for the observer's
%   motion when transmitting to a target is the inverse of the
%   usual stellar aberration correction applied to the light-time
%   corrected position of the target as seen by the observer.
%
%   Below is the description of the stellar aberration correction
%   used in the Mice routine cspice_stelab (with the notation changed
%   slightly):
%
%      Let `r' be the vector from the observer to the object, and `v' be
%      the velocity of the observer with respect to the Solar System
%      barycenter. Let `w' be the angle between them. The aberration
%      angle `phi' is given by
%
%         sin(phi) = v * sin(w) / C
%
%      Let `h' be the vector given by the cross product
%
%         h = r x v
%
%      Rotate `r' by `phi' radians about `h' to obtain the apparent position
%      of the object.
%
%   This routine applies the inverse correction, so here the rotation
%   about `h' is by -phi radians.
%
%-Exceptions
%
%   1)  If the velocity of the observer is greater than or equal to
%       the speed of light, an error is signaled by a routine in the
%       call tree of this routine. The outputs are undefined.
%
%   2)  If any of the input arguments, `pobj' or `vobs', is undefined,
%       an error is signaled by the Matlab error handling system.
%
%   3)  If any of the input arguments, `pobj' or `vobs', is not of the
%       expected type, or it does not have the expected dimensions and
%       size, an error is signaled by the Mice interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   None.
%
%-Required_Reading
%
%   MICE.REQ
%
%-Literature_References
%
%   [1]  W. Owen, "The Treatment of Aberration in Optical Navigation",
%        JPL IOM #314.8-524, 8 February 1985.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   stellar aberration for transmission case
%
%-&
function [corpos] = cspice_stlabx( pobj, vobs )

   switch nargin
      case 2

         pobj = zzmice_dp(pobj);
         vobs = zzmice_dp(vobs);

      otherwise

         error ( 'Usage: [corpos(3)] = cspice_stlabx( pobj(3), vobs(3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [corpos] = mice('stlabx_c', pobj, vobs);
   catch spiceerr
      rethrow(spiceerr)
   end
