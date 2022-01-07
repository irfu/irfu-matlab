%-Abstract
%
%   CSPICE_STELAB corrects the apparent position of an object for stellar
%   aberration.
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
%      pobj     the position (x, y, z, km) of an object with respect to the
%               observer, possibly corrected for light time.
%
%               [3,1] = size(pobj); double = class(pobj)
%
%      vobs     the velocity (dx/dt, dy/dt, dz/dt, km/sec) of the observer
%               with respect to the Solar System barycenter.
%
%               [3,1] = size(vobs); double = class(vobs)
%
%   the call:
%
%      [appobj] = cspice_stelab( pobj, vobs )
%
%   returns:
%
%      appobj   the apparent position of the object relative to the observer,
%               corrected for stellar aberration.
%
%               [3,1] = size(appobj); double = class(appobj)
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
%      Earth, corrected for one light-time and stellar aberration,
%      given the geometric state of the Earth relative to the Solar
%      System Barycenter, and the difference between the stellar
%      aberration corrected and uncorrected position vectors, taking
%      several steps.
%
%      First, compute the light-time corrected state of the Moon body
%      as seen by the Earth, using its geometric state. Then apply
%      the correction for stellar aberration to the light-time
%      corrected state of the target body.
%
%      The code in this example could be replaced by a single call
%      to cspice_spkpos:
%
%         [pos, lt] = cspice_spkpos( 'MOON',  et,     ...
%                                    'J2000', 'LT+S', ...
%                                    'EARTH'          );
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: stelab_ex1.tm
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
%      function stelab_ex1()
%
%         %
%         % Assign an observer, Earth, target, Moon, time of interest and
%         % reference frame for returned vectors.
%         %
%         idobs  = 399;
%         idtarg = 301;
%         utcstr = 'July 4 2004';
%         reffrm = 'J2000';
%
%         %
%         % Load the needed kernels.
%         %
%         cspice_furnsh( 'stelab_ex1.tm' );
%
%         %
%         % Convert the time string to ephemeris time, J2000.
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
%         % body `idtarg' as seen by the observer.
%         %
%         [pos, lt] = cspice_spkapo( idtarg, et, reffrm, sobs, 'LT' );
%
%         %
%         % Output the uncorrected vector.
%         %
%         fprintf( 'Uncorrected position vector\n' )
%         fprintf( '   %18.6f %18.6f %18.6f\n', pos(1), pos(2), pos(3) )
%
%         %
%         % Apply the correction for stellar aberration to the
%         % light-time corrected position of the target body.
%         %
%         [pcorr] = cspice_stelab( pos, sobs(4:6) );
%
%         %
%         % Output the corrected position vector and the apparent
%         % difference from the uncorrected vector.
%         %
%         fprintf( '\n' )
%         fprintf( 'Corrected position vector\n' )
%         fprintf( '   %18.6f %18.6f %18.6f\n', ...
%                  pcorr(1), pcorr(2), pcorr(3) )
%
%         %
%         % Apparent difference.
%         %
%         appdif = pos - pcorr;
%         fprintf( '\n' )
%         fprintf( 'Apparent difference\n' )
%         fprintf( '   %18.6f %18.6f %18.6f\n',    ...
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
%              201738.725087     -260893.141602     -147722.589056
%
%      Corrected position vector
%              201765.929516     -260876.818077     -147714.262441
%
%      Apparent difference
%                 -27.204429         -16.323525          -8.326615
%
%
%-Particulars
%
%   Let r be the vector from the observer to the object, and v be
%       -                                                    -
%   the velocity of the observer with respect to the Solar System
%   barycenter. Let w be the angle between them. The aberration
%   angle phi is given by
%
%        sin(phi) = v sin(w) / c
%
%   Let h be the vector given by the cross product
%       -
%
%         h = r X v
%         -   -   -
%
%   Rotate r by phi radians about h to obtain the apparent position
%          -                      -
%   of the object.
%
%-Exceptions
%
%   1)  If the velocity of the observer is greater than or equal
%       to the speed of light, the error SPICE(VALUEOUTOFRANGE)
%       is signaled by a routine in the call tree of this routine.
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
%   stellar aberration
%
%-&
function [appobj] = cspice_stelab( pobj, vobs )

   switch nargin
      case 2

         pobj = zzmice_dp(pobj);
         vobs = zzmice_dp(vobs);

      otherwise

         error ( 'Usage: [appobj(3)] = cspice_stelab( pobj(3), vobs(3) )' )

   end

   %
   % Call the MEX library.
   %
   try
      [appobj] = mice('stelab_c', pobj, vobs);
   catch spiceerr
      rethrow(spiceerr)
   end
