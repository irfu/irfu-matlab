%-Abstract
%
%   CSPICE_TWOVXF finds the state transformation from a base frame to the
%   right-handed frame defined by two state vectors: one state
%   vector defining a specified axis and a second state vector
%   defining a specified coordinate plane.
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
%      axdef    a "generalized" state vector defining one of the principal
%               axes of a reference frame.
%
%               [6,1] = size(axdef); double = class(axdef)
%
%               This vector consists of three components of a vector-valued
%               function of one independent variable `t' followed by the
%               derivatives of the components with respect to that variable:
%
%                  ( a, b, c, da/dt, db/dt, dc/dt )
%
%               This routine treats the input states as unitless, but
%               in most applications the input states represent
%               quantities that have associated units. The first three
%               components must have the same units, and the units of
%               the last three components must be compatible with
%               those of the first three: if the first three
%               components of `axdef'
%
%                  ( a, b, c )
%
%               have units U and `t' has units T, then the units of
%               `axdef' normally would be
%
%                  ( U, U, U, U/T, U/T, U/T )
%
%               Note that the direction and angular velocity defined
%               by `axdef' are actually independent of U, so scaling
%               `axdef' doesn't affect the output of this routine.
%
%               `axdef' could represent position and velocity; it could
%               also represent velocity and acceleration. `axdef' could
%               for example represent the velocity and acceleration of
%               a time-dependent position vector ( x(t), y(t), z(t) ),
%               in which case `axdef' would be defined by
%
%                  a     = dx/dt
%                  b     = dy/dt
%                  c     = dz/dt
%
%                           2      2
%                  da/dt = d x / dt
%
%                           2      2
%                  db/dt = d y / dt
%
%                           2      2
%                  dc/dt = d z / dt
%
%               Below, we'll call the normalized (unit length) version
%               of
%
%                  ( a, b, c )
%
%               the "direction" of `axdef'.
%
%               We call the frame relative to which `axdef' is specified
%               the "base frame." The input state `plndef' must be
%               specified relative to the same base frame.
%
%      indexa   the index of the reference frame axis that is parallel to the
%               direction of `axdef'.
%
%               [1,1] = size(indexa); int32 = class(indexa)
%
%                  indexa   Axis
%                  ------   ----
%                     1       X
%                     2       Y
%                     3       Z
%
%      plndef   a state vector defining (with `axdef') a principal plane of
%               the reference frame.
%
%               [6,1] = size(plndef); double = class(plndef)
%
%               This vector consists of three components followed by their
%               derivatives with respect to the independent variable `t'
%               associated with `axdef', so `plndef' is
%
%                  ( e, f, g, de/dt, df/dt, dg/dt )
%
%               Below, we'll call the unitized version of
%
%                  ( e, f, g )
%
%               the "direction" of `plndef'.
%
%               The second axis of the principal plane containing the
%               direction vectors of `axdef' and `plndef' is perpendicular
%               to the first axis and has positive dot product with
%               the direction vector of `plndef'.
%
%               The first three components of `plndef' must have the
%               same units, and the units of the last three components
%               must be compatible with those of the first three: if
%               the first three components of `plndef'
%
%                  ( e, f, g )
%
%               have units U2 and `t' has units T, then the units of
%               `plndef' normally would be
%
%                  ( U2, U2, U2, U2/T, U2/T, U2/T )
%
%               Note that ***for meaningful results, the angular
%               velocities defined by `axdef' and `plndef' must both have
%               units of 1/T.***
%
%               As with `axdef', scaling `plndef' doesn't affect the
%               output of this routine.
%
%               `axdef' and `plndef' must be specified relative to a
%               common reference frame, which we call the "base
%               frame."
%
%      indexp   the index of second axis of the principal frame determined by
%               `axdef' and `plndef'.
%
%               [1,1] = size(indexp); int32 = class(indexp)
%
%               The association of integer values and axes is the same as
%               for `indexa'.
%
%   the call:
%
%      [xform] = cspice_twovxf( axdef, indexa, plndef, indexp )
%
%   returns:
%
%      xform    the 6x6 matrix that transforms states from the frame relative
%               to which `axdef' and `plndef' are specified (the "base
%               frame") to the frame whose axes and derivative are determined
%               by `axdef', `plndef', `indexa' and `indexp'.
%
%               [6,6] = size(xform); double = class(xform)
%
%               The matrix `xform' has the structure shown below:
%
%                  .-              -.
%                  |        :       |
%                  |    r   :   0   |
%                  |        :       |
%                  | .......:.......|
%                  |        :       |
%                  |  dr/dt :   r   |
%                  |        :       |
%                  `-              -'
%
%               where `r' is a rotation matrix that is a function of
%               the independent variable associated with `axdef' and
%               `plndef', and where dr/dt is the derivative of `r'
%               with respect to that independent variable.
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
%   1) The time-dependent Sun-Canopus reference frame associated with
%      a spacecraft uses the spacecraft-sun state to define the Z axis
%      and the Canopus direction to define the X-Z plane.
%
%      Find the apparent position of the Earth as seen from the Mars
%      Reconnaissance Orbiter spacecraft (MRO) at a specified time,
%      relative to the Sun-Canopus reference frame associated with
%      MRO.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File: twovxf_ex1.tm
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
%            naif0012.tls                     Leapseconds
%            de430.bsp                        Planetary ephemeris
%            mro_psp4_ssd_mro95a.bsp          MRO ephemeris
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'de430.bsp',
%                                'mro_psp4_ssd_mro95a.bsp' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function twovxf_ex1()
%
%         %
%         % Local parameters
%         %
%         META =   'twovxf_ex1.tm';
%
%         %
%         % Local variables
%         %
%         stcano = zeros(6,1);
%
%         %
%         % Define the Right Ascension and Declination, and the
%         % proper motion in both coordinates, of Canopus, relative
%         % to the J2000 frame at J2000 epoch, in degrees and
%         % arcsecond/yr respectively. Note that the values used here
%         % may not be suitable for real applications.
%         %
%         RAJ2K  =   90.3991968556;
%         DECJ2K =   -52.6956610556;
%         PMRA   =   19.93e-3;
%         PMDEC  =   23.24e-3;
%
%         %
%         % Load kernel files via the meta-kernel.
%         %
%         cspice_furnsh( META );
%
%         %
%         % Convert the TDB input time string to seconds past
%         % J2000, TDB.
%         %
%         [et] = cspice_str2et( '2007 SEP 30 00:00:00 TDB' );
%
%         %
%         % Define an approximate "state vector" for Canopus using
%         % the J2000-relative, unit direction vector toward Canopus
%         % at a specified time `et' (time is needed to compute proper
%         % motion) as position and the zero vector as velocity.
%         %
%         [rpmra]  = cspice_convrt( PMRA, 'ARCSECONDS', 'RADIANS' );
%         [rpmdec] = cspice_convrt( PMDEC, 'ARCSECONDS', 'RADIANS' );
%
%         ra       = RAJ2K  * cspice_rpd + rpmra  * et/cspice_jyear;
%         dec      = DECJ2K * cspice_rpd + rpmdec * et/cspice_jyear;
%
%         [pcano]  = cspice_radrec( 1.0, ra, dec );
%
%         %
%         % Compute MRO geometric velocity w.r.t. the Solar System
%         % Barycenter, and use it to correct the Canopus direction
%         % for stellar aberration.
%         %
%         [state, lt] = cspice_spkezr( 'MRO', et, 'J2000', 'NONE', 'SSB' );
%
%         [stelab]    = cspice_stelab( pcano, state(4:6) );
%         stcano(1:3) = stelab';
%
%         stcano(4:6) = [ 0.0, 0.0, 0.0 ]';
%
%         %
%         % Let `stsun' be the J2000-relative apparent state of the Sun
%         % relative to the spacecraft at `et'.
%         %
%         [stsun, lt] = cspice_spkezr( 'SUN', et, 'J2000', 'CN+S', 'MRO' );
%
%         %
%         % The matrix `xfisc' transforms states from J2000 frame
%         % to the Sun-Canopus reference frame at `et'.
%         %
%         [xfisc] = cspice_twovxf( stsun, 3, stcano, 1 );
%
%         %
%         % Compute the apparent state of the Earth as seen from MRO
%         % in the J2000 frame at `et' and transform that vector into
%         % the Sun-Canopus reference frame.
%         %
%         [state, lt] = cspice_spkezr( 'EARTH', et, 'J2000', 'CN+S', 'MRO' );
%
%         sterth = xfisc * state;
%
%         %
%         % Display the results.
%         %
%         fprintf( [ 'Earth as seen from MRO in Sun-Canopus frame (km',    ...
%                    ' and km/s):\n' ]                                  )
%         fprintf( '   position: %15.3f %15.3f %15.3f\n', sterth(1:3) )
%         fprintf( '   velocity: %15.3f %15.3f %15.3f\n', sterth(4:6) )
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
%      Earth as seen from MRO in Sun-Canopus frame (km and km/s):
%         position:   -16659764.322    97343706.915   106745539.738
%         velocity:           2.691         -10.345          -7.877
%
%
%-Particulars
%
%   Given two linearly independent state vectors `axdef' and `plndef',
%   define vectors `dir1' and `dir2' by
%
%      dir1 = ( axdef(1),   axdef(2),   axdef(3)  )
%      dir2 = ( plndef(1),  plndef(2),  plndef(3) )
%
%   Then there is a unique right-handed reference frame `f' having:
%
%      `dir1' lying along the `indexa' axis.
%
%      `dir2' lying in the indexa-indexp coordinate plane, such that
%      the dot product of `dir2' with the positive `indexp' axis is
%      positive.
%
%   This routine determines the 6x6 matrix that transforms states
%   from the base frame used to represent the input vectors to the
%   the frame `f' determined by `axdef' and `plndef'. Thus a state vector
%
%      s       = ( x, y, z, dx/dt, dy/dt, dz/dt )
%       base
%
%   in the input reference frame will be transformed to
%
%      s       = xform * s
%       f                 base
%
%   in the frame `f' determined by `axdef' and `plndef'.
%
%-Exceptions
%
%   1)  If `indexa' or `indexp' is not in the set {1,2,3}, the error
%       SPICE(BADINDEX) is signaled by a routine in the call tree of
%       this routine.
%
%   2)  If `indexa' and `indexp' are the same, the error
%       SPICE(UNDEFINEDFRAME) is signaled by a routine in the call
%       tree of this routine.
%
%   3)  If the cross product of the vectors `axdef' and `plndef' is zero,
%       the error SPICE(DEPENDENTVECTORS) is signaled by a routine in
%       the call tree of this routine.
%
%   4)  If any of the input arguments, `axdef', `indexa', `plndef' or
%       `indexp', is undefined, an error is signaled by the Matlab
%       error handling system.
%
%   5)  If any of the input arguments, `axdef', `indexa', `plndef' or
%       `indexp', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
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
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%
%-Version
%
%   -Mice Version 1.0.0, 08-FEB-2021 (JDR)
%
%-Index_Entries
%
%   define a state transformation matrix from two states
%
%-&
function [xform] = cspice_twovxf( axdef, indexa, plndef, indexp )

   switch nargin
      case 4

         axdef  = zzmice_dp(axdef);
         indexa = zzmice_int(indexa);
         plndef = zzmice_dp(plndef);
         indexp = zzmice_int(indexp);

      otherwise

         error ( [ 'Usage: [xform(6,6)] = '                                 ...
                   'cspice_twovxf( axdef(6), indexa, plndef(6), indexp )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [xform] = mice('twovxf_c', axdef, indexa, plndef, indexp);
   catch spiceerr
      rethrow(spiceerr)
   end
