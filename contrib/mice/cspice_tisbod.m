%-Abstract
%
%   CSPICE_TISBOD returns a 6x6 matrix that transforms states in inertial
%   coordinates to states in body-equator-and-prime-meridian
%   coordinates.
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
%      ref      the NAIF name for an inertial reference frame.
%
%               [1,c1] = size(ref); char = class(ref)
%
%                  or
%
%               [1,1] = size(ref); cell = class(ref)
%
%               Acceptable names include:
%
%                  Name       Description
%                  --------   --------------------------------
%                  'J2000'    Earth mean equator, dynamical
%                             equinox of J2000
%
%                  'B1950'    Earth mean equator, dynamical
%                             equinox of B1950
%
%                  'FK4'      Fundamental Catalog (4)
%
%                  'DE-118'   JPL Developmental Ephemeris (118)
%
%                  'DE-96'    JPL Developmental Ephemeris ( 96)
%
%                  'DE-102'   JPL Developmental Ephemeris (102)
%
%                  'DE-108'   JPL Developmental Ephemeris (108)
%
%                  'DE-111'   JPL Developmental Ephemeris (111)
%
%                  'DE-114'   JPL Developmental Ephemeris (114)
%
%                  'DE-122'   JPL Developmental Ephemeris (122)
%
%                  'DE-125'   JPL Developmental Ephemeris (125)
%
%                  'DE-130'   JPL Developmental Ephemeris (130)
%
%                  'GALACTIC' Galactic System II
%
%                  'DE-200'   JPL Developmental Ephemeris (200)
%
%                  'DE-202'   JPL Developmental Ephemeris (202)
%
%               See the Frames Required Reading frames.req for a full
%               list of inertial reference frame names built into
%               SPICE.
%
%               The output `tsipm' will give the transformation
%               from this frame to the bodyfixed frame specified by
%               `body' at the epoch specified by `et'.
%
%      body     the integer ID code of the body for which the state
%               transformation matrix is requested.
%
%               [1,1] = size(body); int32 = class(body)
%
%               Bodies are numbered according to the standard NAIF
%               numbering scheme. The numbering scheme is explained in the
%               NAIF IDs Required Reading naif_ids.req.
%
%      et       the epoch at which the state transformation matrix is
%               requested.
%
%               [1,1] = size(et); double = class(et)
%
%               (This is typically the epoch of observation minus the
%               one-way light time from the observer to the body at the epoch
%               of observation.)
%
%   the call:
%
%      [tsipm] = cspice_tisbod( ref, body, et )
%
%   returns:
%
%      tsipm    a 6x6 transformation matrix.
%
%               [6,6] = size(tsipm); double = class(tsipm)
%
%               It is used to transform states from inertial coordinates to
%               body fixed (also called equator and prime meridian --- PM)
%               coordinates.
%
%               Given a state `s' in the inertial reference frame
%               specified by `ref', the corresponding bodyfixed state
%               is given by the matrix vector product:
%
%                  tsipm * s
%
%               The X axis of the PM system is directed  to the
%               intersection of the equator and prime meridian.
%               The Z axis points along  the spin axis and points
%               towards the same side of the invariable plane of
%               the solar system as does earth's north pole.
%
%               NOTE: The inverse of `tsipm' is NOT its transpose.
%               The matrix, `tsipm', has a structure as shown below:
%
%                  .-            -.
%                  |       :      |
%                  |   r   :  0   |
%                  | ......:......|
%                  |       :      |
%                  | dr/dt :  r   |
%                  |       :      |
%                  `-            -'
%
%               where `r' is a time varying rotation matrix and dr/dt is
%               its derivative. The inverse of this matrix is:
%
%                  .-              -.
%                  |     T  :       |
%                  |    r   :  0    |
%                  | .......:.......|
%                  |        :       |
%                  |      T :   T   |
%                  | dr/dt  :  r    |
%                  |        :       |
%                  `-              -'
%
%               The Mice routine cspice_invstm is available for producing
%               this inverse.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for these examples may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%   1) Calculate the matrix to transform a state vector from the
%      J2000 frame to the Saturn fixed frame at a specified
%      time, and use it to compute the geometric position and
%      velocity of Titan in Saturn's body-fixed frame.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: tisbod_ex1.tm
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
%            sat375.bsp                    Saturn satellite ephemeris
%            pck00010.tpc                  Planet orientation and
%                                          radii
%            naif0012.tls                  Leapseconds
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'sat375.bsp',
%                                'pck00010.tpc',
%                                'naif0012.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function tisbod_ex1()
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'tisbod_ex1.tm' );
%
%         %
%         % The body ID for Saturn.
%         %
%         satid = 699;
%
%         %
%         % Retrieve the transformation matrix at some time.
%         %
%         [et]    = cspice_str2et( 'Jan 1 2005' );
%         [tsipm] = cspice_tisbod( 'J2000', satid, et );
%
%         %
%         % Retrieve the state of Titan as seen from Saturn
%         % in the J2000 frame at `et'.
%         %
%         [state, lt] = cspice_spkezr( 'TITAN', et,       'J2000',         ...
%                                      'NONE',  'SATURN'           );
%
%         fprintf( 'Titan as seen from Saturn (J2000 frame):\n' )
%         fprintf( '   position   (km): %12.3f %12.3f %12.3f\n', state(1:3) )
%         fprintf( '   velocity (km/s): %12.3f %12.3f %12.3f\n', state(4:6) )
%
%         %
%         % Rotate the 6-vector `state' into the
%         % Saturn body-fixed reference frame.
%         %
%         satvec = tsipm * state;
%
%         fprintf( 'Titan as seen from Saturn (IAU_SATURN frame):\n' )
%         fprintf( '   position   (km): %12.3f %12.3f %12.3f\n', satvec(1:3) )
%         fprintf( '   velocity (km/s): %12.3f %12.3f %12.3f\n', satvec(4:6) )
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
%      Titan as seen from Saturn (J2000 frame):
%         position   (km):  1071928.661  -505781.970   -60383.976
%         velocity (km/s):        2.404        5.176       -0.560
%      Titan as seen from Saturn (IAU_SATURN frame):
%         position   (km):   401063.338 -1116965.364    -5408.806
%         velocity (km/s):     -177.547      -63.745        0.028
%
%
%      Note that the complete example could be replaced by a single
%      cspice_spkezr call:
%
%         [state, lt] = cspice_spkezr( 'TITAN',      et,                ...
%                                      'IAU_SATURN', 'NONE',            ...
%                                      'SATURN'              );
%
%
%   2) Use cspice_tisbod is used to compute the angular velocity vector (with
%      respect to the J2000 inertial frame) of the specified body at
%      given time.
%
%      Use the meta-kernel from Example 1 above.
%
%
%      Example code begins here.
%
%
%      function tisbod_ex2()
%
%         %
%         % Local variables
%         %
%         dtipm = zeros(3,3);
%         tipm  = zeros(3,3);
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'tisbod_ex1.tm' );
%
%         %
%         % The body ID for Saturn.
%         %
%         satid = 699;
%
%         %
%         % First get the state transformation matrix.
%         %
%         [et]    = cspice_str2et( 'Jan 1 2005' );
%         [tsipm] = cspice_tisbod( 'J2000', satid, et );
%
%         %
%         % This matrix has the form:
%         %
%         %      .-            -.
%         %      |       :      |
%         %      | tipm  :  0   |
%         %      | ......:......|
%         %      |       :      |
%         %      | dtipm : tipm |
%         %      |       :      |
%         %      `-            -'
%         %
%         % We extract `tipm' and `dtipm'
%         %
%         for i=1:3
%
%            for j=1:3
%
%               tipm(i,j)  = tsipm(i,j);
%               dtipm(i,j) = tsipm(i+3,j);
%
%            end
%
%         end
%
%         %
%         % The transpose of `tipm' and `dtipm', (`tpmi' and `dtpmi'), gives
%         % the transformation from bodyfixed coordinates to inertial
%         % coordinates.
%         %
%         % Here is a fact about the relationship between angular
%         % velocity associated with a time varying rotation matrix
%         % that gives the orientation of a body with respect to
%         % an inertial frame.
%         %
%         %    The angular velocity vector can be read from the off
%         %    diagonal components of the matrix product:
%         %
%         %                            t
%         %    omega =     dtpmi * tpmi
%         %
%         %                     t
%         %          =     dtipm * tipm
%         %
%         %    the components of the angular velocity `v' will appear
%         %    in this matrix as:
%         %
%         %        .-                   -.
%         %        |                     |
%         %        |   0    -v(3)  v(2)  |
%         %        |                     |
%         %        |  v(3)    0   -v(1)  |
%         %        |                     |
%         %        | -v(2)   v(1)   0    |
%         %        |                     |
%         %        `-                   -'
%         %
%         %
%         omega = dtipm' * tipm;
%
%         v = [ omega(3,2), omega(1,3), omega(2,1) ]';
%
%         %
%         % Display the results.
%         %
%         fprintf( 'Angular velocity (km/s):\n' )
%         fprintf( '%16.9f %15.9f %15.9f\n', v )
%
%         %
%         % It is possible to compute the angular velocity using
%         % a single call to cspice_xf2rav.
%         %
%         [rot, av] = cspice_xf2rav( tsipm );
%
%         fprintf( 'Angular velocity using cspice_xf2rav (km/s):\n' )
%         fprintf( '%16.9f %15.9f %15.9f\n', av )
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
%      Angular velocity (km/s):
%           0.000014001     0.000011995     0.000162744
%      Angular velocity using cspice_xf2rav (km/s):
%           0.000014001     0.000011995     0.000162744
%
%
%-Particulars
%
%   Note: NAIF recommends the use of cspice_spkezr with the appropriate
%   frames kernels when possible over cspice_tisbod.
%
%   The matrix for transforming inertial states to bodyfixed
%   states is the 6x6 matrix shown below as a block structured
%   matrix.
%
%      .-            -.
%      |       :      |
%      | tipm  :  0   |
%      | ......:......|
%      |       :      |
%      | dtipm : tipm |
%      |       :      |
%      `-            -'
%
%   This can also be expressed in terms of Euler angles
%   `phi', `delta' and `w'. The transformation from inertial to
%   bodyfixed coordinates is represented in the SPICE kernel
%   pool as:
%
%      tipm   = [w]  [delta]  [phi]
%                  3        1      3
%   Thus
%
%      dtipm  = d[w] /dt [delta]  [phi]
%                   3           1      3
%
%             + [w]  d[delta] /dt  [phi]
%                  3             1      3
%
%             + [w]  [delta]  d[phi] /dt
%                  3        1           3
%
%
%   If a binary PCK file record can be used for the time and
%   body requested, it will be used. The most recently loaded
%   binary PCK file has first priority, followed by previously
%   loaded binary PCK files in backward time order. If no
%   binary PCK file has been loaded, the text P_constants
%   kernel file is used.
%
%   If there is only text PCK kernel information, it is
%   expressed in terms of `ra', `dec' and `w', where
%
%      ra  = phi - cspice_halfpi
%      dec = cspice_halfpi - delta
%      w   = w
%
%   The angles `ra', `dec', and `w' are defined as follows in the
%   text PCK file:
%
%                                    2    .-----
%                    ra1*t      ra2*t      \
%      ra  = ra0  + -------  + -------   +  )  a(i) * sin( theta(i) )
%                      T          2        /
%                                T        '-----
%                                            i
%
%                                     2   .-----
%                    dec1*t     dec2*t     \
%      dec = dec0 + -------- + --------  +  )  d(i) * cos( theta(i) )
%                      T           2       /
%                                 T       '-----
%                                            i
%
%                                   2     .-----
%                     w1*t      w2*t       \
%      w   = w0   +  ------  + -------   +  )  w(i) * sin( theta(i) )
%                      d          2        /
%                                d        '-----
%                                            i
%
%
%   where `d' is in seconds/day; T in seconds/Julian century;
%   a(i), d(i), and w(i) arrays apply to satellites only; and
%   theta(i), defined as
%
%                              theta1(i)*t
%      theta(i) = theta0(i) + -------------
%                                   T
%
%   are specific to each planet.
%
%   These angles ---typically nodal rates--- vary in number and
%   definition from one planetary system to the next.
%
%   Thus
%
%                                 .-----
%               ra1     2*ra2*t    \   a(i)*theta1(i)*cos(theta(i))
%   dra/dt   = ----- + --------- +  ) ------------------------------
%                T          2      /                 T
%                          T      '-----
%                                    i
%
%                                   .-----
%               dec1     2*dec2*t    \   d(i)*theta1(i)*sin(theta(i))
%    ddec/dt = ------ + ---------- -  ) ------------------------------
%                 T          2       /                 T
%                           T       '-----
%                                      i
%
%                               .-----
%               w1     2*w2*t    \   w(i)*theta1(i)*cos(theta(i))
%    dw/dt   = ---- + -------- +  ) ------------------------------
%               d         2      /                 T
%                        d      '-----
%                                  i
%
%-Exceptions
%
%   1)  If data required to define the body-fixed frame associated
%       with `body' are not found in the binary PCK system or the kernel
%       pool, the error SPICE(FRAMEDATANOTFOUND) is signaled by a
%       routine in the call tree of this routine. In the case of IAU
%       style body-fixed frames, the absence of prime meridian
%       polynomial data (which are required) is used as an indicator
%       of missing data.
%
%   2)  If the test for exception (1) passes, but in fact requested
%       data are not available in the kernel pool, an error is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If the kernel pool does not contain all of the data required
%       to define the number of nutation precession angles
%       corresponding to the available nutation precession
%       coefficients, the error SPICE(INSUFFICIENTANGLES) is
%       signaled by a routine in the call tree of this routine.
%
%   4)  If the reference frame `ref' is not recognized, an error is
%       signaled by a routine in the call tree of this routine.
%
%   5)  If the specified body code `body' is not recognized, an error is
%       signaled by a routine in the call tree of this routine.
%
%   6)  If, for a given body, both forms of the kernel variable names
%
%          BODY<body ID>_CONSTANTS_JED_EPOCH
%          BODY<body ID>_CONSTS_JED_EPOCH
%
%       are found in the kernel pool, the error
%       SPICE(COMPETINGEPOCHSPEC) is signaled by a routine in the call
%       tree of this routine. This is done regardless of whether the
%       values assigned to the kernel variable names match.
%
%   7)  If, for a given body, both forms of the kernel variable names
%
%          BODY<body ID>_CONSTANTS_REF_FRAME
%          BODY<body ID>_CONSTS_REF_FRAME
%
%       are found in the kernel pool, the error
%       SPICE(COMPETINGFRAMESPEC) is signaled by a routine in the call
%       tree of this routine. This is done regardless of whether the
%       values assigned to the kernel variable names match.
%
%   8)  If the central body associated with the input `body', whether a
%       system barycenter or `body' itself, has associated phase angles
%       (aka nutation precession angles), and the kernel variable
%       BODY<body ID>_MAX_PHASE_DEGREE for the central body is present
%       but has a value outside the range 1:3, the error
%       SPICE(DEGREEOUTOFRANGE) is signaled by a routine in the call
%       tree of this routine.
%
%   9)  If any of the input arguments, `ref', `body' or `et', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   10) If any of the input arguments, `ref', `body' or `et', is not
%       of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  The kernel pool must be loaded with the appropriate
%       coefficients (from a text or binary PCK file) prior to calling
%       this routine.
%
%-Required_Reading
%
%   FRAMES.REQ
%   MICE.REQ
%   NAIF_IDS.REQ
%   PCK.REQ
%   ROTATION.REQ
%   TIME.REQ
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
%   -Mice Version 1.0.0, 07-AUG-2021 (JDR)
%
%-Index_Entries
%
%   transformation from inertial state to bodyfixed
%
%-&
function [tsipm] = cspice_tisbod( ref, body, et )

   switch nargin
      case 3

         ref  = zzmice_str(ref);
         body = zzmice_int(body);
         et   = zzmice_dp(et);

      otherwise

         error ( 'Usage: [tsipm(6,6)] = cspice_tisbod( `ref`, body, et )' )

   end

   %
   % Call the MEX library.
   %
   try
      [tsipm] = mice('tisbod_c', ref, body, et);
   catch spiceerr
      rethrow(spiceerr)
   end
