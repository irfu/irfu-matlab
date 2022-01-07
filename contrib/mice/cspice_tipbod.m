%-Abstract
%
%   CSPICE_TIPBOD returns a 3x3 matrix that transforms positions in inertial
%   coordinates to positions in body-equator-and-prime-meridian
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
%               The output `tipm' will give the transformation
%               from this frame to the bodyfixed frame specified by
%               `body' at the epoch specified by `et'.
%
%      body     the integer ID code of the body for which the position
%               transformation matrix is requested.
%
%               [1,1] = size(body); int32 = class(body)
%
%               Bodies are numbered according to the standard NAIF
%               numbering scheme. The numbering scheme is explained in the
%               NAIF IDs Required Reading naif_ids.req.
%
%      et       the epoch at which the position transformation matrix is
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
%      [tipm] = cspice_tipbod( ref, body, et )
%
%   returns:
%
%      tipm     a 3x3 coordinate transformation matrix.
%
%               [3,3] = size(tipm); double = class(tipm)
%
%               It is used to transform positions from inertial coordinates
%               to body fixed (also called equator and prime meridian --- PM)
%               coordinates.
%
%               Given a position P in the inertial reference frame
%               specified by `ref', the corresponding bodyfixed
%               position is given by the matrix vector product:
%
%                  tipm * s
%
%               The X axis of the PM system is directed to the
%               intersection of the equator and prime meridian.
%               The Z axis points along  the spin axis and points
%               towards the same side of the invariable plane of
%               the solar system as does earth's north pole.
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
%   1) Calculate the matrix to rotate a position vector from the
%      J2000 frame to the Saturn fixed frame at a specified
%      time, and use it to compute the position of Titan in
%      Saturn's body-fixed frame.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: tipbod_ex1.tm
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
%      function tipbod_ex1()
%
%         %
%         % Load the kernels.
%         %
%         cspice_furnsh( 'tipbod_ex1.tm' );
%
%         %
%         % The body ID for Saturn.
%         %
%         satid = 699;
%
%         %
%         % Retrieve the transformation matrix at some time.
%         %
%         [et]   = cspice_str2et( 'Jan 1 2005' );
%         [tipm] = cspice_tipbod( 'J2000', satid, et );
%
%         %
%         % Retrieve the position of Titan as seen from Saturn
%         % in the J2000 frame at `et'.
%         %
%         [pos, lt] = cspice_spkpos( 'TITAN', et,       'J2000',           ...
%                                    'NONE',  'SATURN'           );
%
%         fprintf( 'Titan as seen from Saturn:\n' )
%         fprintf( '   in J2000 frame     : %12.3f %12.3f %12.3f\n', pos )
%
%         %
%         % Rotate the position 3-vector `pos' into the
%         % Saturn body-fixed reference frame.
%         %
%         satvec = tipm * pos;
%
%         fprintf( '   in IAU_SATURN frame: %12.3f %12.3f %12.3f\n', satvec )
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
%      Titan as seen from Saturn:
%         in J2000 frame     :  1071928.661  -505781.970   -60383.976
%         in IAU_SATURN frame:   401063.338 -1116965.364    -5408.806
%
%
%      Note that the complete example could be replaced by a single
%      cspice_spkpos call:
%
%         [pos, lt] = cspice_spkpos( 'TITAN', et,       'IAU_SATURN',   ...
%                                    'NONE',  'SATURN'                );
%
%-Particulars
%
%   cspice_tipbod takes PCK information as input, either in the
%   form of a binary or text PCK file. High precision
%   binary files are searched for first (the last loaded
%   file takes precedence); then it defaults to the text
%   PCK file. If binary information is found for the
%   requested body and time, the Euler angles are
%   evaluated and the transformation matrix is calculated
%   from them. Using the Euler angles `phi', `delta' and `w'
%   we compute
%
%      tipm = [w]  [delta]  [phi]
%                3        1      3
%
%   If no appropriate binary PCK files have been loaded,
%   the text PCK file is used. Here information is found
%   as `ra', `dec' and `w' (with the possible addition of nutation
%   and libration terms for satellites). Again, the Euler
%   angles are found, and the transformation matrix is
%   calculated from them. The transformation from inertial to
%   body-fixed coordinates is represented as:
%
%      tipm = [w]  [cspice_halfpi-dec]  [ra+cspice_halfpi]
%                3                    1                   3
%
%   These are basically the Euler angles, `phi', `delta' and `w':
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
%-Exceptions
%
%   1)  If the kernel pool does not contain all of the data required
%       for computing the transformation matrix, `tipm', the error
%       SPICE(INSUFFICIENTANGLES) is signaled by a routine in the call
%       tree of this routine.
%
%   2)  If the reference frame, `ref', is not recognized, an error is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If the specified body code, `body', is not recognized, an error
%       is signaled by a routine in the call tree of this routine.
%
%   4)  If any of the input arguments, `ref', `body' or `et', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   5)  If any of the input arguments, `ref', `body' or `et', is not
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
%       coefficients (from a text or binary PCK file) prior to
%       calling this routine.
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
%   -Mice Version 1.0.0, 07-SEP-2020 (JDR)
%
%-Index_Entries
%
%   transformation from inertial position to bodyfixed
%
%-&
function [tipm] = cspice_tipbod( ref, body, et )

   switch nargin
      case 3

         ref  = zzmice_str(ref);
         body = zzmice_int(body);
         et   = zzmice_dp(et);

      otherwise

         error ( 'Usage: [tipm(3,3)] = cspice_tipbod( `ref`, body, et )' )

   end

   %
   % Call the MEX library.
   %
   try
      [tipm] = mice('tipbod_c', ref, body, et);
   catch spiceerr
      rethrow(spiceerr)
   end
