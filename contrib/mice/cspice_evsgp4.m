%-Abstract
%
%   CSPICE_EVSGP4 evaluates NORAD two-line element data for earth orbiting
%   spacecraft. This evaluator uses algorithms as described
%   in Vallado 2006 [4].
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
%      et       the epoch in seconds past ephemeris epoch J2000 at which a
%               state should be produced from the input elements.
%
%               [1,1] = size(et); double = class(et)
%
%      geophs   a collection of 8 geophysical constants needed for computing
%               a state.
%
%               [8,1] = size(geophs); double = class(geophs)
%
%               The order of these constants must be:
%
%                  geophs(1) = j2 gravitational harmonic for Earth.
%                  geophs(2) = j3 gravitational harmonic for Earth.
%                  geophs(3) = j4 gravitational harmonic for Earth.
%
%               These first three constants are dimensionless.
%
%                  geophs(4) = KE: Square root of the GM for Earth where
%                              GM is expressed in Earth radii cubed per
%                              minutes squared.
%
%                  geophs(5) = QO: High altitude bound for atmospheric
%                              model in km.
%
%                  geophs(6) = SO: Low altitude bound for atmospheric
%                              model in km.
%
%                  geophs(7) = RE: Equatorial radius of the earth in km.
%
%                  geophs(8) = AE: Distance units/earth radius
%                              (normally 1)
%
%               Below are currently recommended values for these
%               items:
%
%                  J2 =    1.082616e-3
%                  J3 =   -2.53881e-6
%                  J4 =   -1.65597e-6
%
%               The next item is the square root of GM for the Earth
%               given in units of earth-radii**1.5/Minute
%
%                  KE =    7.43669161e-2
%
%               The next two items define the top and bottom of the
%               atmospheric drag model used by the type 10 ephemeris
%               type. Don't adjust these unless you understand the full
%               implications of such changes.
%
%                  QO =  120.0e0
%                  SO =   78.0e0
%
%               The ER value is the equatorial radius in km of the Earth
%               as used by NORAD.
%
%                  ER = 6378.135e0
%
%               The value of AE is the number of distance units per
%               Earth radii used by the NORAD state propagation
%               software. The value should be 1 unless you've got a very
%               good understanding of the NORAD routine SGP4 and the
%               affect of changing this value.
%
%                  AE =    1.0e0
%
%      elems    an array containing two-line element data as prescribed
%               below.
%
%               [10,1] = size(elems); double = class(elems)
%
%               The elements NDD6O and BSTAR must already be scaled by the
%               proper exponent stored in the two line elements set.
%               Moreover, the various items must be converted to the units
%               shown here.
%
%                  elems(  1 ) = NDT20 in radians/minute**2
%                  elems(  2 ) = NDD60 in radians/minute**3
%                  elems(  3 ) = BSTAR
%                  elems(  4 ) = INCL  in radians
%                  elems(  5 ) = NODE0 in radians
%                  elems(  6 ) = ECC
%                  elems(  7 ) = OMEGA in radians
%                  elems(  8 ) = M0    in radians
%                  elems(  9 ) = N0    in radians/minute
%                  elems( 10 ) = `epoch' of the elements in seconds
%                                past ephemeris epoch J2000.
%
%   the call:
%
%      [state] = cspice_evsgp4( et, geophs, elems )
%
%   returns:
%
%      state    the state produced by evaluating the input elements at the
%               input epoch `et'.
%
%               [6,1] = size(state); double = class(state)
%
%               Units are km and km/sec relative to the TEME reference
%               frame.
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
%   1) Suppose you have a set of two-line elements for the LUME-1
%      cubesat. This example shows how you can use this routine
%      together with the routine cspice_getelm to propagate a state to an
%      epoch of interest.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: evsgp4_ex1.tm
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
%            File name           Contents
%            ---------           ------------------------------------
%            naif0012.tls        Leapseconds
%            geophysical.ker     geophysical constants for evaluation
%                                of two-line element sets.
%
%         The geophysical.ker is a PCK file that is provided with the
%         Mice toolkit under the "/data" directory.
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                'geophysical.ker'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function evsgp4_ex1()
%
%         %
%         % Local parameters.
%         %
%         TIMSTR =   '2020-05-26 02:25:00';
%
%         %
%         % The LUME-1 cubesat is an Earth orbiting object; set
%         % the center ID to the Earth ID.
%         %
%         CENTER =   399;
%
%         %
%         % Local variables.
%         %
%         geophs = zeros(8,1);
%
%         %
%         % These are the variables that will hold the constants
%         % required by cspice_evsgp4. These constants are available
%         % from the loaded PCK file, which provides the actual
%         % values and units as used by NORAD propagation model.
%         %
%         %    Constant   Meaning
%         %    --------   ------------------------------------------
%         %    J2         J2 gravitational harmonic for Earth.
%         %    J3         J3 gravitational harmonic for Earth.
%         %    J4         J4 gravitational harmonic for Earth.
%         %    KE         Square root of the GM for Earth.
%         %    QO         High altitude bound for atmospheric model.
%         %    SO         Low altitude bound for atmospheric model.
%         %    ER         Equatorial radius of the Earth.
%         %    AE         Distance units/earth radius.
%         %
%         noadpn = {'J2','J3','J4','KE','QO','SO','ER','AE'};
%
%         %
%         % Define the Two-Line Element set for LUME-1.
%         %
%         tle = [ '1 43908U 18111AJ  20146.60805006  .00000806'            ...
%                                  '  00000-0  34965-4 0  9999',
%                 '2 43908  97.2676  47.2136 0020001 220.6050 '            ...
%                                  '139.3698 15.24999521 78544' ];
%
%         %
%         % Load the MK file that includes the PCK file that provides
%         % the geophysical constants required for the evaluation of
%         % the two-line elements sets and the LSK, as it is required
%         % by cspice_getelm to perform time conversions.
%         %
%         cspice_furnsh( 'evsgp4_ex1.tm' );
%
%         %
%         % Retrieve the data from the kernel, and place it on
%         % the `geophs' array.
%         %
%         for i=1:8
%
%            [geophs(i)] = cspice_bodvcd( CENTER, noadpn(i), 1 );
%
%         end
%
%         %
%         % Convert the Two Line Elements lines to the element sets.
%         % Set the lower bound for the years to be the beginning
%         % of the space age.
%         %
%         [epoch, elems] = cspice_getelm( 1957, tle );
%
%         %
%         % Now propagate the state using cspice_evsgp4 to the epoch
%         % of interest.
%         %
%         [et]    = cspice_str2et( TIMSTR );
%         [state] = cspice_evsgp4( et, geophs, elems );
%
%         %
%         % Display the results.
%         %
%         fprintf( 'Epoch   : %s\n', TIMSTR )
%         fprintf( 'Position: %15.8f %15.8f %15.8f\n',                     ...
%                         state(1), state(2), state(3) )
%         fprintf( 'Velocity: %15.8f %15.8f %15.8f\n',                     ...
%                         state(4), state(5), state(6) )
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Matlab due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a PC/Linux/Matlab9.x/32-bit
%      platform, the output was:
%
%
%      Epoch   : 2020-05-26 02:25:00
%      Position:  -4644.60403398  -5038.95025539   -337.27141116
%      Velocity:     -0.45719025      0.92884817     -7.55917355
%
%
%-Particulars
%
%   This routine evaluates any NORAD two-line element sets for
%   near-earth orbiting satellites using the algorithms described in
%   Vallado 2006 [4].
%
%-Exceptions
%
%   1)  No checks are made on the reasonableness of the inputs.
%
%   2)  If a problem occurs when evaluating the elements, an
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   3)  If any of the input arguments, `et', `geophs' or `elems', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   4)  If any of the input arguments, `et', `geophs' or `elems', is
%       not of the expected type, or it does not have the expected
%       dimensions and size, an error is signaled by the Mice
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
%   [1]  F. Hoots and R. Roehrich, "Spacetrack Report #3: Models for
%        Propagation of the NORAD Element Sets," U.S. Air Force
%        Aerospace Defense Command, Colorado Springs, CO, 1980.
%
%   [2]  F. Hoots, "Spacetrack Report #6: Models for Propagation of
%        Space Command Element Sets,"  U.S. Air Force Aerospace
%        Defense Command, Colorado Springs, CO, 1986.
%
%   [3]  F. Hoots, P. Schumacher and R. Glover, "History of Analytical
%        Orbit Modeling in the U. S. Space Surveillance System,"
%        Journal of Guidance, Control, and Dynamics. 27(2):174-185,
%        2004.
%
%   [4]  D. Vallado, P. Crawford, R. Hujsak and T. Kelso, "Revisiting
%        Spacetrack Report #3," paper AIAA 2006-6753 presented at the
%        AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO.,
%        August 21-24, 2006.
%
%-Author_and_Institution
%
%   M. Costa Sitja      (JPL)
%
%-Version
%
%   -Mice Version 1.0.0, 02-NOV-2021 (MCS)
%
%-Index_Entries
%
%   Evaluate NORAD two-line element data using SGP4.
%
%-&
function [state] = cspice_evsgp4( et, geophs, elems )

   switch nargin
      case 3

         et     = zzmice_dp(et);
         geophs = zzmice_dp(geophs);
         elems  = zzmice_dp(elems);

      otherwise

         error ( ['Usage: [_starg(6)_] = ' ...
                  'cspice_spkezr( et, _geophs(8)_, _elems(10_) )'] )


   end

   %
   % Call the MEX library.
   %
   try
      [state] = mice( 'evsgp4_c', et, geophs, elems );
   catch spiceerr
      rethrow(spiceerr)
   end
