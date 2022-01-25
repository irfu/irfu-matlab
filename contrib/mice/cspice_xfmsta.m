%-Abstract
%
%   CSPICE_XFMSTA transforms a state between coordinate systems.
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
%      istate   state vector(s) in the input `icosys' coordinate system
%               representing position and velocity.
%
%               [6,n] = size(istate); double = class(istate)
%
%               All angular measurements must be in radians.
%
%               Note: body radii values taken from the kernel
%               pool are used when converting to or from geodetic or
%               planetographic coordinates. It is the user's
%               responsibility to verify the distance inputs are in
%               the same units as the radii in the kernel pool,
%               typically kilometers.
%
%      icosys   the name of the coordinate system that the input state vector
%               `istate' is currently in.
%
%               [1,c1] = size(icosys); char = class(icosys)
%
%                  or
%
%               [1,1] = size(icosys); cell = class(icosys)
%
%               `icosys' may be any of the following:
%
%                  'RECTANGULAR'
%                  'CYLINDRICAL'
%                  'LATITUDINAL'
%                  'SPHERICAL'
%                  'GEODETIC'
%                  'PLANETOGRAPHIC'
%
%               Leading spaces, trailing spaces, and letter case
%               are ignored. For example, ' cyLindRical  ' would be
%               accepted.
%
%      ocosys   the name of the coordinate system that the state should be
%               converted to.
%
%               [1,c2] = size(ocosys); char = class(ocosys)
%
%                  or
%
%               [1,1] = size(ocosys); cell = class(ocosys)
%
%               Please see the description of `icosys' for details.
%
%      body     the name or NAIF ID of the body associated with the
%               planetographic or geodetic coordinate system.
%
%               [1,c3] = size(body); char = class(body)
%
%                  or
%
%               [1,1] = size(body); cell = class(body)
%
%               If neither of the coordinate system choices are
%               geodetic or planetographic, `body' is ignored. It may
%               be a blank string.
%
%               Examples of accepted body names or IDs are:
%
%                  'Earth'
%                  '399'
%
%               Leading spaces, trailing spaces, and letter case are
%               ignored.
%
%   the call:
%
%      [ostate] = cspice_xfmsta( istate, icosys, ocosys, body )
%
%   returns:
%
%      ostate   the state vector(s) that ha(s|ve) been converted to the output
%               coordinate system `ocosys'.
%
%               [6,n] = size(ostate); double = class(ostate)
%
%               `ostate' returns with the same vectorization measure, N,
%               as `istate'.
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
%
%   1) Find the apparent state of Phoebe as seen by CASSINI in the
%      J2000 frame at three times starting at 2004 Jun 11 19:32:00.
%      Transform the state from rectangular to latitudinal coordinates.
%      For verification, transform the state back from latitudinal to
%      rectangular coordinates.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: xfmsta_ex1.tm
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
%            pck00010.tpc                     Planet orientation and
%                                             radii
%            naif0012.tls                     Leapseconds
%            041014R_SCPSE_01066_04199.bsp    CASSINI, planetary and
%                                             Saturn Satellite
%                                             ephemeris
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'naif0012.tls',
%                                '041014R_SCPSE_01066_04199.bsp',
%                                'pck00010.tpc'                 )
%
%         \begintext
%
%         End of meta-kernel.
%
%
%      Example code begins here.
%
%
%      function xfmsta_ex1()
%
%         metakr = 'xfmsta_ex1.tm';
%         form = 'YYYY-MM-DD HR:MN:SC (TDB) ::TDB';
%
%         %
%         %   Load the meta kernel.
%         %
%         cspice_furnsh ( metakr );
%
%         %
%         %   Calculate the state at 2004 Jun 11 19:32:00 UTC.
%         %
%         times = [ '2004-JUN-11-19:32:00';
%                   '2004-JUN-11-19:40:00';
%                   '2004-JUN-11-19:48:00'];
%         et = cspice_str2et ( times );
%
%         %
%         %   Calculate the apparent state of Phoebe as seen by
%         %   CASSINI in the J2000 frame.
%         %
%         [state_rec, lt] = cspice_spkezr ( 'phoebe', et, 'iau_phoebe', ...
%                                           'lt+s', 'cassini' );
%
%         %
%         %   Transform the state from rectangular to latitudinal.
%         %   Notice that since neither the input nor output
%         %   coordinate frames are 'geodetic' or 'planetographic',
%         %   the input for the body name is a blank string.
%         %
%         state_lat = cspice_xfmsta ( state_rec, 'rectangular', ...
%                                     'latitudinal', ' ' );
%
%         %
%         %   Transform the state back to rectangular from latitudinal. The
%         %   result should be very close to 'state_rec'.
%         %
%         state_rec2 = cspice_xfmsta ( state_lat, 'latitudinal', ...
%                                      'rectangular', ' ');
%
%         %
%         %   Report the results.
%         %
%         fprintf('\nPhoebe as seen by Cassini - rectangular\n')
%         for i = 1:length(et)
%             fprintf('Time: %s\n', cspice_timout ( et(i), form))
%             fprintf('  Position: %16.6f %16.6f %16.6f\n', state_rec(1:3,i))
%             fprintf('  Velocity: %16.6f %16.6f %16.6f\n', state_rec(4:6,i))
%         end
%
%         fprintf('\nPhoebe as seen by Cassini - latitudinal\n')
%         for i = 1:length(et)
%             fprintf('Time: %s\n', cspice_timout ( et(i), form))
%             fprintf('  Position: %16.6f %16.6f %16.6f\n', state_lat(1:3,i))
%             fprintf('  Velocity: %16.6f %16.6f %16.6f\n', state_lat(4:6,i))
%         end
%
%         fprintf('\nVerification: Phoebe as seen by Cassini - rectangular\n')
%         for i = 1:length(et)
%             fprintf('Time: %s\n', cspice_timout ( et(i), form))
%             fprintf('  Position: %16.6f %16.6f %16.6f\n', state_rec2(1:3,i))
%             fprintf('  Velocity: %16.6f %16.6f %16.6f\n', state_rec2(4:6,i))
%         end
%
%         %
%         %   Unload the kernels.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Phoebe as seen by Cassini - rectangular
%      Time: 2004-06-11 19:33:04 (TDB)
%        Position:     -2059.271283      -942.128329       -95.837672
%        Velocity:         3.910113        -4.228139        -1.526561
%      Time: 2004-06-11 19:41:04 (TDB)
%        Position:      -381.823839     -3135.942729      -828.567239
%        Velocity:         3.064009        -4.893554        -1.526501
%      Time: 2004-06-11 19:49:04 (TDB)
%        Position:       868.679897     -5620.615726     -1561.286108
%        Velocity:         2.133569        -5.438272        -1.526496
%
%      Phoebe as seen by Cassini - latitudinal
%      Time: 2004-06-11 19:33:04 (TDB)
%        Position:      2266.580876        -2.712515        -0.042296
%        Velocity:        -1.730462         0.002416        -0.000706
%      Time: 2004-06-11 19:41:04 (TDB)
%        Position:      3265.953141        -1.691957        -0.256502
%        Velocity:         4.727809         0.001150        -0.000104
%      Time: 2004-06-11 19:49:04 (TDB)
%        Position:      5897.757219        -1.417457        -0.267919
%        Velocity:         5.901077         0.000225         0.000006
%
%      Verification: Phoebe as seen by Cassini - rectangular
%      Time: 2004-06-11 19:33:04 (TDB)
%        Position:     -2059.271283      -942.128329       -95.837672
%        Velocity:         3.910113        -4.228139        -1.526561
%      Time: 2004-06-11 19:41:04 (TDB)
%        Position:      -381.823839     -3135.942729      -828.567239
%        Velocity:         3.064009        -4.893554        -1.526501
%      Time: 2004-06-11 19:49:04 (TDB)
%        Position:       868.679897     -5620.615726     -1561.286108
%        Velocity:         2.133569        -5.438272        -1.526496
%
%
%   2) Transform a given state from cylindrical to planetographic
%      coordinates with respect to Earth.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model and orientation data for the Earth.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function xfmsta_ex2()
%
%         %   Initialize the cylindrical state.
%         %
%         state_cyl = [1 0.5 0.5 0.2 0.1 -0.2]';
%
%         %
%         %   Load kernels.
%         %
%         cspice_furnsh( 'pck00010.tpc' );
%
%         %
%         %   Transform the state from cylindrical to planetographic.
%         %   Note that since one of the coordinate systems is
%         %   planetographic, the body name must be input.
%         %
%         state_plan = cspice_xfmsta ( state_cyl, 'cylindrical',  ...
%                                      'planetographic', 'earth' );
%
%         %
%         %   Transform the state back to cylindrical from
%         %   planetographic. The result should be very close to 'state_cyl'.
%         %
%         state_cyl2 = cspice_xfmsta ( state_plan, 'planetographic',...
%                                      'cylindrical', 'earth' );
%
%         %
%         %   Report the results.
%         %
%         fprintf('\nCylindrical State\n')
%         fprintf('  Position: %9.3f %9.3f %9.3f\n', state_cyl(1:3))
%         fprintf('  Velocity: %9.3f %9.3f %9.3f\n', state_cyl(4:6))
%
%         fprintf('\nPlanetographic State\n')
%         fprintf('  Position: %9.3f %9.3f %9.3f\n', state_plan(1:3))
%         fprintf('  Velocity: %9.3f %9.3f %9.3f\n', state_plan(4:6))
%
%         fprintf('\nVerification: Cylindrical State\n')
%         fprintf('  Position: %9.3f %9.3f %9.3f\n', state_cyl2(1:3))
%         fprintf('  Velocity: %9.3f %9.3f %9.3f\n', state_cyl2(4:6))
%
%         %
%         %   Unload the kernels.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Cylindrical State
%        Position:     1.000     0.500     0.500
%        Velocity:     0.200     0.100    -0.200
%
%      Planetographic State
%        Position:     0.500     1.548 -6356.240
%        Velocity:     0.100    -0.005    -0.195
%
%      Verification: Cylindrical State
%        Position:     1.000     0.500     0.500
%        Velocity:     0.200     0.100    -0.200
%
%
%-Particulars
%
%   Input Order
%   -----------
%
%   The input and output states will be structured by the
%   following descriptions.
%
%   For rectangular coordinates, the state vector is the following
%   in which `x', `y', and `z' are the rectangular position components and
%   `dx', `dy', and `dz' are the time derivatives of each position
%   component.
%
%           istate = ( x, y, z, dx, dy, dz );
%
%   For cylindrical coordinates, the state vector is the following
%   in which `r' is the radius, `long' is the longitudes, `z' is the
%   height, and `dr', `dlong', and `dz' are the time derivatives of each
%   position component.
%
%           istate = ( r, long, z, dr, dlong, dz );
%
%   For latitudinal coordinates, the state vector is the following
%   in which `r' is the radius, `long' is the longitude, `lat' is the
%   latitude, and `dr', `dlong', and `dlat' are the time derivatives of
%   each position component.
%
%           istate = ( r, long, lat, dr, dlong, dlat );
%
%   For spherical coordinates, the state vector is the following in
%   which `r' is the radius, `colat' is the colatitude, `long' is the
%   longitude, and `dr', `dcolat', and `dlong' are the time derivatives of
%   each position component.
%
%           istate = ( r, colat, long, dr, dcolat, dlong );
%
%   For geodetic coordinates, the state vector is the following in
%   which `long' is the longitude, `lat' is the latitude, `alt' is the
%   altitude, and `dlong', `dlat', and `dalt' are the time derivatives of
%   each position component.
%
%           istate = ( long, lat, alt, dlong, dlat, dalt );
%
%   For planetographic coordinates, the state vector is the
%   following in which `long' is the longitude, `lat' is the latitude,
%   `alt' is the altitude, and `dlong', `dlat', and `dalt' are the time
%   derivatives of each position component.
%
%           istate = ( long, lat, alt, dlong, dlat, dalt );
%
%
%   Input Boundaries
%   ----------------
%
%   There are intervals the input angles must fall within if
%   the input coordinate system is not rectangular. These
%   intervals are provided below.
%
%      Input variable    Input meaning   Input interval [rad]
%      --------------    -------------   ------------------------
%          long           Longitude        0     <= long  <  2*pi
%          lat            Latitude        -pi/2  <= lat   <= pi/2
%          colat          Colatitude       0     <= colat <= pi
%
%-Exceptions
%
%   1)  If either the input or output coordinate system is not
%       recognized, the error SPICE(COORDSYSNOTREC) is signaled by a
%       routine in the call tree of this routine.
%
%   2)  If the input body name cannot be converted to a NAIF ID
%       (applies to geodetic and planetographic coordinate systems),
%       the error SPICE(IDCODENOTFOUND) is signaled by a routine in
%       the call tree of this routine.
%
%   3)  If the input state `istate' is not valid, meaning the position
%       but not the velocity is along the Z-axis, the error
%       SPICE(INVALIDSTATE) is signaled by a routine in the call tree
%       of this routine.
%
%       Note: If both the input position and velocity are along
%       the Z-axis and the output coordinate system is not
%       rectangular, the velocity can still be calculated even
%       though the Jacobian is undefined. This case will not
%       signal an error. An example of the input position and
%       velocity along the Z-axis is below.
%
%                     Term    Value
%                     -----   ------
%                       x       0
%                       y       0
%                       z       z
%                     dx/dt     0
%                     dy/dt     0
%                     dz/dt   dz_dt
%
%   4)  If either the input or output coordinate system is geodetic or
%       planetographic and at least one of the body's radii is less
%       than or equal to zero, the error SPICE(INVALIDRADIUS) is
%       signaled by a routine in the call tree of this routine.
%
%   5)  If either the input or output coordinate system is geodetic or
%       planetographic and the difference of the equatorial and polar
%       radii divided by the equatorial radius would produce numeric
%       overflow, the error SPICE(INVALIDRADIUS) is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If the product of the Jacobian and velocity components may
%       lead to numeric overflow, the error SPICE(NUMERICOVERFLOW) is
%       signaled by a routine in the call tree of this routine.
%
%   7)  If radii for `body' are not found in the kernel pool, an error
%       is signaled by a routine in the call tree of this routine.
%
%   8)  If the size of the `body' radii kernel variable is not three,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   9)  If any of the three `body' radii is less-than or equal to zero,
%       an error is signaled by a routine in the call tree of this
%       routine.
%
%   10) If body's equatorial radii are not equal and either the input
%       or output coordinate system is geodetic or planetographic, the
%       error SPICE(NOTSUPPORTED) is signaled by a routine in the call
%       tree of this routine.
%
%   11) If any of the input arguments, `istate', `icosys', `ocosys' or
%       `body', is undefined, an error is signaled by the Matlab error
%       handling system.
%
%   12) If any of the input arguments, `istate', `icosys', `ocosys' or
%       `body', is not of the expected type, or it does not have the
%       expected dimensions and size, an error is signaled by the Mice
%       interface.
%
%-Files
%
%   SPK, PCK, CK, and FK kernels may be required.
%
%   If the input or output coordinate systems are either geodetic or
%   planetographic, a PCK providing the radii of the body
%   name `body' must be loaded via cspice_furnsh.
%
%   Kernel data are normally loaded once per program run, NOT every
%   time this routine is called.
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
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 01-NOV-2021 (EDW) (JDR)
%
%       Changed argument names "input_state", "input_coord_sys",
%       "output_coord_sys" and "output_state" to "istate", "icosys", "ocosys"
%       and "ostate" for consistency with other routines.
%
%       Edited the header to comply with NAIF standard.
%
%       Updated Examples' kernels set to use PDS archived data.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections. Edited
%       the header to comply with NAIF standard.
%
%       Added square brackets to output argument in function declaration.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 28-AUG-2012 (SCK) (EDW)
%
%-Index_Entries
%
%   state transformation between coordinate systems
%   convert state
%
%-&

function [ostate] = cspice_xfmsta( istate, icosys, ocosys, body )

   switch nargin
      case 4

         istate = zzmice_dp (istate);
         icosys = zzmice_str(icosys);
         ocosys = zzmice_str(ocosys);
         body   = zzmice_str(body);

      otherwise

         error ( ['Usage: [_ostate(6)_] = cspice_xfmsta( `_istate(6)_`, ' ...
                  '`icosys`, `ocosys`, `body`)'] )

   end

   %
   % Call the MEX library.
   %
   try
      [ostate] = mice('xfmsta_c', istate, icosys, ocosys, body);
   catch spiceerr
      rethrow(spiceerr)
   end
