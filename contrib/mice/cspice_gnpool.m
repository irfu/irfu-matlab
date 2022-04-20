%-Abstract
%
%   CSPICE_GNPOOL returns the names of kernel variables matching a
%   specified template.
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
%      name     the template to use when searching for variable names
%               in the kernel pool.
%
%               [1,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
%
%               The characters '*' and '%' are used for the wild card string
%               and wild card character respectively. For details of string
%               pattern matching see the header of the CSPICE routine matchi_c.
%
%      start    value for the index indicating the first component of the data
%               vector assigned to `name' for return (index 1 for all
%               elements).
%
%               [1,1] = size(start); int32 = class(start)
%
%      room     value specifying the maximum number of components that can
%               return for `name'.
%
%               [1,1] = size(room); int32 = class(room)
%
%   the call:
%
%      [cvals, found] = cspice_gnpool( name, start, room )
%
%   returns:
%
%      cvals    the values assigned to `name' beginning at index `start'.
%
%               [n,c2] = size(cvals); char = class(cvals)
%
%               `cvals' returns empty if variables described by `name' does
%               not exist in the kernel pool.
%
%               n <= room, p = length of longest string in return array
%               `cvals'.
%
%      found    the flag indicating true if variables matching the `name'
%               template exist in the kernel pool, false otherwise.
%
%               [1,1] = size(found); logical = class(found)
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
%   1) Load a PCK kernel, create a template for Jupiter kernel
%      variables, and after performing a query for them, output all the
%      variable names found in the kernel pool that match that template.
%
%      Use the PCK kernel below to load the triaxial ellipsoidal shape
%      model and orientation data for Jupiter.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function gnpool_ex1()
%
%         %
%         % Load a PCK kernel.
%         %
%         cspice_furnsh( 'pck00010.tpc' )
%
%         %
%         % A template for Jupiter kernel variables.
%         %
%         VAR = 'BODY599*';
%
%         %
%         % Query for the variable name, return all matches from
%         % index 1.
%         %
%         INDEX  = 1;
%         ROOM   = 10;
%
%         [kervar, found] = cspice_gnpool( VAR, INDEX, ROOM );
%
%         if( found )
%
%            n_elements = size(kervar, 1);
%
%            %
%            % Output the returned variable names.
%            %
%            for n=1: n_elements
%               txt = sprintf( 'Variable %d matching %s: %s', n, VAR, ...
%                              kervar(n,:));
%               disp( txt )
%            end
%
%         else
%            txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%            disp( txt )
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
%      platform, the output was:
%
%
%      Variable 1 matching BODY599*: BODY599_PM
%      Variable 2 matching BODY599*: BODY599_LONG_AXIS
%      Variable 3 matching BODY599*: BODY599_RADII
%      Variable 4 matching BODY599*: BODY599_NUT_PREC_DEC
%      Variable 5 matching BODY599*: BODY599_NUT_PREC_PM
%      Variable 6 matching BODY599*: BODY599_POLE_RA
%      Variable 7 matching BODY599*: BODY599_POLE_DEC
%      Variable 8 matching BODY599*: BODY599_NUT_PREC_RA
%
%
%   2) Obtain from the kernel pool the names of the first
%      10 variables stored. Use the * wildcard character
%      as a template to indicate a request for all kernel
%      variables.
%
%      Use the PCK kernel below to load the triaxial ellipsoidal shape
%      model and orientation data for all the Solar System planets.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function gnpool_ex2()
%
%         %
%         % Load a PCK kernel.
%         %
%         cspice_furnsh( 'pck00010.tpc' )
%         %
%         % Return all matches from 'INDEX' 1.
%         %
%         INDEX  = 1;
%         ROOM   = 10;
%         VAR    = '*';
%
%         [kervar, found] = cspice_gnpool( VAR, INDEX, ROOM );
%
%         if ( found )
%
%            n_elements = size(kervar, 1);
%
%            %
%            % Output the returned variable names.
%            %
%            for n=1: n_elements
%               txt = sprintf( 'Variable %d matching %s: %s', n, VAR, ...
%                              kervar(n,:));
%               disp( txt )
%            end
%
%         else
%            txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%            disp( txt )
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in MATLAB due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Variable 1 matching *: BODY611_LONG_AXIS
%      Variable 2 matching *: BODY4_NUT_PREC_ANGLES
%      Variable 3 matching *: BODY604_POLE_RA
%      Variable 4 matching *: BODY605_POLE_DEC
%      Variable 5 matching *: BODY399_N_GEOMAG_CTR_DIPOLE_LAT
%      Variable 6 matching *: BODY399_POLE_RA
%      Variable 7 matching *: BODY703_NUT_PREC_PM
%      Variable 8 matching *: BODY708_LONG_AXIS
%      Variable 9 matching *: BODY501_NUT_PREC_RA
%      Variable 10 matching *: BODY710_NUT_PREC_PM
%
%
%      Note, the seemingly random order of the output list reflects the
%      order used by the SPICE kernel subsystem to store/lookup the
%      variable names.
%
%-Particulars
%
%   This routine provides the user interface for retrieving the names
%   of kernel pool variables. This interface allows you to retrieve
%   the names matching a template via multiple accesses. Under some
%   circumstances this alleviates the problem of having to know in
%   advance the maximum amount of space needed to accommodate all
%   matching names.
%
%   However, this method of access does come with a price. It is
%   always more efficient to retrieve all of the data associated with
%   a kernel pool variable in one call than it is to retrieve it in
%   sections.
%
%-Exceptions
%
%   1)  If the value of `room' is less than one, the error
%       SPICE(BADARRAYSIZE) is signaled by a routine in the call tree
%       of this routine.
%
%   2)  If any of the input arguments, `name', `start' or `room', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `name', `start' or `room', is
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
%   KERNEL.REQ
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
%   -Mice Version 1.2.0, 26-NOV-2021 (EDW) (JDR)
%
%       Changed the output argument name "kvars" to "cvals" for consistency
%       with other routines. Corrected typo in usage message.
%
%       Edited the header to comply with NAIF standard. Added examples'
%       problem statement and a reference to the required PCK.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.0, 12-MAR-2012 (EDW) (SCK)
%
%       "logical" call replaced with "zzmice_logical."
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 15-DEC-2006 (EDW)
%
%-Index_Entries
%
%   return names of kernel pool variables matching a template
%
%-&

function [cvals,found] = cspice_gnpool( name, start, room )

   switch nargin
      case 3

         name  = zzmice_str(name);
         start = zzmice_int(start);
         room  = zzmice_int(room);

      otherwise

         error ( ['Usage: [cvals(), found] = ' ...
                  'cspice_gnpool( `name`, start, room )' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [cvals, found] = mice( 'gnpool_c', name, start, room );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end



