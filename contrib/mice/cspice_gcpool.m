%-Abstract
%
%   CSPICE_GCPOOL returns the value of a string kernel variable
%   (scalar or array) from the kernel pool.
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
%      name     name of a pool variable associated to string values.
%
%               [1,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,1] = size(name); cell = class(name)
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
%      [cvals, found] = cspice_gcpool( name, start, room )
%
%   returns:
%
%      cvals    the values assigned to `name' beginning at index `start'.
%
%               [n,c2] = size(cvals); char = class(cvals)
%
%               `cvals' returns empty if the variable `name' does not exist in
%               the kernel pool.
%
%               `cvals' has a size of `room' or less (n <= room).
%
%      found    the flag indicating true if `name' exists in the kernel pool
%               and has character type, false if it is not.
%
%               [1,1] = size(found); logical = class(found)
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
%   1) The following code example demonstrates how the data stored
%      in a kernel pool variable can be retrieved in pieces.
%
%      Use the kernel shown below to load the kernel pool with the
%      variables used within the example.
%
%
%         KPL/MK
%
%         File name: gcpool_ex1.tm
%
%         This kernel is intended to support operation of SPICE
%         example programs.
%
%         \begindata
%
%            CTEST_VAL = ('LARRY', 'MOE', 'CURLY' )
%
%            ITEST_VAL = ( 3141, 186, 282 )
%
%            DTEST_VAL = ( 3.1415, 186. , 282.397 )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function gcpool_ex1()
%
%         %
%         % Load the test data.
%         %
%         cspice_furnsh( 'gcpool_ex1.tm' )
%
%         %
%         % Retrieve up-to 'ROOM' character entries for
%         % kernel pool variable named 'CTEST_VAL' to
%         % the array named 'cvals'. The first index to return,
%         % 'START', has value 1 (this returns all strings).
%         %
%         VAR    = 'CTEST_VAL';
%         ROOM   = 25;
%         START  = 1;
%
%         %
%         % cspice_gcpool returns an empty array if the variable
%         % does not exist in the kernel pool.
%         %
%         [cvals, found] = cspice_gcpool( VAR, START, ROOM );
%
%         if ( found )
%
%            txt = sprintf( 'Found %s in the kernel pool', VAR );
%            disp(txt)
%
%            n_elements = size( cvals, 1 );
%
%            %
%            % Retrieve the number of elements returned in 'cvals' from the
%            % first element returned from "size".
%            %
%            for n=1:n_elements
%               fprintf( '   Element %d of %s: ``%s``\n', ...
%                         n, VAR, cvals(n,:) );
%            end
%
%         else
%
%            txt = sprintf( 'Failed to find %s in the kernel pool', VAR );
%            disp(txt)
%
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
%      Found CTEST_VAL in the kernel pool
%         Element 1 of CTEST_VAL: ``LARRY``
%         Element 2 of CTEST_VAL: ``MOE  ``
%         Element 3 of CTEST_VAL: ``CURLY``
%
%
%-Particulars
%
%   This routine provides the user interface to retrieving
%   character data stored in the kernel pool. This interface
%   allows you to retrieve the data associated with a variable
%   in multiple accesses. Under some circumstances this alleviates
%   the problem of having to know in advance the maximum amount
%   of space needed to accommodate all kernel variables.
%
%   However, this method of access does come with a price. It is
%   always more efficient to retrieve all of the data associated
%   with a kernel pool data in one call than it is to retrieve
%   it in sections.
%
%   See also the routines cspice_gdpool and cspice_gipool.
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
%       Edited the header to comply with NAIF standard. Added example's input
%       data and problem statement.
%
%       Added -Parameters, -Particulars, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.1.1, 01-JUN-2016 (EDW)
%
%       Corrected typo in example code. Ouput loop showed use of "i" rather
%       than "n."
%
%   -Mice Version 1.1.0, 12-MAR-2012 (EDW) (SCK)
%
%       "logical" call replaced with "zzmice_logical."
%
%       -I/O descriptions edits to conform to Mice documentation format.
%
%   -Mice Version 1.0.0, 06-MAR-2007 (EDW)
%
%-Index_Entries
%
%   RETURN the character value of a pooled kernel variable
%   RETURN the string value of a pooled kernel variable
%
%-&

function [cvals, found] = cspice_gcpool( name, start, room )

   switch nargin
      case 3

         name  = zzmice_str(name);
         start = zzmice_int(start);
         room  = zzmice_int(room);

      otherwise

         error ( ['Usage: [cvals(), found] = ' ...
                  'cspice_gcpool( `name`, start, room)' ] )

   end

   %
   % Call the MEX library.
   %
   try
      [cvals, found] = mice( 'gcpool_c', name, start, room );

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end


