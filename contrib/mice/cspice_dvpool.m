%-Abstract
%
%   CSPICE_DVPOOL deletes a variable from the kernel pool.
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
%      name     name(s) of a pool variable(s) to delete from the kernel pool.
%
%               [n,c1] = size(name); char = class(name)
%
%                  or
%
%               [1,n] = size(name); cell = class(name)
%
%               The name and associated values are removed from the kernel
%               pool, freeing the occupied space.
%
%               If watches are set on the variable(s) designated by `name',
%               the corresponding agents are placed on the list of agents
%               to notify of a kernel variable update.
%
%   the call:
%
%      cspice_dvpool( name )
%
%   performs the delete operation.
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
%   1) Query the kernel pool for variables matching a pattern,
%      returning 10 or less matches starting from index #1. Then
%      remove those variables from the kernel pool, and check that
%      indeed they have been deleted.
%
%      Use the LSK kernel below to load the variables used within
%      the example.
%
%         naif0009.tls
%
%
%      Example code begins here.
%
%
%      function dvpool_ex1()
%
%         %
%         % Load an LSK.
%         %
%         cspice_furnsh( 'naif0009.tls' )
%
%         %
%         % A template for the leapseconds kernel variables.
%         %
%         VAR = 'DELTET*';
%
%         %
%         % Query for the variable name, return 10 or less matches from
%         % index 1.
%         %
%         INDEX  = 1;
%         ROOM   = 10;
%
%
%         txt = sprintf( 'Kernel pool state after load.' );
%         disp( txt )
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
%               txt = sprintf( 'Variable %d matching %s: %s', ...
%                                           n, VAR, kervar(n,:) );
%               disp( txt )
%            end
%
%         else
%            txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%            disp( txt )
%         end
%
%
%         %
%         % Delete the kernel pool variables returned from cspice_gnpool.
%         %
%         cspice_dvpool( kervar )
%
%         txt = sprintf( '\nKernel pool state after deletion.' );
%         disp( txt )
%
%         %
%         % Confirm the variables were deleted from the pool.
%         %
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
%               txt = sprintf( 'Variable %d matching %s: %s', ...
%                                           n, VAR, kervar(n,:) );
%               disp( txt )
%            end
%
%         else
%            txt = sprintf( ['Failed to find  ' VAR ' in the kernel pool.'] );
%            disp( txt )
%         end
%
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
%      Kernel pool state after load.
%      Variable 1 matching DELTET*: DELTET/DELTA_AT
%      Variable 2 matching DELTET*: DELTET/K
%      Variable 3 matching DELTET*: DELTET/M
%      Variable 4 matching DELTET*: DELTET/EB
%      Variable 5 matching DELTET*: DELTET/DELTA_T_A
%
%      Kernel pool state after deletion.
%      Failed to find  DELTET* in the kernel pool.
%
%
%-Particulars
%
%   This routine enables users to programmatically remove variables
%   from the kernel pool, as opposed to having to clear the pool and
%   reload it.
%
%   Note that it is not necessary to remove kernel variables in order
%   to simply update them; this routine should be used only when
%   variables are to be removed.
%
%-Exceptions
%
%   1)  If the specified variable is not present in the kernel pool,
%       this routine simply returns. No error is signaled.
%
%   2)  If the input argument `name' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `name' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
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
%   -Mice Version 1.1.0, 24-AUG-2021 (EDW) (JDR)
%
%       Fixed bug: the return value from the zzmice_str call was assigned
%       to the wrong variable name, causing input cell strings to break the
%       interface.
%
%       Edited the header to comply with NAIF standard. Added
%       example's meta-kernel.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 12-MAR-2012 (EDW) (SCK)
%
%-Index_Entries
%
%   delete a kernel pool variable
%
%-&

function cspice_dvpool(name)

   switch nargin
      case 1

         name = zzmice_str(name);

      otherwise

         error ( 'Usage: cspice_dvpool(_`name`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('dvpool_c', name);
   catch spiceerr
      rethrow(spiceerr)
   end


