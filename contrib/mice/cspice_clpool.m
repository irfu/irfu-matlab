%-Abstract
%
%   CSPICE_CLPOOL removes all kernel variables from the kernel pool. Watches
%   on kernel variables are retained.
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
%   The call:
%
%      cspice_clpool
%
%   deletes all variable assignments loaded into the kernel
%   pool.
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
%   1) This code example demonstrates how to assing values to kernel
%      pool variables, how to check for the existence of kernel pool
%      variables and how to clear the kernel pool, i.e. how to delete
%      all variable assignments loaded into the kernel pool.
%
%      Place a value into the kernel pool and check for the variable
%      to which the value has been assigned. Clear the kernel pool
%      and check for that variable again.
%
%      Example code begins here.
%
%
%      function clpool_ex1()
%
%         %
%         % Place a value into the kernel pool. Recall
%         % the routines for direct insertion
%         % of pool assignments have arrays for input,
%         % but in MATLAB a scalar is a 1x1 array.
%         %
%         cspice_pdpool( 'TEST_VAR', -666. )
%
%         %
%         % Check for the variable assignment to TEST_VAR.
%         % cspice_gdpool returns an empty array if the variable
%         % does not exist in the kernel pool.
%         %
%         dvals = cspice_gdpool( 'TEST_VAR', 0, 1 );
%
%         disp( 'First call to cspice_gdpool:' )
%         if ( ~isempty(dvals) )
%            disp( sprintf( '   TEST_VAR value: %f', dvals ) )
%         end
%
%         %
%         % Now clear the kernel pool.
%         %
%         cspice_clpool
%
%         %
%         % Again, check for the TEST_VAR assignment.
%         %
%         dvals = cspice_gdpool( 'TEST_VAR', 0, 1 );
%
%         disp( 'Second call to cspice_gdpool:' )
%         if ( isempty(dvals)  )
%            disp( '   TEST_VAR not in kernel pool' )
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      First call to cspice_gdpool:
%         TEST_VAR value: -666.000000
%      Second call to cspice_gdpool:
%         TEST_VAR not in kernel pool
%
%
%-Particulars
%
%   cspice_clpool clears the pool of kernel variables maintained by
%   the kernel POOL subsystem. All the variables in the pool are
%   deleted. However, all watcher information is retained.
%
%   Each watched variable will be regarded as having been updated.
%   Any agent associated with that variable will have a notice
%   posted for it indicating that its watched variable has been
%   updated.
%
%   Note, cspice_clpool deletes ALL pool assignments, including those
%   from the cspice_pipool, cspice_pdpool, and cspice_pcpool
%   set. Use cspice_unload to remove the assignments loaded from a
%   particular kernel or cspice_kclear to completely clear the kernel
%   pool.
%
%-Exceptions
%
%   1)  All known agents (those established through the SPICELIB
%       routine SWPOOL) will be "notified" that their watched
%       variables have been updated whenever cspice_clpool is called.
%
%-Files
%
%   None.
%
%-Restrictions
%
%   1)  This routine should not be used to unload kernels that
%       have been loaded via cspice_furnsh.
%
%-Required_Reading
%
%   KERNEL.REQ
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   J. Diaz del Rio     (ODC Space)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 21-JUL-2020 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard.
%       Reformatted example's output and added problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   CLEAR the pool of kernel variables
%
%-&

function cspice_clpool

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: cspice_clpool' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('clpool_c');
   catch spiceerr
      rethrow(spiceerr)
   end

