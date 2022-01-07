%-Abstract
%
%   CSPICE_KCLEAR clears the KEEPER system: unload all kernels, clears
%   the kernel pool, and re-initialize the system.
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
%   The call:
%
%      cspice_kclear
%
%      Re-initialize the KEEPER system.
%
%-Parameters
%
%   None.
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%    platforms as the results depend on the SPICE kernels used as input
%    and the machine specific arithmetic implementation.
%
%    1) Load a meta-kernel containing three kernels, and count the
%       number of files in the kernel pool before and after calling
%       cspice_kclear.
%
%       Use the meta-kernel shown below to load the required SPICE
%       kernels.
%
%
%          KPL/MK
%
%          File name: kclear_ex1.tm
%
%          This meta-kernel is intended to support operation of SPICE
%          example programs. The kernels shown here should not be
%          assumed to contain adequate or correct versions of data
%          required by SPICE-based user applications.
%
%          In order for an application to use this meta-kernel, the
%          kernels referenced here must be present in the user's
%          current working directory.
%
%          The names and contents of the kernels referenced
%          by this meta-kernel are as follows:
%
%             File name                     Contents
%             ---------                     --------
%             de421.bsp                     Planetary ephemeris
%             pck00008.tpc                  Planet orientation and
%                                           radii
%             naif0009.tls                  Leapseconds
%
%
%          \begindata
%
%             KERNELS_TO_LOAD = ( 'de421.bsp',
%                                 'pck00008.tpc',
%                                 'naif0009.tls'  )
%
%          \begintext
%
%          End of meta-kernel
%
%
%       Example code begins here.
%
%
%       function kclear_ex1()
%
%          %
%          % Load the standard meta kernel, retrieve the number of
%          % loaded kernels.
%          %
%          cspice_furnsh( 'kclear_ex1.tm' )
%
%          n   = cspice_ktotal( 'ALL' );
%          txt = sprintf(['Count of loaded kernels before ', ...
%                         'cspice_kclear call: %d'], n     );
%          disp( txt )
%
%          %
%          % Clear the KEEPER system, retrieve the number of loaded
%          % after the clear.
%          %
%          cspice_kclear
%
%          n   = cspice_ktotal( 'ALL' );
%          txt = sprintf(['Count of loaded kernels after ', ...
%                          'cspice_kclear call:  %d'], n   );
%          disp( txt )
%
%
%       When this program was executed on a Mac/Intel/Octave6.x/64-bit
%       platform, the output was:
%
%
%       Count of loaded kernels before cspice_kclear call: 4
%       Count of loaded kernels after cspice_kclear call:  0
%
%
%-Particulars
%
%   This routine allows you re-initialize the KEEPER system with
%   a single call. The KEEPER system is the kernel management system
%   underlying the set of Mice APIs
%
%      cspice_furnsh
%      cspice_ktotal
%      cspice_kdata
%      cspice_kinfo
%      cspice_kclear
%      cspice_unload
%
%   This routine unloads all kernels from their kernel-type-specific
%   kernel management subsystems (SPKBSR, CKBSR, etc.), clears the
%   kernel pool, clears KEEPER's internal file database, and re-sets
%   the watch status for the kernel variables used to load kernels
%   via meta-kernels.
%
%   This capability, though implemented in Fortran, is particularly
%   relevant to SPICE implementations such as Mice, for which the
%   state of the KEEPER system persists after any Mice-based MATLAB
%   script is run. Successive runs of Mice-based scripts may perform
%   in unexpected ways when scripts access data loaded during runs of
%   previous scripts.
%
%   Cleaning up after such programs using explicit unload_c commands is
%   tedious and error-prone. One call to this routine sets the
%   KEEPER system to its initial state, preventing unintentional
%   interaction between scripts via KEEPER's state.
%
%-Exceptions
%
%   1)  If an error occurs when setting a kernel pool watch or
%       checking watched variables, the error is signaled by a routine
%       in the call tree of this routine.
%
%-Files
%
%   See -Particulars.
%
%-Restrictions
%
%   1)  Calling this routine will wipe out any kernel pool data
%       inserted via the Mice API routines to put data into the
%       kernel pool (cspice_pcpool, cspice_pdpool and cspice_pipool).
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
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 1.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Edited the -Examples section to comply with NAIF standard. Added
%       example's problem statement and meta-kernel. Merged the existing
%       code fragments into a complete example.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.0.0, 01-DEC-2006 (EDW)
%
%-Index_Entries
%
%   Re-initialize the keeper system
%   Clear the keeper system
%   Unload all kernels
%
%-&

function cspice_kclear

   switch nargin
      case 0
         ;
      otherwise

         error ( 'Usage: cspice_kclear' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('kclear_c');
   catch spiceerr
      rethrow(spiceerr)
   end
