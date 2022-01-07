%-Abstract
%
%   CSPICE_KINFO returns information about a loaded kernel
%   specified by name.
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
%      file     the name of a kernel file for which descriptive
%               information is desired.
%
%               [1,c1] = size(file); char = class(file)
%
%                  or
%
%               [1,1] = size(file); cell = class(file)
%
%   the call:
%
%      [filtyp, srcfil, handle, found] = cspice_kinfo( file )
%
%   returns:
%
%      filtyp   the type name of the kernel specified by `file'.
%               `filtyp' will be empty if file is not on the list of kernels
%               loaded via cspice_furnsh.
%
%               [1,c2] = size(file); char = class(file)
%
%      srcfil   the name of the source file used to
%               specify `file' as one to load. If `file' was loaded
%               directly via a call to cspice_furnsh, `srcfil' will be empty.
%               If file is not on the list of kernels loaded via
%               cspice_furnsh, `srcfil' will be empty.
%
%               [1,c3] = size(srcfil); char = class(srcfil)
%
%      handle   the integer handle attached to 'file' if it is a binary
%               kernel. If file is a text kernel or meta-text kernel
%               handle will be zero. If file is not on the list of
%               kernels loaded via cspice_furnsh, 'handle' has value zero.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%      found    returns true if the specified file exists.
%               If there is no such file, 'found' will be set to
%               false.
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
%   1) Load a meta-kernel listing a path to an SPK file, and verify
%      that the kernel system loaded the SPK file of interest.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: kinfo_ex1.tm
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
%            de421.bsp                     Planetary ephemeris
%
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function kinfo_ex1()
%
%         %
%         % Load a meta kernel listing a path to an SPK file.
%         %
%         cspice_kclear
%         cspice_furnsh( 'kinfo_ex1.tm' )
%
%         %
%         % Use cspice_kinfo to ensure the kernel system loaded
%         % the SPK file of interest.
%         %
%         file = 'de421.bsp';
%
%         [filtyp, srcfil, handle, found ] = cspice_kinfo( file );
%
%         %
%         % Take appropriate action depending on the returned
%         % state of found. If found has value false, then
%         % `file' is not loaded.
%         %
%         if ( found )
%            disp( [ 'File type: ' filtyp ] )
%            disp( [ 'Source   : ' srcfil ] )
%         else
%            disp( [ 'Kernel not loaded: ' file ] )
%         end
%
%         %
%         % It's always good form to unload kernels after use,
%         % particularly in Mice due to data persistence.
%         %
%         cspice_kclear
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      File type: SPK
%      Source   : kinfo_ex1.tm
%
%
%-Particulars
%
%   This routine allows you to request information directly
%   for a specific SPICE kernel.
%
%-Exceptions
%
%   1)  If the specified file is not on the list of files that
%       are currently loaded via the interface cspice_furnsh, `found'
%       will be false, `handle' will be set to zero and `filtyp'
%       and `srcfil' will be set to blanks.
%
%   2)  If the input argument `file' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `file' is not of the expected type, or
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
%   DSK.REQ
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
%   -Mice Version 1.3.0, 10-AUG-2021 (EDW) (JDR)
%
%       Changed the output argument name "source" to "srcfil" for
%       consistency with other routines.
%
%       Edited -Examples section to comply with NAIF standard. Added
%       example's problem statement.
%
%       Added -Parameters, -Exceptions, -Files, -Restrictions,
%       -Literature_References and -Author_and_Institution sections, and
%       completed -Particulars section.
%
%       Eliminated use of "lasterror" in rethrow.
%
%       Removed reference to the function's corresponding CSPICE header from
%       -Required_Reading section.
%
%   -Mice Version 1.2.1, 01-DEC-2014 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.2.0, 10-MAY-2011 (EDW)
%
%       "logical" call replaced with "zzmice_logical."
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 01-DEC-2006 (EDW)
%
%-Index_Entries
%
%   Fetch information about a loaded SPICE kernel
%
%-&

function [filtyp, srcfil, handle, found] = cspice_kinfo( file )

   switch nargin
      case 1

         file = zzmice_str(file);

      otherwise

         error( [ 'Usage: [ `filtyp`, `srcfil`, handle, found ] = ' ...
                                             'cspice_kinfo( `file` )']  )

   end

   %
   % Call the MEX library.
   %
   try
      [filtyp, srcfil, handle, found]  = mice('kinfo_c', file);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end

