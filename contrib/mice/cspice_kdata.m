%-Abstract
%
%   CSPICE_KDATA returns data for the nth kernel that is among a list of
%   specified kernel types.
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
%      which    the number of the kernel to fetch (matching the type
%               specified by `kind') from the list of kernels that have been
%               loaded through the routine cspice_furnsh but that have not
%               been unloaded through the routine cspice_unload.
%
%               [1,1] = size(which); int32 = class(which)
%
%               The range of `which' is 1 to `count', where `count' is the
%               number of kernels loaded via cspice_furnsh of type `kind'.
%               This count may be obtained by calling cspice_ktotal. See the
%               -Examples section for an illustrative example.
%
%      kind     a list of types of kernels to be considered when fetching
%               kernels from the list of loaded kernels.
%
%               [1,c1] = size(kind); char = class(kind)
%
%                  or
%
%               [1,1] = size(kind); cell = class(kind)
%
%               `kind' should consist of words from list of kernel types
%               given below.
%
%                  SPK  --- All SPK files are counted in the total.
%                  CK   --- All CK files are counted in the total.
%                  DSK  --- All DSK files are counted in the total.
%                  PCK  --- All binary PCK files are counted in the
%                           total.
%                  EK   --- All EK files are counted in the total.
%                  TEXT --- All text kernels that are not meta-text
%                           kernels are included in the total.
%                  META --- All meta-text kernels are counted in the
%                           total.
%                  ALL  --- Every type of kernel is counted in the
%                           total.
%
%               `kind' is case insensitive. If a word appears in `kind'
%               that is not one of those listed above, it is ignored.
%
%               When `kind' consists of multiple words, the words must
%               be separated by blanks. Examples of valid lists are the
%               strings
%
%                  'SPK CK TEXT'
%                  'SPK CK text'
%                  'PCK DSK'
%                  'CK'
%                  'ALL'
%
%               See the routine cspice_ktotal for examples of the use of
%               `kind'.
%
%   the call:
%
%      [file, filtyp, srcfil, handle, found] = cspice_kdata( which, kind )
%
%   returns:
%
%      file     the name of the file having index `which' in the sequence of
%               files of type `kind' that is currently loaded via
%               cspice_furnsh.
%
%               [1,c2] = size(file); char = class(file)
%
%               `file' will be empty if there is not such kernel loaded.
%
%      filtyp   the type of the kernel specified by `file'.
%
%               [1,c3] = size(filtyp); char = class(filtyp)
%
%               `file' will be empty if there is no file matching the
%               specification of `which' and `kind'.
%
%      srcfil   the name of the source file that was used to specify `file'
%               as one to load.
%
%               [1,c4] = size(srcfil); char = class(srcfil)
%
%               If `file' was loaded directly via a call to cspice_furnsh,
%               `srcfil' will be empty. If there is no file matching the
%               specification of `which' and `kind', `srcfil' will be empty.
%
%      handle   the handle attached to `file' if it is a binary kernel.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               If `file' is a text kernel or meta-text kernel `handle'
%               will be zero. If there is no file matching the specification
%               of `which' and `kind', `handle' will be set to zero.
%
%      found    returned true if a `file' matching the specification of
%               `which' and `kind' exists.
%
%               [1,1] = size(found); logical = class(found)
%
%               If there is no such file, `found' will be set to false.
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
%   1) Load a meta-kernel with a PCK, an LSK and an SPK and loop over
%      the loaded kernels, outputting file information for each of
%      them.
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: kdata_ex1.tm
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
%            pck00009.tpc                  Planet orientation and
%                                          radii
%            naif0009.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'pck00009.tpc',
%                                'naif0009.tls'  )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Example code begins here.
%
%
%      function kdata_ex1()
%
%         %
%         % Load several kernel files.
%         %
%         cspice_furnsh( 'kdata_ex1.tm' )
%
%         %
%         % Count the number of loaded kernel files.
%         %
%         count = cspice_ktotal( 'ALL' );
%
%         %
%         % Loop over the count, outputting file information as we loop.
%         % The loop tells us all files loaded via cspice_furnsh, their
%         % type, and how they were loaded.
%         %
%         for i = 1:count+1
%
%            [ file, type, srcfil, handle, found ] = ...
%                                             cspice_kdata( i, 'ALL');
%
%            if ( found )
%               fprintf( 'Index : %d\n',   i     );
%               fprintf( 'File  : %s\n',   file  );
%               fprintf( 'Type  : %s\n',   type  );
%               fprintf( 'Source: %s\n\n', srcfil);
%
%            else
%
%               fprintf( 'No kernel found with index: %d\n', i );
%
%            end
%
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
%      Index : 1
%      File  : kdata_ex1.tm
%      Type  : META
%      Source:
%
%      Index : 2
%      File  : de421.bsp
%      Type  : SPK
%      Source: kdata_ex1.tm
%
%      Index : 3
%      File  : pck00009.tpc
%      Type  : TEXT
%      Source: kdata_ex1.tm
%
%      Index : 4
%      File  : naif0009.tls
%      Type  : TEXT
%      Source: kdata_ex1.tm
%
%      No kernel found with index: 5
%
%
%-Particulars
%
%   This routine allows you to determine which kernels have been
%   loaded via cspice_furnsh and to obtain information sufficient to directly
%   query those files.
%
%-Exceptions
%
%   1)  If a file is not loaded matching the specification of `which' and
%       `kind', `found' will be false, `file', `filtyp', and `srcfil' will be
%       empty and `handle' will be set to zero.
%
%   2)  If any of the input arguments, `which' or `kind', is
%       undefined, an error is signaled by the Matlab error handling
%       system.
%
%   3)  If any of the input arguments, `which' or `kind', is not of
%       the expected type, or it does not have the expected dimensions
%       and size, an error is signaled by the Mice interface.
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
%   KERNEL.REQ
%   MICE.REQ
%
%-Literature_References
%
%   None.
%
%-Author_and_Institution
%
%   N.J. Bachman        (JPL)
%   J. Diaz del Rio     (ODC Space)
%   S.C. Krening        (JPL)
%   E.D. Wright         (JPL)
%
%-Version
%
%   -Mice Version 2.1.0, 13-AUG-2021 (EDW) (JDR)
%
%       Changed output argument name "source" to "srcfil" for consistency
%       with other routines.
%
%       Edited the header to comply with NAIF standard. Added example's
%       problem statement.
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
%       Updated -I/O description of input argument "kind", to illustrate
%       use of multi-word lists. Added kernel.req and removed dsk.req to the
%       list of required readings.
%
%   -Mice Version 2.0.0, 20-JAN-2016 (EDW) (NJB)
%
%       Corrected "Usage" string to include 'found'.
%
%       Header update to expand argument descriptions and
%       reflect support for use of DSKs
%
%   -Mice Version 1.2.0, 12-MAR-2012 (EDW) (SCK)
%
%       "logical" call replaced with "zzmice_logical."
%
%       -I/O descriptions edits to parallel the Icy version.
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%       Edits to Example section, proper description of "standard.tm"
%       meta kernel.
%
%   -Mice Version 1.0.1, 06-MAY-2009 (EDW)
%
%       Added mice.req reference to the Required Reading section.
%
%   -Mice Version 1.0.0, 30-MAR-2007 (EDW)
%
%-Index_Entries
%
%   Retrieve information on loaded SPICE kernels
%
%-&

function [ file, filtyp, srcfil, handle, found ] = cspice_kdata( which, kind )

   switch nargin
      case 2

         which = zzmice_int(which);
         kind  = zzmice_str(kind);

      otherwise

         error( [ 'Usage: [ `file`, `filtyp`, `srcfil`, handle, found ] = ' ...
                                         'cspice_kdata( which, `kind` )']  )

   end

   %
   % Call the MEX library.
   %
   try
      [ file, filtyp, srcfil, handle, found ]  = mice('kdata_c', which, kind);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end
