%-Abstract
%
%   CSPICE_UNLOAD unloads a SPICE kernel file (of any type)
%   from MATLAB.
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
%      file     the name of the file(s) to unload.
%
%               [n,c1] = size(file); char = class(file)
%
%                  or
%
%               [1,n] = size(file); cell = class(file)
%
%               This file should be one loaded through the interface
%               cspice_furnsh. If the file is not on the list of loaded
%               kernels no action is taken.
%
%               Note that if `file' is a meta-text kernel, all of
%               the files loaded as a result of loading the meta-text
%               kernel will be unloaded.
%
%   the call:
%
%      cspice_unload( file )
%
%   returns:
%
%      None.
%
%      It removes the file and all associated data from the kernel
%      sub-system. If `file' is a meta-text kernel, the sub-system
%      unloads all files listed in the kernel.
%
%      Note: a cspice_unload call deletes ALL kernel variables except
%      those loaded into the kernel pool via a cspice_furnsh kernel
%      load  call, i.e. cspice_unload erases kernel variables placed
%      in the pool by the pool functions: cspice_pipool, cspice_pdpool,
%      and cspice_pcpool.
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
%   1) Load a PCK in the kernel pool and look up the kernel variable
%      BODY399_RADII, then unload the kernel and look the variable up
%      again. In the later case, the POOL subsystem shall throw
%      an error indicating that the kernel POOL variable does
%      not exist.
%
%      Use the PCK kernel below to load the required triaxial
%      ellipsoidal shape model for the Earth, which uses the kernel
%      variable BODY399_RADII to store the Earth's radii data.
%
%         pck00010.tpc
%
%
%      Example code begins here.
%
%
%      function unload_ex1()
%
%         %
%         %  Load a PCK kernel.
%         %
%         cspice_furnsh( 'pck00010.tpc' )
%
%         %
%         % When the kernel variable
%         %
%         %    BODY399_RADII
%         %
%         % is present in the kernel pool---normally because a PCK
%         % defining this variable has been loaded (as is the case
%         % here)---the call
%         %
%         disp( 'Calling cspice_bodvrd after loading the PCK:' )
%         try
%            values = cspice_bodvrd( 'EARTH', 'RADII', 3);
%            disp('   Expected result, found kernel data')
%         catch
%            disp('   ERROR: Unexpected result, no kernel data found')
%         end
%
%         %
%         %  Now unload the kernel and try again.
%         %
%         cspice_unload( 'pck00010.tpc' )
%
%         disp( 'Calling cspice_bodvrd after unloading the PCK:' )
%         try
%            values = cspice_bodvrd( 'EARTH', 'RADII', 3);
%            disp('   ERROR: Unexpected result, found kernel data')
%         catch
%            disp('   Expected result, no kernel data found')
%         end
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Calling cspice_bodvrd after loading the PCK:
%         Expected result, found kernel data
%      Calling cspice_bodvrd after unloading the PCK:
%         Expected result, no kernel data found
%
%
%   2) Load a meta-kernel with a PCK, an LSK and an SPK, and
%      separately a text kernel and a binary PCK. Loop over the
%      loaded kernels, outputting file information for each of
%      them.
%
%      Then unload the text kernels, check that they have been
%      unloaded, and finally unload the meta-kernel.
%
%
%      Use the meta-kernel shown below to load the required SPICE
%      kernels.
%
%
%         KPL/MK
%
%         File name: unload_ex2.tm
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
%            naif0012.tls                  Leapseconds
%
%         \begindata
%
%            KERNELS_TO_LOAD = ( 'de421.bsp',
%                                'naif0012.tls',
%                                'pck00009.tpc' )
%
%         \begintext
%
%         End of meta-kernel
%
%
%      Use the PCK kernel below as the binary PCK required for the
%      example.
%
%         earth_latest_high_prec.bpc
%
%
%      Use the FK kernel below as the text kernel required for the
%      example.
%
%         RSSD0002.TF
%
%
%      Example code begins here.
%
%
%      function unload_ex2()
%
%         %
%         % Load several kernel files.
%         %
%         cspice_furnsh( 'unload_ex2.tm' );
%         cspice_furnsh( 'RSSD0002.TF' );
%         cspice_furnsh( 'earth_latest_high_prec.bpc' );
%
%         %
%         % Count the number of loaded kernel files.
%         %
%         [count] = cspice_ktotal( 'ALL' );
%
%         fprintf( [ 'The total number of kernels after final',            ...
%                    ' cspice_furnsh:  %1d\n' ], count          )
%         fprintf( ' \n' )
%
%         %
%         % Unload the text kernels.
%         %
%         [count] = cspice_ktotal( 'TEXT' );
%
%         fprintf( ' \n' )
%         fprintf( 'Unloading %1d text kernels...\n', count )
%         fprintf( ' \n' )
%
%         while ( count != 0 )
%
%            [file,   filtyp, srcfil,                                      ...
%             handle, found]          = cspice_kdata( 1, 'TEXT' );
%
%            %
%            % If the kernel is found in the pool, unload it.
%            %
%            if ( found )
%
%               cspice_unload( file );
%
%               %
%               % Check if the file has been unloaded.
%               %
%               [filtyp, srcfil, handle, found] = cspice_kinfo( file );
%
%               if ( found )
%
%                  fprintf( '  Error'    )
%
%               else
%
%                   fprintf( '  Success' )
%
%               end
%
%               fprintf( ' unloading %s\n', file )
%
%               %
%               % Something is not working. Inform NAIF.
%               %
%            else
%
%               fprintf( [ ' ERROR: No kernel found but cspice_ktotal',    ...
%                          ' returns  %d\n' ], count                    )
%
%            end
%
%            %
%            % Check if we have more text kernels to unload from
%            % the kernel pool. Note that unloading a text kernel
%            % or meta-kernel implies that the kernel pool is
%            % cleared, and any kernel(s) that were not to be
%            % unloaded are re-loaded. Therefore the `count' value
%            % changes, and the indexing of the files within the
%            % kernel pool too.
%            %
%            [count] = cspice_ktotal( 'TEXT' );
%
%         end
%
%         [count] = cspice_ktotal( 'ALL' );
%
%         fprintf( ' \n' )
%         fprintf( [ 'The total number of kernels after cspice_unload',    ...
%                    ' calls:  %1d\n' ], count                          )
%
%         %
%         % Unload the meta-kernel and retrieve the number of loaded
%         % after the clear.
%         %
%         cspice_unload( 'unload_ex2.tm' );
%
%         [count] = cspice_ktotal( 'ALL' );
%
%         fprintf( ' \n' )
%         fprintf( [ 'The total number of kernels after final',            ...
%                    ' cspice_unload:  %1d\n' ], count          )
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
%      The total number of kernels after final cspice_furnsh:  6
%
%
%      Unloading 3 text kernels...
%
%        Success unloading naif0012.tls
%        Success unloading pck00009.tpc
%        Success unloading RSSD0002.TF
%
%      The total number of kernels after cspice_unload calls:  3
%
%      The total number of kernels after final cspice_unload:  1
%
%
%   3) Redo the previous example, using the cspice_unload capability
%      of unloading kernels by listing them in an array of strings.
%
%      Use the meta-kernel and kernels of Example 2.
%
%
%      Example code begins here.
%
%
%      function unload_ex3()
%
%         %
%         % Load several kernel files.
%         %
%         cspice_furnsh( 'unload_ex2.tm' );
%         cspice_furnsh( 'RSSD0002.TF' );
%         cspice_furnsh( 'earth_latest_high_prec.bpc' );
%
%         %
%         % Count the number of loaded kernel files.
%         %
%         [count] = cspice_ktotal( 'ALL' );
%
%         fprintf( [ 'The total number of kernels after final',            ...
%                    ' cspice_furnsh:  %1d\n' ], count          )
%         fprintf( ' \n' )
%
%         %
%         % Unload the text kernels.
%         %
%         [count] = cspice_ktotal( 'TEXT' );
%
%         fprintf( ' \n' )
%         fprintf( 'Unloading %1d text kernels...\n', count )
%         fprintf( ' \n' )
%
%         %
%         % Create an empty array of strings, to hold the
%         % names of the kernels to unload.
%         %
%         kernels = cell(count,1);
%         for i = 1: count
%
%            [file,   filtyp, srcfil,                                      ...
%             handle, found]          = cspice_kdata( i, 'TEXT' );
%
%            %
%            % If the kernel is found in the pool, add it to the array.
%            %
%            if ( found )
%
%               kernels(i) = file;
%               fprintf( ' %s will be unloaded.\n', file );
%
%            %
%            % Something is not working. Inform NAIF.
%            %
%            else
%
%               fprintf( ' ERROR: No kernel found with index %d\n', i )
%
%            end
%
%         end
%
%         %
%         % Unload the kernels present in the `kernels' variable, and
%         % retrieve the number of remaining loaded kernels.
%         %
%         cspice_unload( kernels );
%
%         [count] = cspice_ktotal( 'ALL' );
%
%         fprintf( ' \n' )
%         fprintf( [ 'The total number of kernels after cspice_unload',    ...
%                    ' call :  %1d\n' ], count                          )
%
%         %
%         % Unload the meta-kernel and retrieve the number of loaded
%         % after the clear.
%         %
%         cspice_unload( 'unload_ex2.tm' );
%
%         [count] = cspice_ktotal( 'ALL' );
%
%         fprintf( ' \n' )
%         fprintf( [ 'The total number of kernels after final',            ...
%                    ' cspice_unload:  %1d\n' ], count          )
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
%      The total number of kernels after final cspice_furnsh:  6
%
%
%      Unloading 3 text kernels...
%
%       naif0012.tls will be unloaded.
%       pck00009.tpc will be unloaded.
%       RSSD0002.TF will be unloaded.
%
%      The total number of kernels after cspice_unload call :  3
%
%      The total number of kernels after final cspice_unload:  1
%
%
%-Particulars
%
%   The call
%
%      cspice_unload ( file );
%
%   has the effect of "erasing" the last previous call:
%
%      cspice_furnsh ( file );
%
%   This interface allows you to unload binary and text kernels.
%   Moreover, if you used a meta-text kernel to set up your
%   working environment, you can unload all of the kernels loaded
%   through the meta-kernel by unloading the meta-kernel.
%
%
%   Unloading Text Kernels or Meta-Kernels
%   --------------------------------------
%
%   Part of the action of unloading text (or meta-text kernels) is
%   clearing the kernel pool and re-loading any kernels that were not in
%   the specified set of kernels to unload. Since loading of text
%   kernels is not a very fast process, unloading text kernels takes
%   considerably longer than unloading binary kernels. Moreover, since
%   the kernel pool is cleared, any kernel pool variables you have set
%   from your program by using one of the interfaces cspice_pcpool,
%   cspice_pdpool, cspice_pipool, or cspice_lmpool will be removed from the
%   kernel pool. For this reason, if you plan to use this feature in your
%   program, together with one of the routines specified above, you will need
%   to take special precautions to make sure kernel pool variables required
%   by your program do not inadvertently disappear.
%
%-Exceptions
%
%   1)  If the specified kernel is not on the list of loaded kernels
%       no action is taken.
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
%   1)  See the note regarding the unloading of Text and meta-text
%       Kernels.
%
%-Required_Reading
%
%   MICE.REQ
%   KERNEL.REQ
%   PCK.REQ
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
%   -Mice Version 1.2.0, 24-AUG-2021 (EDW) (JDR)
%
%       Fixed bug in vectorized unloading of kernels.
%
%       Edited the header to comply with NAIF standard. Added
%       example's problem statement and a reference to the required PCK.
%       Extended the example's output and added the solution to the
%       section.
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
%   -Mice Version 1.1.1, 13-FEB-2015 (EDW)
%
%       Edited -I/O section to conform to NAIF standard for Mice
%       documentation.
%
%   -Mice Version 1.1.0, 17-DEC-2008 (EDW)
%
%       Added zzmice_str call on input `file' to convert string cells to
%       character arrays if `file' has type string cells. Properly
%       identified `file' as a vectorizable string/character array.
%
%   -Mice Version 1.0.0, 22-NOV-2005 (EDW)
%
%-Index_Entries
%
%   Unload a SPICE kernel
%
%-&

function cspice_unload(file)

   switch nargin
      case 1

         file = zzmice_str( file );

      otherwise

         error ( 'Usage: cspice_unload(_`file`_)' )

   end

   %
   % Call the MEX library.
   %
   try
      mice('unload_c', file)
   catch spiceerr
      rethrow(spiceerr)
   end

