%-Abstract
%
%   CSPICE_DAFOPW opens a DAF for subsequent write requests.
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
%      fname    the string name of a DAF to open for write access.
%
%               [1,c1] = size(fname); char = class(fname)
%
%                  or
%
%               [1,1] = size(fname); cell = class(fname)
%
%   the call:
%
%      [handle] = cspice_dafopw( fname )
%
%   returns:
%
%      handle   the file handle associated with the file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%               This handle is used to identify the file in subsequent
%               calls to other DAF routines.
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
%   1) Delete the entire comment area of a DAF file. Note that this
%      action should only be performed if fresh new comments are to
%      be placed within the DAF file.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         earthstns_itrf93_201023.bsp
%
%
%      Example code begins here.
%
%
%      function dafopw_ex1()
%
%         %
%         % Local parameters
%         %
%         KERNEL = 'earthstns_itrf93_201023.bsp';
%         BUFFSZ = 10;
%         LINLEN = 1000;
%
%         %
%         % Open a DAF for write. Return a `handle' referring to the
%         % file.
%         %
%         [handle] = cspice_dafopw( KERNEL );
%
%         %
%         % Print the first 10 lines of comments from the DAF file.
%         %
%         fprintf( 'Comment area of input DAF file (max. 10 lines): \n' )
%         fprintf( ['--------------------------------',                    ...
%                   '--------------------------------\n'] )
%
%         [buffer, done] = cspice_dafec( handle, BUFFSZ, LINLEN );
%
%         for i=1:size(buffer,1)
%            fprintf( '%s\n', buffer(i,:) )
%         end
%
%         fprintf( ['--------------------------------',                    ...
%                   '--------------------------------\n'] )
%         fprintf( ' \n' )
%         fprintf( 'Deleting entire comment area...\n' )
%
%         %
%         % Delete all the comments from the DAF file.
%         %
%         cspice_dafdc( handle );
%
%         %
%         % Close the DAF file and re-open it for read
%         % access to work around the cspice_dafec restriction
%         % on comments not to be modified while they are
%         % being extracted.
%         %
%         cspice_dafcls( handle );
%
%         [handle] = cspice_dafopr( KERNEL );
%
%         %
%         % Check if the comments have indeed been deleted.
%         %
%         [buffer, done] = cspice_dafec( handle, BUFFSZ, LINLEN );
%
%         if ( done && numel(buffer) == 0 )
%            fprintf( ' \n' )
%            fprintf( '   Successful operation.\n' )
%         else
%            fprintf( ' \n' )
%            fprintf( '   Operation failed.\n' )
%         end
%
%         %
%         % Safely close the DAF.
%         %
%         cspice_dafcls( handle );
%
%
%      When this program was executed on a Mac/Intel/Octave6.x/64-bit
%      platform, the output was:
%
%
%      Comment area of input DAF file (max. 10 lines):
%      ----------------------------------------------------------------
%
%         SPK for DSN Station Locations
%         =====================================================================
%
%         Original file name:                   earthstns_itrf93_201023.bsp
%         Creation date:                        2020 October 28 12:30
%         Created by:                           Nat Bachman  (NAIF/JPL)
%
%
%         Introduction
%      ----------------------------------------------------------------
%
%      Deleting entire comment area...
%
%         Successful operation.
%
%
%-Particulars
%
%   Most DAFs require only read access. If you do not need to
%   change the contents of a file, you should open it with cspice_dafopr.
%   Use cspice_dafopw when you need to
%
%      -- change (update) one or more summaries, names, or
%         arrays within a file; or
%
%      -- add new arrays to a file.
%
%   Use cspice_dafcls to close files opened by this routine.
%
%-Exceptions
%
%   1)  If the specified file has already been opened, either by the
%       DAF routines or by other code, an error is signaled by a
%       routine in the call tree of this routine. Note that this
%       response is not paralleled by cspice_dafopr, which allows you to open
%       a DAF for reading even if it is already open for reading.
%
%   2)  If the specified file cannot be opened without exceeding the
%       maximum number of files, the error SPICE(DAFFTFULL) is
%       signaled by a routine in the call tree of this routine.
%
%   3)  If the attempt to read the file's file record fails, the error
%       SPICE(FILEREADFAILED) is signaled by a routine in the call
%       tree of this routine.
%
%   4)  If the specified file is not a DAF file, an error is signaled
%       by a routine in the call tree of this routine.
%
%   5)  If no logical units are available, an error is signaled by a
%       routine in the call tree of this routine.
%
%   6)  If the file does not exist, an error is signaled by a routine
%       in the call tree of this routine.
%
%   7)  If an I/O error occurs in the process of opening the file, the
%       error is signaled by a routine in the call tree of this
%       routine.
%
%   8)  If the file name is blank or otherwise inappropriate, an error
%       is signaled by a routine in the call tree of this routine.
%
%   9)  If the file was transferred improperly via FTP, an error is
%       signaled by a routine in the call tree of this routine.
%
%   10) If the file utilizes a non-native binary file format, an error
%       is signaled by a routine in the call tree of this routine.
%
%   11) If the input argument `fname' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   12) If the input argument `fname' is not of the expected type, or
%       it does not have the expected dimensions and size, an error is
%       signaled by the Mice interface.
%
%-Files
%
%   See argument `fname'.
%
%-Restrictions
%
%   1)  Only files of the native binary file format may be opened
%       with this routine.
%
%   2)  Files opened using this routine must be closed with cspice_dafcls.
%
%-Required_Reading
%
%   DAF.REQ
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
%   -Mice Version 1.1.0, 25-NOV-2021 (EDW) (JDR)
%
%       Edited the header to comply with NAIF standard. Added complete
%       code example.
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
%   -Mice Version 1.0.0, 10-JUL-2012 (EDW)
%
%-Index_Entries
%
%   open existing DAF for write
%
%-&

function [handle] = cspice_dafopw( fname )

   switch nargin
      case 1

         fname  = zzmice_str(fname);

      otherwise

         error ( 'Usage: [handle] = cspice_dafopw(`fname`)' )

   end

   %
   % Call the MEX library.
   %
   try
      [handle] = mice( 'dafopw_c', fname );
   catch spiceerr
      rethrow(spiceerr)
   end




