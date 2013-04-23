%-Abstract
%
%   CSPICE_KDATA returns data for the nth kernel among a list of specified 
%   kernel types. 
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
%      which   the scalar integer index of the kernel to fetch 
%              (matching the type specified by kind) from the list 
%              of kernels loaded using cspice_furnsh 
%              but not unloaded using cspice_unload or cleared
%              by cspice_klear 
%
%              The range of 'which' is 1 to 'count', where 'count' is 
%              the number of kernels loaded via cspice_furnsh. Retrieve
%              this value from a cspice_ktotal call.  See the 
%              Examples section for an illustrative code fragment.
%
%     kind     a scalar string list of types of kernels to consider when 
%              fetching kernels from the list of loaded kernels. 'kind'
%              should consist of a list of words of kernels to 
%              examine.  Recognized types are 
% 
%                 SPK
%                 CK 
%                 PCK 
%                 EK 
%                 TEXT
%                 META
%                 ALL 
% 
%              'kind' lacks case sensitivity. The cspice_kdata algorithm
%              ignores words in 'kind' if not one of those listed above.
%
%              See the routine cspice_ktotal for example use 'kind'. 
%
%   the call:
%
%      [file, filtyp, source, handle, found] = cspice_kdata(which, kind)
%
%   returns:
%
%      file     the scalar string name of the file having index 'which' in the 
%               sequence of files of type 'kind' currently loaded via 
%               cspice_furnsh. 'file' returns empty if no loaded kernels 
%               match the specification of 'which' and 'kind'.
%
%      filtyp   the scalar string name of the type of kernel specified 
%               by 'file'. 'filtyp' will be empty if no loaded kernels
%               match the specification of 'which' and 'kind'.
%
%      source   the scalar string name of the source file used to 
%               specify file as one to load.  If file was loaded 
%               directly via a call to cspice_furnsh, 'source' will be empty. 
%               If there is no file matching the specification of 
%               which and kind, 'source' will be empty. 
%       
%      handle   the scalar integer handle attached to file if it is a binary 
%               kernel.  If file is a text kernel or meta-text kernel 
%               handle will be zero.  If there is no file matching 
%               the specification of 'which' and 'kind', 'handle' will be 
%               set to zero.
%
%      found    returns true if a file matching the specification 
%               of 'which' and 'kind' exists. If there is no such file, 
%               'found' will return false. 
%
%-Examples
%
%   Any numerical results shown for this example may differ between
%   platforms as the results depend on the SPICE kernels used as input
%   and the machine specific arithmetic implementation.
%
%      %
%      % Load several kernel files.
%      %
%      cspice_kclear
%      cspice_furnsh( 'standard.tm' )
%   
%      % 
%      % Count the number of loaded kernel files. 
%      %
%      count = cspice_ktotal( 'ALL' );
%   
%      %
%      % Loop over the count, outputting file information as we loop.
%      % The loop tells us all files loaded via cspice_furnsh, their
%      % type, and how they were loaded.
%      %
%      for i = 1:count+1
%   
%         [ file, type, source, handle, found ] = ...
%                                          cspice_kdata( i, 'ALL');
%   
%         if ( found )
%            fprintf( 'Index : %d\n', i     );            
%            fprintf( 'File  : %s\n', file  );
%            fprintf( 'Type  : %s\n', type  );
%            fprintf( 'Source: %s\n\n', source);
%
%         else
%
%            fprintf( 'No kernel found with index: %d\n', i );
%            
%         end
%   
%      end
%
%      %
%      % It's always good form to unload kernels after use,
%      % particularly in Mice due to data persistence.
%      %
%      cspice_kclear
%
%   MATLAB outputs:
%
%      Index : 1
%      File  : standard.tm
%      Type  : META
%      Source: 
%
%      Index : 2
%      File  : /kernels/gen/lsk/naif0008.tls
%      Type  : TEXT
%      Source: standard.tm
%       
%      Index : 3
%      File  : /kernels/gen/spk/de405_2000-2050.bsp
%      Type  : SPK
%      Source: standard.tm
%
%      Index : 4
%      File  : /kernels/gen/pck/pck00008.tpc
%      Type  : TEXT
%      Source: standard.tm
%
%      No kernel found with index: 5
%
%-Particulars
%
%   None.
%
%-Required Reading
%
%   For important details concerning this module's function, please refer to
%   the CSPICE routine kdata_c.
%
%   MICE.REQ
%
%-Version
%
%   -Mice Version 1.0.1, 06-MAY-2009, EDW (JPL)
%
%      Added MICE.REQ reference to the Required Reading section.
%      
%   -Mice Version 1.0.0, 30-MAR-2007, EDW (JPL)
%
%-Index_Entries
% 
%   Retrieve information on loaded SPICE kernels 
% 
%-&

function [ file, filtyp, source, handle, found ] = cspice_kdata( which, kind )

   switch nargin
      case 2

         which = zzmice_int(which);
         kind  = zzmice_str(kind);

      otherwise

         error( [ 'Usage: [ `file`, `filtyp`, `source`, handle ] = ' ...
                                         'cspice_kdata( which, `kind` )']  )

   end

   %
   % Call the MEX library.
   %
   try
      [ file, filtyp, source, handle, found ]  = mice('kdata_c', which, kind);

      %
      % Convert the integer flags to MATLAB logicals for return to
      % the caller.
      %
      found = logical(found);
   catch
      rethrow(lasterror)
   end


