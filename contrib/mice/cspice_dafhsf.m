%-Abstract
%
%   CSPICE_DAFHSF returns the summary format associated with a handle.
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
%      handle   the handle associated with a previously opened DAF file.
%
%               [1,1] = size(handle); int32 = class(handle)
%
%   the call:
%
%      [nd, ni] = cspice_dafhsf( handle )
%
%   returns:
%
%      nd,
%      ni       the numbers of double precision and integer components,
%               respectively, in each array summary in the specified file.
%
%               [1,1] = size(nd); int32 = class(nd)
%               [1,1] = size(ni); int32 = class(ni)
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
%   1) Find the number of d.p. `words' in a DAF having an
%      arbitrary summary format.
%
%      Use the SPK kernel below as input DAF file for the program.
%
%         de421.bsp
%
%
%      Example code begins here.
%
%
%      function dafhsf_ex1()
%
%         %
%         % Local variables.
%         %
%         daf = 'de421.bsp';
%
%         %
%         % Open the `daf' and find the summary format.
%         %
%         [handle] = cspice_dafopr( daf );
%         [nd, ni] = cspice_dafhsf( handle );
%
%         %
%         % Start a forward search and examine each array in
%         % turn.
%         %
%         cspice_dafbfs( handle );
%         [found] = cspice_daffna;
%
%         n       = 0;
%         while found
%
%            %
%            % Obtain the array summary, unpack it, and get
%            % the initial and final array addresses from
%            % the integer descriptor component.
%            %
%            [dc, ic] = cspice_dafgs( nd, ni );
%
%            ia      =  ic( ni - 1);
%            fa      =  ic( ni );
%
%            n       =  fa - ia + 1 + n;
%
%            [found] = cspice_daffna;
%
%         end
%
%         fprintf( 'Number of d.p. words is   %d\n', n )
%
%
%      When this program was executed on a Mac/Intel/Octave5.x/64-bit
%      platform, the output was:
%
%
%      Number of d.p. words is   2098004
%
%
%-Particulars
%
%   The summary format must be known in order to pack or unpack
%   an array summary. See the DAF Required Reading for a discussion
%   of summary formats.
%
%-Exceptions
%
%   1)  If the specified handle does not belong to any file that is
%       currently known to be open, the error SPICE(DAFNOSUCHHANDLE)
%       is signaled by a routine in the call tree of this routine.
%
%   2)  If the input argument `handle' is undefined, an error is
%       signaled by the Matlab error handling system.
%
%   3)  If the input argument `handle' is not of the expected type, or
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
%
%-Version
%
%   -Mice Version 1.0.0, 09-AUG-2021 (JDR)
%
%-Index_Entries
%
%   handle to DAF summary format
%
%-&
function [nd, ni] = cspice_dafhsf( handle )

   switch nargin
      case 1

         handle = zzmice_int(handle);

      otherwise

         error ( 'Usage: [nd, ni] = cspice_dafhsf( handle )' )

   end

   %
   % Call the MEX library.
   %
   try
      [nd, ni] = mice('dafhsf_c', handle);
   catch spiceerr
      rethrow(spiceerr)
   end
